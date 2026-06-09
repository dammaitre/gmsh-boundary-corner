# BoundaryCorner — guide d'implémentation pour agent Claude

## Contexte

Ajout d'un nouveau champ `BoundaryCorner` au code source de gmsh (C++).
Ce champ gère la transition entre une couche limite structurée et un point
critique où le profil intersecte l'axe de symétrie à 90°, cas typique des
maillages axisymétriques OpenFOAM wedge.

Repo : `https://gitlab.onelab.info/gmsh/gmsh.git`
Branche de travail : `feature/boundary-corner`

---

## Fichiers à modifier ou créer

```
Field/Field.cpp          ← implémentation principale (modifier)
Field/Field.h            ← déclaration de la classe (modifier)
Mesh/meshGModel.cpp      ← branchement dans le pipeline (modifier)
api/gmsh.h               ← exposition API C (modifier)
api/gmsh.cpp             ← exposition API C (modifier)
```

Aucun autre fichier à toucher.

---

## 1. Field/Field.h — déclaration

Ajouter la déclaration de la classe après la déclaration de `FieldBoundaryLayer` :

```cpp
class FieldBoundaryCorner : public Field {
public:
  FieldBoundaryCorner();
  ~FieldBoundaryCorner() override = default;

  // Implémentation de l'interface Field
  double operator()(double x, double y, double z,
                    GEntity *ge = nullptr) override;
  const char *getName() override { return "BoundaryCorner"; }
  std::string getDescription() override;

  // Construction des colonnes de quads dans la GFace
  void buildCornerColumns(GModel *gm);

private:
  // Paramètres exposés (enregistrés via options dans le constructeur)
  std::vector<int> curvesList_;    // CurvesList  : tags des GEdge du profil
  double axisPoint_[2];            // AxisPoint   : (x_a, 0.0)
  double startPoint_[2];           // StartPoint  : (s, y(s))
  double h1_;                      // Size        : hauteur première rangée BL
  double ratio_;                   // Ratio       : raison progression BL
  int    nbLayers_;                // NbLayers    : nombre de rangées BL
  int    nbCornerColumns_;         // NbCornerColumns : n_BC
  double delta1_;                  // Delta1      : longueur dernière cellule
  double omega_;                   // Omega       : facteur hauteur (défaut 1.0)

  // Données internes calculées par buildCornerColumns()
  std::vector<SPoint2> ptsProfil_; // coordonnées des pts sur le profil
  std::vector<SPoint2> ptsAxe_;    // coordonnées des pts sur l'axe
  std::vector<SPoint2> ptsBLTop_;  // coordonnées des pts sur l'offset BL
  double eps_;                     // raison progression longitudinale
  double hTotal_;                  // hauteur totale BL

  // Méthodes privées
  void   computeParameters();
  void   distributeProfilePoints(GEdge *ge);
  SPoint2 normalAtPoint(GEdge *ge, double t); // normale unitaire en paramètre t
  double  arcLengthToParam(GEdge *ge, double s); // abscisse curviligne → t
};
```

---

## 2. Field/Field.cpp — implémentation

### 2a. Constructeur — enregistrement des paramètres

Le pattern est identique à `FieldBoundaryLayer`. Chercher dans `Field.cpp`
la définition du constructeur de `FieldBoundaryLayer` et reproduire le même
mécanisme `options[...] = new FieldOption...` :

```cpp
FieldBoundaryCorner::FieldBoundaryCorner()
  : h1_(1e-3), ratio_(1.15), nbLayers_(30),
    nbCornerColumns_(10), delta1_(-1.0), omega_(1.0)
{
  // delta1 = -1 → sera remplacé par h1 dans computeParameters()

  options["CurvesList"] = new FieldOptionList(
    curvesList_, "Tags des courbes du profil", &update_needed);

  options["AxisPoint"] = new FieldOptionList(
    /* stocker dans un vecteur temporaire converti dans computeParameters() */
    "Point critique sur l'axe [x_a, 0.0]", &update_needed);

  options["StartPoint"] = new FieldOptionList(
    "Point de depart du BC sur le profil [x, y]", &update_needed);

  options["Size"] = new FieldOptionDouble(
    h1_, "Hauteur premiere rangee BL", &update_needed);

  options["Ratio"] = new FieldOptionDouble(
    ratio_, "Raison progression geometrique BL", &update_needed);

  options["NbLayers"] = new FieldOptionInt(
    nbLayers_, "Nombre de rangees BL", &update_needed);

  options["NbCornerColumns"] = new FieldOptionInt(
    nbCornerColumns_, "Nombre de colonnes dans la zone BC", &update_needed);

  options["Delta1"] = new FieldOptionDouble(
    delta1_, "Longueur derniere cellule BC (-1 = egal a Size)", &update_needed);

  options["Omega"] = new FieldOptionDouble(
    omega_, "Facteur d'echelle des hauteurs BC", &update_needed);
}
```

### 2b. computeParameters()

Appelé en début de `buildCornerColumns()` et de `operator()`.

```cpp
void FieldBoundaryCorner::computeParameters()
{
  // delta1 par défaut = h1
  if(delta1_ < 0.0) delta1_ = h1_;

  // hauteur totale BL
  if(std::abs(ratio_ - 1.0) < 1e-10)
    hTotal_ = h1_ * nbLayers_;
  else
    hTotal_ = h1_ * (std::pow(ratio_, nbLayers_) - 1.0) / (ratio_ - 1.0);

  // taille tangentielle courante l : estimée depuis la discrétisation
  // existante de la première courbe de curvesList_
  // (longueur de la GEdge / nombre de nœuds mesh existants)
  // → utiliser ge->mesh_vertices.size() et ge->length()
  GEdge *ge = GModel::current()->getEdgeByTag(curvesList_[0]);
  double l = ge->length() / std::max(1, (int)ge->mesh_vertices.size());

  // raison de progression longitudinale
  // l * eps^(nbCornerColumns_-1) = delta1_
  if(nbCornerColumns_ > 1)
    eps_ = std::pow(delta1_ / l, 1.0 / (nbCornerColumns_ - 1));
  else
    eps_ = 1.0;
}
```

### 2c. buildCornerColumns()

C'est la méthode centrale. Elle crée les `MVertex` et `MQuadrangle`
directement dans la `GFace` du domaine fluide, avant que `meshGFaces()`
ne la traite.

Logique en pseudo-code C++ :

```cpp
void FieldBoundaryCorner::buildCornerColumns(GModel *gm)
{
  computeParameters();

  GEdge *ge = gm->getEdgeByTag(curvesList_[0]);

  // 1. Trouver le paramètre t_start correspondant à startPoint_
  double t_start = arcLengthToParam(ge, /* abscisse de startPoint_ */);
  double t_end   = arcLengthToParam(ge, /* abscisse de axisPoint_  */);

  // 2. Distribuer nbCornerColumns_+1 points sur le profil
  //    selon la progression eps_ entre t_start et t_end
  ptsProfil_.resize(nbCornerColumns_ + 1);
  ptsProfil_[0] = SPoint2(startPoint_[0], startPoint_[1]);
  for(int i = 1; i <= nbCornerColumns_; i++) {
    double L_i = h1_ * std::pow(eps_, i - 1); // longueur du segment i
    // avancer de L_i le long de la courbe depuis la position actuelle
    // utiliser GEdge::reparamOnFace() ou intégration de la métrique
    ptsProfil_[i] = advanceAlongCurve(ge, ptsProfil_[i-1], L_i);
  }
  ptsProfil_[nbCornerColumns_] = SPoint2(axisPoint_[0], axisPoint_[1]);

  // 3. Calculer pts sur l'axe (projection y=0)
  ptsAxe_.resize(nbCornerColumns_ + 1);
  for(int i = 0; i <= nbCornerColumns_; i++)
    ptsAxe_[i] = SPoint2(ptsProfil_[i].x(), 0.0);

  // 4. Calculer pts sur l'offset BL (normale * hTotal_)
  ptsBLTop_.resize(nbCornerColumns_ + 1);
  for(int i = 0; i <= nbCornerColumns_; i++) {
    double t    = arcLengthToParam(ge, ptsProfil_[i].x());
    SPoint2 n   = normalAtPoint(ge, t);
    ptsBLTop_[i] = SPoint2(ptsProfil_[i].x() + n.x() * hTotal_,
                            ptsProfil_[i].y() + n.y() * hTotal_);
  }

  // 5. Identifier la GFace cible
  // La GFace est celle qui contient ge comme GEdge de bord
  GFace *gf = /* gm->getFaceContaining(ge) — voir GModel::getFaceByTag() */ ;

  // 6. Créer les MVertex et MQuadrangle dans gf
  // Pour chaque colonne i = 0..nbCornerColumns_-1 :
  //   Nappe inférieure (profil → axe) : n_BL+1 nœuds en hauteur
  //   Nappe supérieure (profil → offset) : n_BL+1 nœuds en hauteur
  // Utiliser le même pattern que BoundaryLayers.cpp :
  //   MVertex *v = new MVertex(x, y, 0.0, gf);
  //   gf->mesh_vertices.push_back(v);
  //   gf->quadrangles.push_back(new MQuadrangle(v0,v1,v2,v3));

  for(int i = 0; i < nbCornerColumns_; i++) {
    // interpoler les nœuds internes en hauteur (progression ratio_*omega_)
    for(int k = 0; k <= nbLayers_; k++) {
      double hk = (k == 0) ? 0.0
                : h1_ * omega_ * (std::pow(ratio_, k) - 1.0) / (ratio_ - 1.0);

      // nœuds nappe inférieure : interpolation linéaire profil → axe
      double xInf = ptsProfil_[i].x()
                  + (ptsAxe_[i].x() - ptsProfil_[i].x()) * hk / hTotal_;
      double yInf = ptsProfil_[i].y() * (1.0 - hk / hTotal_);
      MVertex *vInf = new MVertex(xInf, yInf, 0.0, gf);
      gf->mesh_vertices.push_back(vInf);

      // nœuds nappe supérieure : interpolation profil → offset BL
      // (offset dans la direction normale)
    }
    // assembler les MQuadrangle depuis les nœuds créés
  }
}
```

### 2d. operator()

Retourne la taille de maille locale pour guider `surf_ext` au voisinage du BC.
Utilisé par le champ `Min` global.

```cpp
double FieldBoundaryCorner::operator()(double x, double y, double z,
                                        GEntity *ge)
{
  computeParameters();

  // Distance au point axisPoint_
  double dx = x - axisPoint_[0];
  double dy = y - axisPoint_[1];
  double dist = std::sqrt(dx*dx + dy*dy);

  // Transition linéaire de delta1_ (près du coin) à lc_far (loin)
  // lc_far n'est pas connu ici → retourner delta1_ dans la zone BC
  // et laisser le Min avec les autres champs gérer la transition
  double zoneSize = hTotal_ * 5.0; // rayon d'influence approximatif
  if(dist < zoneSize)
    return delta1_ + (h1_ - delta1_) * dist / zoneSize;
  return 1e22; // hors zone : pas d'influence
}
```

### 2e. Enregistrement dans le FieldManager

À la fin du bloc d'enregistrement des champs existants dans `Field.cpp`
(chercher `GMSH_ADD_FIELD_TYPE` ou `fieldFactory.insert`) :

```cpp
fieldFactory["BoundaryCorner"] =
  []() -> Field * { return new FieldBoundaryCorner(); };
```

---

## 3. Mesh/meshGModel.cpp — branchement dans le pipeline

Chercher dans `meshGModel.cpp` l'appel à `meshGFaces(gm, ...)`.
Insérer **avant** cet appel :

```cpp
// Injection des colonnes BoundaryCorner avant le maillage des surfaces
if(CTX::instance()->mesh.boundaryCornerField > 0) {
  FieldManager *fields = gm->getFields();
  Field *f = fields->get(CTX::instance()->mesh.boundaryCornerField);
  if(f) {
    // Itérer sur tous les champs BoundaryCorner actifs
    // (il peut y en avoir plusieurs : un par point critique)
    for(auto &[tag, field] : *fields) {
      FieldBoundaryCorner *bc =
        dynamic_cast<FieldBoundaryCorner *>(field);
      if(bc) bc->buildCornerColumns(gm);
    }
  }
}

// Appel existant
meshGFaces(gm, ...);
```

Ajouter l'option dans `CTX` (`Common/Context.h`, struct `contextMeshOptions`) :

```cpp
int boundaryCornerField; // tag du champ BoundaryCorner actif (0 = désactivé)
```

Et l'initialiser à 0 dans `Common/Context.cpp`.

---

## 4. Exposition API Python

Dans `api/gmsh.h` et `api/gmsh.cpp`, aucune fonction nouvelle n'est nécessaire :
`BoundaryCorner` est un `Field` standard, accessible via l'API existante :

```python
gmsh.model.mesh.field.add("BoundaryCorner")        # déjà fonctionnel
gmsh.model.mesh.field.setNumber(tag, "Size", h1)   # déjà fonctionnel
```

Seule l'option globale nécessite un ajout dans l'API :

```python
gmsh.option.setNumber("Mesh.BoundaryCornerField", tag)
```

Ce qui correspond dans `api/gmsh.cpp` à exposer `CTX::instance()->mesh.boundaryCornerField`
via le mécanisme existant `setNumber`/`getNumber` — chercher comment
`Mesh.BoundaryLayerField` est exposé et reproduire exactement le même pattern.

---

## Ordre d'implémentation recommandé

- [x] 1. `Field.h` — déclaration de la classe  
- [x] 2. `Field.cpp` — constructeur + `computeParameters()` + `operator()` +
         `buildCornerColumns()` + `arcLengthToParam()` + `normalAtPoint()` +
         enregistrement dans `map_type_name`  
- [x] 3. `Common/Context.h` + `Context.cpp` — ajout de `boundaryCornerField`  
- [x] 4. `Mesh/Generator.cpp` — branchement POST-maillage (voir note ci-dessous)  
- [x] 5. Build complet + test Python : 20 quads + 1983 triangles — pipeline OK  
- [ ] 6. Validation qualité mesh — raccordement frontière quads/triangles

## Progression — état au 2026-06-09

### ✅ Field.h (terminé)

Stub cassé remplacé par la déclaration complète de `BoundaryCornerField`.
Forward declarations ajoutées : `GModel`, `GEdge`, `GFace`, `SPoint2`.

Membres privés — options exposées via `FieldOption*` :
- `curvesList_` (`std::list<int>`) — `CurvesList`
- `axisPointList_`, `startPointList_` (`std::list<double>`) — `AxisPoint`, `StartPoint`
- `h1_`, `ratio_`, `delta1_`, `omega_` (`double`) — `Size`, `Ratio`, `Delta1`, `Omega`
- `nbLayers_`, `nbCornerColumns_` (`int`) — `NbLayers`, `NbCornerColumns`

Membres privés — données calculées : `axisPoint_[2]`, `startPoint_[2]`, `eps_`, `hTotal_`.

Méthodes publiques : `BoundaryCornerField()`, `operator()`, `getName()`,
`getDescription()`, `buildCornerColumns(GModel*)`.

Méthodes privées : `computeParameters()`, `arcLengthToParam(GEdge*,x,y)`,
`normalAtPoint(GEdge*,t)`.

### ✅ Field.cpp (terminé)

Includes : `SPoint2.h`, `GEdge.h`, `GFace.h`, `MVertex.h`, `MTriangle.h`,
`MQuadrangle.h`.

Toutes les méthodes sont implémentées et le champ est enregistré :
```cpp
map_type_name["BoundaryCorner"] = new FieldFactoryT<BoundaryCornerField>();
```

`Field.cpp` compile sans erreur C++.

#### buildCornerColumns — approche post-traitement (décision architecturale)

La spec initiale proposait d'injecter les quads **avant** `meshGFaces()`. Cette
approche ne fonctionne pas : `meshGenerator()` appelle `gf->GFace::deleteMesh()`
à la ligne 1690 de `meshGFace.cpp`, ce qui efface tous les éléments pré-injectés.

**Solution retenue** : appel **après** le maillage 2D complet, dans `Generator.cpp`.
`buildCornerColumns` fonctionne en post-processeur :
1. Calcule la grille structurée `grid[i][k]` (colonnes × hauteurs)
2. Supprime les triangles dont le centroïde est dans la zone corner
   (x ∈ [xMin, xMax], y ∈ [0, hTotal_ × 1.05])
3. Purge les sommets orphelins
4. Injecte les `N × nbLayers_` `MQuadrangle`

Limitation actuelle : la frontière entre la zone quads et la zone triangles
n'est pas raccordée (gap) — les sommets des quads ne coïncident pas avec les
nœuds des triangles voisins. À résoudre en étape 6.

#### mmg3d désactivé

Le linker échouait sur des définitions multiples dans `contrib/mmg3d` (bug
préexistant de la version 4.0 bundlée, sans lien avec nos modifications).
Désactivé avec `-DENABLE_MMG3D=OFF` dans CMake. Sans impact sur notre feature
(mmg3d = raffinement volumique 3D anisotrope).

Build : `cd build && cmake .. -DENABLE_MMG3D=OFF && make -j$(nproc)`

### ✅ Context.h / Context.cpp (terminé)

`mesh.boundaryCornerField` ajouté comme `int` dans `contextMeshOptions`,
initialisé à 0.

### ✅ Options.cpp / Options.h / DefaultOptions.h (terminé)

`Mesh.BoundaryCornerField` enregistré dans le système d'options gmsh via
`opt_mesh_boundary_corner_field()`. Accessible depuis Python :
```python
gmsh.option.setNumber("Mesh.BoundaryCornerField", tag)
```

### ✅ Generator.cpp (terminé)

Hook inséré **après** la boucle de maillage 2D, avant
`Msg::SetNumThreads(prevNumThreads)`. Itère sur tous les champs actifs,
appelle `bc->buildCornerColumns(m)` pour chaque `BoundaryCornerField` trouvé.

### ✅ Test Python (terminé — `test_bc.py`)

Carré unitaire, profil sur y=0, coin en (1,0) :
- 20 quads injectés (4 colonnes × 5 couches)
- 1983 triangles dans le reste de la surface
- Pipeline de bout en bout sans crash

Note : `printf("coucou\n")` dans la sortie — debug préexistant dans
`Mesh/automaticMeshSizeField.cpp:138`, sans lien avec notre code.

### ⏳ Étape 6 — Raccordement frontière quads/triangles

Problème : les sommets des quads injectés ne coïncident pas avec les nœuds des
triangles voisins → gap dans le maillage, invalide pour CFD.

Approches envisagées :
- **Snap** : après injection des quads, déplacer les nœuds triangles du bord
  vers les sommets de quads les plus proches
- **Retriangulation locale** : supprimer une bande plus large de triangles et
  retrianguler en fixant les nœuds des quads comme contraintes
- **Contraintes Delaunay** : passer les arêtes du bord de la grille quad comme
  `embedded edges` dans la GFace avant maillage, forçant le Delaunay à les
  respecter — nécessite de revenir à une approche pré-maillage partielle

---

## Fonctions utilitaires gmsh à réutiliser

Ces fonctions existent déjà dans le code source et sont utilisées par
`BoundaryLayers.cpp` et `Field.cpp` — ne pas les réécrire :

| Besoin | Fonction existante | Fichier |
|---|---|---|
| Longueur d'une GEdge | `ge->length()` | `Geo/GEdge.h` |
| Normale à une GEdge en t | `ge->firstDer(t)` → rotation 90° | `Geo/GEdge.h` |
| Paramètre depuis point | `ge->closestPoint(p, t)` | `Geo/GEdge.h` |
| Créer un nœud dans GFace | `new MVertex(x,y,z,gf)` | `Mesh/MVertex.h` |
| Créer un quad | `new MQuadrangle(v0,v1,v2,v3)` | `Mesh/MElement.h` |
| Accéder aux champs | `gm->getFields()->get(tag)` | `Field/Field.h` |

---

## Test minimal après implémentation

```python
import gmsh

gmsh.initialize()
gmsh.model.add("test_bc")

# Carré simple avec un profil sur le bord bas
p1 = gmsh.model.geo.addPoint(0, 0, 0, 0.1)
p2 = gmsh.model.geo.addPoint(1, 0, 0, 0.1)
p3 = gmsh.model.geo.addPoint(1, 1, 0, 0.1)
p4 = gmsh.model.geo.addPoint(0, 1, 0, 0.1)
l1 = gmsh.model.geo.addLine(p1, p2)   # profil (sur l'axe y=0)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
sf = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize()

# Champ BoundaryCorner sur le coin (1, 0)
bc = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc, "CurvesList",      [l1])
gmsh.model.mesh.field.setNumbers(bc, "AxisPoint",       [1.0, 0.0])
gmsh.model.mesh.field.setNumbers(bc, "StartPoint",      [0.7, 0.0])
gmsh.model.mesh.field.setNumber (bc, "Size",            0.01)
gmsh.model.mesh.field.setNumber (bc, "Ratio",           1.15)
gmsh.model.mesh.field.setNumber (bc, "NbLayers",        5)
gmsh.model.mesh.field.setNumber (bc, "NbCornerColumns", 4)
gmsh.model.mesh.field.setAsBackgroundMesh(bc)

gmsh.option.setNumber("Mesh.BoundaryCornerField", 1)

gmsh.model.mesh.generate(2)

# Vérification minimale
elem_types, elem_tags, _ = gmsh.model.mesh.getElements(2)
for t, tags in zip(elem_types, elem_tags):
    quals = gmsh.model.mesh.getElementQualities(tags, "minSICN")
    print(f"type {t} — n={len(tags)} — qualité min={min(quals):.3f}")
# Attendu : qualité min > 0.3 sur tous les types

gmsh.write("test_bc.msh")
gmsh.finalize()
```

---

## À ne pas faire

- Ne pas modifier `Parser/Gmsh.y` — l'API Python suffit, la syntaxe `.geo`
  peut être ajoutée plus tard si nécessaire
- Ne pas toucher à `FieldBoundaryLayer` — `BoundaryCorner` est une classe
  indépendante qui réutilise ses utilitaires mais ne l'hérite pas
- Ne pas créer de nouvelles `GEdge` ou `GFace` géométriques — injecter
  uniquement des `MVertex` et `MQuadrangle` dans la `GFace` existante
- Ne pas appeler `synchronize()` dans `buildCornerColumns()` — la géométrie
  est déjà synchronisée à ce stade du pipeline