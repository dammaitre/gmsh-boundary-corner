# GMSH FORK / BOUNDARY CORNER

Le but de ce fork est de créer une nouvelle classe BoundaryCornerField. Son fonctionnement sera analogue à la classe BoundaryLayerField de Fields.cpp.

## Motivation

Contexte CFD 2d axisymétrique. GMSH ne sait pas gérer les intersections entre couches limites et axe y^{-} à 90°. Par exemple pour un profil elliptique coupé par l'axe de symétrie. Le but ici est de gérer une interface entre l'axe y=0 et une couche limite.

## Descriptions de la stratégie de mesh

Notations : 
 + **P :** la courbe du profile s = (x ; y(x))
 + **A :** l'axe y=0
 + **Point C = (x_C ; 0) :** Le point criique au bord d'attaque ou au bord de fuite, intersection perpendiculaire entre P et A.
 + **h1, r et n_BL :** Paramètres de maillage de couche limite : hauteur du premier layer, nombre de layer, ratio de croissance géometrique de la hauteur du layer
 + **l_BL :** longueur des cellules de la couche limite

 + **Point S = s(x_S) :** Point de départ du BC. Entre S et C, le BC sera appliqué ; C marque la fin du domaine ; on supposera toujours qu'une couche limite est maillée avant S
 + **delta_1 :** longueur de la dernière cellule du BC. Cette cellule est positionnée au coin formé par P, A au point C.
 + **n_BC :** Nombre de colonnes du BL en comptant la dernière au coin.

À partir du point S, et de manière analogue au maillage de couche limite, on créé n_BL layers x n_BC colonnes de cellules, les layers étant empilés dans la direction normale au profile.

Première grande différence : décroissance géométrique de la longueur des colonnes. La première cellule (au coin S) a une longueur **l_BL** ; la dernière cellule a une longueur **delta_1**. Les cellules du k-ieme layer du BC et du BL ont la même hauteur.

Deuxième grande différence : gestion des interfaces BL // BC et BC // A. Afin d'éviter de créer des interstices aux interfaces, les points de la première colonne du BC situés au dessus du point S sont ceux du BL adjacents. Idem : la dernière colonne de points du BC au niveau de C sont tous sur l'axe y=0 puisque le profil le tape à 90°. Il faut donc découper l'axe en n_BL segments verticaux.

## Analogie avec BL

Lire la documentation du BoundaryLayerField dans local/BL_docs.md


