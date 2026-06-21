// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef BOUNDARYTOOLBOX_H
#define BOUNDARYTOOLBOX_H

#include "Field.h"

class GModel;
class GEdge;
class GFace;
class SPoint2;
class MLine;
class MQuadrangle;
class MVertex;

class BoundaryCornerField : public Field {
public:
  BoundaryCornerField();
  ~BoundaryCornerField() override = default;

  double operator()(double x, double y, double z,
                    GEntity *ge = nullptr) override;
  const char *getName() override { return "BoundaryCorner"; }
  std::string getDescription() override;

  void buildCornerColumns(GModel *gm);  // legacy post-mesh entry (kept for reference)

  // Pre-mesh: compute BC geometry for GFace gf.
  // Fills bcQuads (structured quads), verts (all interior+outer BC nodes not
  // on boundary curves), and outerLines (MLine segments along the outer BC
  // boundary, i.e. k=nbLayers_ row).  Does NOT touch gf->triangles.
  // Returns true if this field applies to gf.
  bool buildForFace(GFace *gf,
                    const std::vector<MQuadrangle *> &blQuads,
                    const std::set<MVertex *> &blVerts,
                    std::vector<MQuadrangle *> &bcQuads,
                    std::set<MVertex *> &verts,
                    std::vector<MLine *> &outerLines,
                    std::map<MVertex *, std::vector<MVertex *>> &junctionMap,
                    GEdge *&axisEdgeOut,
                    std::vector<MVertex *> &axisColVertsOut);

private:
  // Options exposées via FieldOption*
  std::list<int>    curvesList_;       // CurvesList      : tags des GEdge du profil
  std::list<double> axisPointList_;    // AxisPoint       : [x_a, 0.0]
  double h1_;                          // Size            : hauteur 1ère rangée BL
  double ratio_;                       // Ratio           : raison progression BL
  int    nbLayers_;                    // NbLayers        : nombre de rangées BL
  int    nbLengthControl_;             // NbLengthControl : geometric-length columns (outer part of corner zone)
  int    nbHeightControl_;             // NbHeightControl : height-blend columns (inner part, at the axis)
  double w0max_;                       // MaxColumnWidth  : arc max colonne (transition width)
  double lBL_;                         // ColWidth : arc dernière colonne à l'AxisPoint ; -1 = Size
  double omega_;                       // Omega           : facteur hauteur (défaut 1.0)
  int    skipAxisColumn_;              // SkipAxisColumn  : 1 = omit the y=0 axis column (for 3D revolve)

  // Données calculées par computeParameters() / buildCornerColumns()
  double axisPoint_[2];
  double lBLeff_;   // lBL_ résolu (> 0 : valeur utilisateur, sinon h1_)
  double hTotal_;

  // Méthodes privées
  void    computeParameters();
  double  arcLengthToParam(GEdge *ge, double x, double y);
  SPoint2 normalAtPoint(GEdge *ge, double t);
};

class BoundaryDoubleCornerField : public Field {
public:
  BoundaryDoubleCornerField();
  ~BoundaryDoubleCornerField() override = default;

  double operator()(double x, double y, double z,
                    GEntity *ge = nullptr) override;
  const char *getName() override { return "BoundaryDoubleCorner"; }
  std::string getDescription() override;

  // Pre-mesh: compute BC geometry for a profile that meets the axis at both ends.
  // Fills bcQuads, verts, outerLines.  Appends up to two {axisEdge, axisColVerts}
  // pairs to axisReclassifyOut (nose end first, then tail end).
  // Returns true if this field applies to gf.
  bool buildForFace(GFace *gf,
                    const std::vector<MQuadrangle *> &blQuads,
                    const std::set<MVertex *> &blVerts,
                    std::vector<MQuadrangle *> &bcQuads,
                    std::set<MVertex *> &verts,
                    std::vector<MLine *> &outerLines,
                    std::map<MVertex *, std::vector<MVertex *>> &junctionMap,
                    std::vector<std::pair<GEdge *,
                      std::vector<MVertex *>>> &axisReclassifyOut);

private:
  std::list<int>    curvesList_;
  std::list<double> nosePointList_;   // NosePoint [x_nose, 0.0]
  std::list<double> tailPointList_;   // TailPoint [x_tail, 0.0]
  double h1_, ratio_, w0max_, lBL_, omega_;
  int    nbLayers_, nbLengthControl_, nbHeightControl_, skipAxisColumn_;

  double nosePoint_[2], tailPoint_[2];
  double lBLeff_, hTotal_;

  void    computeParameters();
  double  arcLengthToParam(GEdge *ge, double x, double y);
  SPoint2 normalAtPoint(GEdge *ge, double t);
  void    subdivideAxisEdge(GEdge *axisEdge, MVertex *axisVert,
                            const std::vector<MVertex *> &axisColVerts,
                            GFace *gf);
};

#endif
