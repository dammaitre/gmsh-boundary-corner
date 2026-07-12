// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include <cstdlib>
#include <limits>
#include <list>
#include <cmath>
#include <algorithm>
#include "GmshConfig.h"
#include "boundarytoolbox.h"
#include "GModel.h"
#include "GmshMessage.h"
#include "Numeric.h"
#include "SPoint2.h"
#include "GEdge.h"
#include "GFace.h"
#include "MVertex.h"
#include "MTriangle.h"
#include "MQuadrangle.h"
#include "MLine.h"
#include "meshGFaceDelaunayInsertion.h"

// ---------------------------------------------------------------------------
// BoundaryCornerField
// ---------------------------------------------------------------------------

std::string BoundaryCornerField::getDescription()
{
  return "Structured quad columns at the corner where a boundary-layer profile "
         "meets the symmetry axis at 90 degrees (typical OpenFOAM wedge mesh).";
}

BoundaryCornerField::BoundaryCornerField()
  : h1_(1e-3), ratio_(1.15), nbLayers_(30),
    nbLengthControl_(5), nbHeightControl_(5), w0max_(0.1), lBL_(-1.0), omega_(1.0),
    skipAxisColumn_(0), farSize_(-1.0), farRatio_(1.2), lBLeff_(1e-3), hTotal_(0.0)
{
  axisPoint_[0] = axisPoint_[1] = 0.0;

  options["CurvesList"] = new FieldOptionList(
    curvesList_, "Tags of the profile curves", &update_needed);
  options["AxisPoint"] = new FieldOptionListDouble(
    axisPointList_, "Critical point on the axis [x_a, 0.0]", &update_needed);
  options["Size"] = new FieldOptionDouble(
    h1_, "Height of the first BL row", &update_needed);
  options["Ratio"] = new FieldOptionDouble(
    ratio_, "Geometric growth ratio between successive BL rows", &update_needed);
  options["NbLayers"] = new FieldOptionInt(
    nbLayers_, "Number of BL rows", &update_needed);
  options["NbLengthControl"] = new FieldOptionInt(
    nbLengthControl_,
    "Number of columns, outer part of the corner zone, over which the column "
    "arc-length transitions geometrically from MaxColumnWidth down to ColWidth.  "
    "Row heights in this region stay BL-geometric (no blending).",
    &update_needed);
  options["NbHeightControl"] = new FieldOptionInt(
    nbHeightControl_,
    "Number of columns, inner part of the corner zone (adjacent to AxisPoint), "
    "over which the column arc-length is held constant at ColWidth.  Row "
    "heights stay BL-geometric all the way to AxisPoint.",
    &update_needed);
  options["MaxColumnWidth"] = new FieldOptionDouble(
    w0max_,
    "Column arc-length in the constant-width section (farthest from AxisPoint) "
    "and maximum arc-length of the compressed corner columns.",
    &update_needed);
  options["ColWidth"] = new FieldOptionDouble(
    lBL_,
    "Arc-length of the innermost BC column at AxisPoint (corner cell).  "
    "-1 = use the height of the last (outermost) BL row",
    &update_needed);
  options["Omega"] = new FieldOptionDouble(
    omega_, "Height scale factor for BC rows (default 1.0)", &update_needed);
  options["SkipAxisColumn"] = new FieldOptionInt(
    skipAxisColumn_,
    "Set to 1 to omit the structured quad column exactly on y=0 (the AxisPoint "
    "column).  Required when the 2D mesh will be revolved for a 3D wedge (e.g. "
    "OpenFOAM axisymmetric); the axis column creates zero-volume cells under "
    "rotation.  Default 0 (column included, correct for 2D axisymmetric).",
    &update_needed);
  options["FarSize"] = new FieldOptionDouble(
    farSize_,
    "Target far-field size that the size field grows toward geometrically "
    "(zone C, beyond MaxColumnWidth) at FarRatio, smoothing the transition to "
    "the background mesh size instead of jumping straight to it.  "
    "<= 0 disables zone C (legacy behaviour: unconstrained beyond zone B).",
    &update_needed);
  options["FarRatio"] = new FieldOptionDouble(
    farRatio_, "Geometric growth ratio used in zone C (FarSize transition)",
    &update_needed);
}

void BoundaryCornerField::computeParameters()
{
  double hLastRow = h1_ * std::pow(ratio_, nbLayers_ - 1);
  lBLeff_ = (lBL_ > 0.0) ? lBL_ : hLastRow;
  if(lBLeff_ < 1e-100) lBLeff_ = 1e-100;

  {
    auto it = axisPointList_.begin();
    axisPoint_[0] = (it != axisPointList_.end()) ? *it++ : 0.0;
    axisPoint_[1] = (it != axisPointList_.end()) ? *it   : 0.0;
  }

  if(std::abs(ratio_ - 1.0) < 1e-10)
    hTotal_ = h1_ * nbLayers_;
  else
    hTotal_ = h1_ * (std::pow(ratio_, nbLayers_) - 1.0) / (ratio_ - 1.0);
}

// Returns parameter t on ge closest to (x, y)
double BoundaryCornerField::arcLengthToParam(GEdge *ge, double x, double y)
{
  double t = 0.0;
  const SPoint3 p(x, y, 0.0);  // const forces the virtual GPoint overload (double &param)
  ge->closestPoint(p, t);       // t is now an output: parameter of closest point
  return t;
}

// Unit inward normal at parameter t (rotate tangent 90° CCW)
SPoint2 BoundaryCornerField::normalAtPoint(GEdge *ge, double t)
{
  SVector3 d = ge->firstDer(t);
  double len = d.norm();
  if(len < 1e-14) return SPoint2(0.0, 1.0);
  return SPoint2(-d.y() / len, d.x() / len);
}

bool BoundaryCornerField::buildForFace(
  GFace *gf,
  const std::vector<MQuadrangle *> &blQuads,
  const std::set<MVertex *> &blVerts,
  std::vector<MQuadrangle *> &bcQuads,
  std::set<MVertex *> &verts,
  std::vector<MLine *> &outerLines,
  std::map<MVertex *, std::vector<MVertex *>> &junctionMap,
  GEdge *&axisEdgeOut,
  std::vector<MVertex *> &axisColVertsOut)
{
  axisEdgeOut = nullptr;
  computeParameters();
  if(curvesList_.empty() || nbLayers_ < 1) return false;

  // Resolve GEdge and check it belongs to gf
  GEdge *ge = nullptr;
  for(int tag : curvesList_) {
    GEdge *e = gf->model()->getEdgeByTag(tag);
    if(!e) continue;
    for(auto *fe : gf->edges()) if(fe == e) { ge = e; break; }
    if(ge) break;
  }
  if(!ge) return false;

  std::vector<MVertex *> baseVerts;
  if(ge->getBeginVertex() && !ge->getBeginVertex()->mesh_vertices.empty())
    baseVerts.push_back(ge->getBeginVertex()->mesh_vertices[0]);
  for(auto *v : ge->mesh_vertices) baseVerts.push_back(v);
  if(ge->getEndVertex() && !ge->getEndVertex()->mesh_vertices.empty())
    baseVerts.push_back(ge->getEndVertex()->mesh_vertices[0]);
  if((int)baseVerts.size() < 2) {
    Msg::Error("BoundaryCorner: GEdge %d has no 1D mesh yet", ge->tag());
    return false;
  }

  // Ensure baseVerts runs from profile-start → axisPoint; reverse if necessary.
  auto sqDist = [](MVertex *v, double x, double y) {
    double dx = v->x() - x, dy = v->y() - y;
    return dx*dx + dy*dy;
  };
  MVertex *vFront = baseVerts.front(), *vBack = baseVerts.back();
  double dFrontA = sqDist(vFront, axisPoint_[0], axisPoint_[1]);
  double dBackA  = sqDist(vBack,  axisPoint_[0], axisPoint_[1]);
  if(dFrontA < dBackA)  // front is closer to axis → reverse so axis end is at back
    std::reverse(baseVerts.begin(), baseVerts.end());

  // Detect how many leading arc_bc segments are consumed by BL.
  // When BL builds a fan at p_start, it may use arc_bc vertices (baseVerts[i])
  // as outer-row (j=2 or j=3) vertices in BL quads.  Such segments are XOR-
  // cancelled by BL and BC must skip them.
  int bcStart = 0;
  {
    std::set<MVertex*> blOuterVerts;
    for(auto *q : blQuads)
      for(int j = 2; j <= 3; j++) blOuterVerts.insert(q->getVertex(j));
    // Walk from baseVerts[1] forward; find the LAST baseVerts[i] that is in blOuterVerts.
    // bcStart = that index (all segments 0..i-1 were consumed).
    for(int i = 1; i < (int)baseVerts.size(); i++) {
      if(blOuterVerts.count(baseVerts[i])) {
        bcStart = i;
      } else break;  // stop at first non-consumed vertex
    }
  }

  // Find the GVertex and GEdge at the axis end (y=0) so the last column can be
  // included as quads.  The axis outer vertices land on y=0; we classify them on
  // the axis GEdge and subdivide its 1D mesh so the XOR bedges cancel cleanly.
  GVertex *axisGV = nullptr;
  for(auto *gv : {ge->getBeginVertex(), ge->getEndVertex()}) {
    if(!gv || gv->mesh_vertices.empty()) continue;
    if(gv->mesh_vertices[0] == baseVerts.back()) { axisGV = gv; break; }
  }
  GEdge *axisEdge = nullptr;
  if(axisGV) {
    for(auto *adj : axisGV->edges()) {
      if(adj == ge) continue;
      GVertex *other = (adj->getBeginVertex() == axisGV) ? adj->getEndVertex()
                                                          : adj->getBeginVertex();
      if(other && std::abs(other->y() - axisPoint_[1]) < 1e-10) { axisEdge = adj; break; }
    }
  }
  // With a known axis edge we include the last column (axis column) as quads.
  // Without it fall back to the old behaviour: skip the last column and let
  // Delaunay fill the small gap with a triangle.
  int N = (axisEdge && !skipAxisColumn_) ? (int)baseVerts.size() - 1 - bcStart
                                        : (int)baseVerts.size() - 2 - bcStart;
  if(N < 1) {
    Msg::Warning("BoundaryCorner: no arc_bc segments remain after BL fan for GEdge %d", ge->tag());
    return false;
  }

  // At the BL/BC junction (baseVerts[bcStart]), the BL fan's last quad has the
  // junction vertex at j=2 and the BL outermost vertex (c1_top) at j=3.
  // Reusing c1_top as grid[0][nbLayers_] closes the gap between the BL outer row
  // and the BC outer row.  Without this the Delaunay fills a ~h_total-wide wedge
  // between the two independently-placed outer vertices with triangles that reach
  // down to the profile.
  MVertex *blOuterAtJunction = nullptr;
  {
    MVertex *jv = baseVerts[bcStart];
    for(auto *q : blQuads)
      if(q->getVertex(2) == jv) { blOuterAtJunction = q->getVertex(3); break; }
  }
  const bool useStitch = (blOuterAtJunction != nullptr) &&
                         (bcStart + 1 < (int)baseVerts.size());

  // Determine whether the normals computed via normalAtPoint point outward.
  // At C=(xC,0) the body meets the axis at 90°; the outward normal must point in
  // the sign(xC) * x̂ direction.  If normalAtPoint gives the opposite sign at C,
  // all normals for this field need to be negated.
  bool flipNormal = false;
  {
    double xC = axisPoint_[0];
    if(std::abs(xC) > 1e-10) {
      MVertex *axisVert = baseVerts.back();
      double tC = arcLengthToParam(ge, axisVert->x(), axisVert->y());
      SPoint2 nC = normalAtPoint(ge, tC);
      flipNormal = (nC.x() * xC < -1e-10);
    }
  }

  // Build grid[i][k]: i indexes columns starting from baseVerts[bcStart].
  // For the last column (i==N, axis column) outer vertices lie on y=0 and are
  // classified on axisEdge so they participate in the edge's 1D mesh.
  // For the first column (i==0, non-axis end), reuse outer vertices from a
  // previously processed BC field if that field already built the same column
  // (junction sharing: two BC arcs meeting at a common GVertex).
  std::vector<std::vector<MVertex *>> grid(N + 1,
    std::vector<MVertex *>(nbLayers_ + 1, nullptr));

  for(int i = 0; i <= N; i++) {
    MVertex *bv = baseVerts[bcStart + i];
    grid[i][0] = bv;
    double ti = arcLengthToParam(ge, bv->x(), bv->y());
    SPoint2 ni = normalAtPoint(ge, ti);
    if(flipNormal) ni = SPoint2(-ni.x(), -ni.y());
    double bx = bv->x(), by = bv->y();

    // Axis column: regardless of curve discretisation the normal must be exactly
    // along x and the base y must be exactly 0 so outer vertices land on y=0.
    if(i == N && axisEdge) {
      double signX = (axisPoint_[0] >= 0.0) ? 1.0 : -1.0;
      ni = SPoint2(signX, 0.0);
      by = axisPoint_[1];
    }
    GEntity *outerEnt = gf;

    // At the non-axis end (i==0), reuse outer vertices if another BC field
    // already owns this junction base vertex — avoids duplicate walls in XOR.
    if(i == 0) {
      auto jit = junctionMap.find(bv);
      if(jit != junctionMap.end() &&
         (int)jit->second.size() == nbLayers_) {
        for(int k = 1; k <= nbLayers_; k++)
          grid[0][k] = jit->second[k - 1];
        continue;
      }
    }

    for(int k = 1; k <= nbLayers_; k++) {
      if(useStitch && i == 0 && k < nbLayers_)
        continue;  // intermediate col-0 vertices unused in stitch path
      double hk = (std::abs(ratio_ - 1.0) < 1e-10)
                      ? h1_ * omega_ * k
                      : h1_ * omega_ * (std::pow(ratio_, k) - 1.0) / (ratio_ - 1.0);
      grid[i][k] = new MVertex(bx + ni.x() * hk, by + ni.y() * hk, 0.0, outerEnt);
    }

    // Register the non-axis end column so a subsequent BC field can share it.
    if(i == 0) {
      std::vector<MVertex *> col(nbLayers_);
      for(int k = 1; k <= nbLayers_; k++) col[k - 1] = grid[0][k];
      junctionMap[bv] = std::move(col);
    }
  }

  // Subdivide the axis GEdge's 1D mesh to include grid[N][1..nbLayers_].
  // grid[N][k] lie on y=0 from axisVert to grid[N][nbLayers_]; they may span
  // several existing axis MLines.  We remove every axis MLine that falls in the
  // range [xNose, xOuter] plus the one that straddles xOuter, then insert the
  // new chain so the XOR bedges from the axis edge and the last-column right
  // side cancel cleanly.
  if(axisEdge) {
    MVertex *axisVert = baseVerts.back();
    double xNose  = axisVert->x();
    double xOuter = grid[N][nbLayers_]->x();
    // xOuter > xNose for a downstream (nose) axis, xOuter < xNose for an
    // upstream (tail) axis — use min/max so both cases find the right lines.
    double xMin = std::min(xNose, xOuter);
    double xMax = std::max(xNose, xOuter);

    MVertex *vReconnect = nullptr;
    std::vector<int> toRemove;
    for(int li = 0; li < (int)axisEdge->lines.size(); li++) {
      MLine *ml   = axisEdge->lines[li];
      MVertex *va = ml->getVertex(0), *vb = ml->getVertex(1);
      bool aIn = (va->x() >= xMin - 1e-14 && va->x() <= xMax + 1e-14);
      bool bIn = (vb->x() >= xMin - 1e-14 && vb->x() <= xMax + 1e-14);
      if(aIn && bIn) {
        toRemove.push_back(li);
      } else if(aIn != bIn) {
        // Straddle: one vertex inside range, one outside (beyond the outer column)
        toRemove.push_back(li);
        vReconnect = aIn ? vb : va;
      }
    }

    if(vReconnect && !toRemove.empty()) {
      int insertPos = toRemove.front();
      for(int i = (int)toRemove.size() - 1; i >= 0; i--) {
        delete axisEdge->lines[toRemove[i]];
        axisEdge->lines.erase(axisEdge->lines.begin() + toRemove[i]);
      }

      // Insert chain running toward axisVert (matching the axis edge direction):
      // vReconnect → grid[N][nbLayers_] → ... → grid[N][1] → axisVert
      axisEdge->lines.insert(axisEdge->lines.begin() + insertPos,
                             new MLine(vReconnect, grid[N][nbLayers_]));
      for(int k = nbLayers_ - 1; k >= 1; k--)
        axisEdge->lines.insert(axisEdge->lines.begin() + ++insertPos,
                               new MLine(grid[N][k + 1], grid[N][k]));
      axisEdge->lines.insert(axisEdge->lines.begin() + ++insertPos,
                             new MLine(grid[N][1], axisVert));

      // Strip and delete the old interior vertices whose MLines we removed.
      // grid[N][k] are classified on gf (not axisEdge) so they are NOT added here.
      auto &mv = axisEdge->mesh_vertices;
      auto rmBegin = std::remove_if(mv.begin(), mv.end(),
        [xMin, xMax](MVertex *v) {
          return v->x() >= xMin - 1e-14 && v->x() <= xMax + 1e-14;
        });
      for(auto vit = rmBegin; vit != mv.end(); ++vit) delete *vit;
      mv.erase(rmBegin, mv.end());
    } else {
      Msg::Warning("BoundaryCorner: failed to subdivide axis GEdge %d near axis point",
                   axisEdge->tag());
    }
  }

  // Export axis edge + axis-column outer vertices so the caller can reclassify
  // them from gf to axisEdge AFTER _deleteUnusedVertices (where reparamOnFace
  // is no longer needed and the vertices still live in gf->mesh_vertices).
  if(axisEdge) {
    axisEdgeOut = axisEdge;
    for(int k = 1; k <= nbLayers_; k++)
      axisColVertsOut.push_back(grid[N][k]);
  }

  // Collect all BC face-interior vertices including the axis column.
  // The axis-column verts stay classified on gf until the caller reclassifies.
  for(int i = 0; i <= N; i++)
    for(int k = 1; k <= nbLayers_; k++)
      if(grid[i][k] && !blVerts.count(grid[i][k]))
        verts.insert(grid[i][k]);

  // Outer boundary MLines (k = nbLayers_ row, used as re-triangulation constraint)
  for(int i = 0; i < N; i++)
    outerLines.push_back(new MLine(grid[i][nbLayers_], grid[i + 1][nbLayers_]));

  if(useStitch) {
    // Replace BC column 0 with a single stitch quad that:
    //   • cancels the BL right-wall edge (blOuterAtJunction → junction)
    //   • cancels the first arc_bc domain segment (junction → baseVerts[bcStart+1])
    //   • bridges the BL outer row to the BC outer row (grid[0][nbLayers_])
    // The remaining pocket (from grid[0][nbLayers_] down to baseVerts[bcStart+1]
    // and back up via BC column 1's left wall) is one column wide and is
    // triangulated cleanly by the Delaunay without touching the junction profile.
    bcQuads.push_back(new MQuadrangle(
      blOuterAtJunction, baseVerts[bcStart],
      baseVerts[bcStart + 1], grid[0][nbLayers_]));
    for(int i = 1; i < N; i++)
      for(int k = 0; k < nbLayers_; k++)
        bcQuads.push_back(new MQuadrangle(
          grid[i][k],         grid[i + 1][k],
          grid[i + 1][k + 1], grid[i][k + 1]));
  }
  else {
    for(int i = 0; i < N; i++)
      for(int k = 0; k < nbLayers_; k++)
        bcQuads.push_back(new MQuadrangle(
          grid[i][k],         grid[i + 1][k],
          grid[i + 1][k + 1], grid[i][k + 1]));
  }

  Msg::Warning("BoundaryCorner (pre-mesh): %d columns x %d layers = %d quads (face %d), axisEdge=%d",
            N, nbLayers_, N * nbLayers_, gf->tag(), axisEdge ? axisEdge->tag() : -1);
  return true;
}

void BoundaryCornerField::buildCornerColumns(GModel *gm)
{
  computeParameters();
  if(curvesList_.empty()) return;
  if(w0max_ <= 0.0 || nbLayers_ < 1) {
    Msg::Error("BoundaryCorner: MaxColumnWidth must be > 0 and NbLayers >= 1");
    return;
  }

  GEdge *ge = gm->getEdgeByTag(curvesList_.front());
  if(!ge) return;

  // Find the GFace that owns this GEdge
  GFace *gf = nullptr;
  for(auto it = gm->firstFace(); it != gm->lastFace(); ++it) {
    for(auto *e : (*it)->edges()) {
      if(e == ge) { gf = *it; break; }
    }
    if(gf) break;
  }
  if(!gf) return;

  // t_end at axisPoint; t_start at the non-axis GVertex endpoint
  double t_end = arcLengthToParam(ge, axisPoint_[0], axisPoint_[1]);

  GVertex *gvBegin = ge->getBeginVertex();
  GVertex *gvEnd   = ge->getEndVertex();
  auto sqd = [&](GVertex *gv) -> double {
    if(!gv) return 1e30;
    double dx = gv->x() - axisPoint_[0], dy = gv->y() - axisPoint_[1];
    return dx*dx + dy*dy;
  };
  SPoint2 ptStart;
  double t_start;
  if(sqd(gvEnd) > sqd(gvBegin) && gvEnd) {
    ptStart = SPoint2(gvEnd->x(), gvEnd->y());
    t_start = arcLengthToParam(ge, gvEnd->x(), gvEnd->y());
  } else if(gvBegin) {
    ptStart = SPoint2(gvBegin->x(), gvBegin->y());
    t_start = arcLengthToParam(ge, gvBegin->x(), gvBegin->y());
  } else {
    t_start = ge->parBounds(0).low();
    GPoint gp0 = ge->point(t_start);
    ptStart = SPoint2(gp0.x(), gp0.y());
  }
  double tSign = (t_end >= t_start) ? 1.0 : -1.0;

  // --- Full profile arc length (numerical integration of |firstDer|) ---
  double S_total = 0.0;
  {
    const int nInt = 200;
    double dti = (t_end - t_start) / nInt;
    for(int i = 0; i < nInt; i++)
      S_total += ge->firstDer(t_start + (i + 0.5) * dti).norm() * std::abs(dti);
  }
  if(S_total < 1e-14) {
    Msg::Error("BoundaryCorner: full profile arc length is zero");
    return;
  }

  // --- Compressed section: NbLengthControl + NbHeightControl columns from
  // axisPoint outward (legacy: no length/height split, kept for reference) ---
  int N_corner = std::max(1, nbLengthControl_ + nbHeightControl_);
  double eps;
  const double q = lBLeff_ / w0max_;
  if(N_corner < 2 || q >= 1.0 - 1e-10)
    eps = 1.0;
  else
    eps = std::pow(q, 1.0 / (N_corner - 1));

  double S_corner_nom = (std::abs(eps - 1.0) < 1e-10)
      ? w0max_ * N_corner
      : w0max_ * (1.0 - std::pow(eps, N_corner)) / (1.0 - eps);

  // --- Constant-width section: remaining arc covered with w0max_ columns ---
  double S_rem = S_total - S_corner_nom;
  int N_constant = (S_rem > 1e-14) ? (int)std::ceil(S_rem / w0max_) : 0;
  int N = N_constant + N_corner;

  // Scale compressed widths to span S_total - N_constant*w0max_ exactly
  double S_corner_actual = S_total - (double)N_constant * w0max_;
  if(S_corner_actual < 1e-14) S_corner_actual = S_total / N_corner;
  double w0_corner = (std::abs(eps - 1.0) < 1e-10)
      ? S_corner_actual / N_corner
      : S_corner_actual * (1.0 - eps) / (1.0 - std::pow(eps, N_corner));

  // --- Build profile points: N_constant uniform columns then N_corner compressed ---
  std::vector<SPoint2> prof(N + 1);
  prof[0] = ptStart;
  for(int i = 1; i <= N; i++) {
    // constant zone: columns 1..N_constant (indices 0..N_constant-1 from start, widest first)
    // compressed zone: columns N_constant+1..N (indices 0..N_corner-1 from constant end, widest first)
    double Li;
    if(i <= N_constant)
      Li = w0max_;
    else
      Li = w0_corner * std::pow(eps, i - N_constant - 1);

    double t0    = arcLengthToParam(ge, prof[i - 1].x(), prof[i - 1].y());
    double speed = ge->firstDer(t0).norm();
    double dt    = (speed > 1e-14) ? Li / speed : 0.0;
    double t1    = (tSign > 0) ? std::min(t0 + tSign * dt, t_end)
                               : std::max(t0 + tSign * dt, t_end);
    GPoint gp    = ge->point(t1);
    prof[i]      = SPoint2(gp.x(), gp.y());
  }
  prof[N] = SPoint2(axisPoint_[0], axisPoint_[1]);

  // --- Build 2D grid: grid[i][k] at column i, height k ---
  // height k=0 on profile, k=nbLayers_ at BL-top offset
  std::vector<std::vector<MVertex *>> grid(N + 1,
    std::vector<MVertex *>(nbLayers_ + 1, nullptr));

  for(int i = 0; i <= N; i++) {
    double ti = arcLengthToParam(ge, prof[i].x(), prof[i].y());
    SPoint2 ni = normalAtPoint(ge, ti);
    for(int k = 0; k <= nbLayers_; k++) {
      double hk = (k == 0) ? 0.0
                : (std::abs(ratio_ - 1.0) < 1e-10)
                    ? h1_ * omega_ * k
                    : h1_ * omega_ * (std::pow(ratio_, k) - 1.0) / (ratio_ - 1.0);
      double x = prof[i].x() + ni.x() * hk;
      double y = prof[i].y() + ni.y() * hk;
      MVertex *v = new MVertex(x, y, 0.0, gf);
      gf->mesh_vertices.push_back(v);
      grid[i][k] = v;
    }
  }

  // --- Exact polygon boundary of the quad block (CCW) ---
  std::vector<SPoint2> blockPoly;
  for(int k = 0; k <= nbLayers_; k++)      // left side: profile → outer
    blockPoly.push_back(SPoint2(grid[0][k]->x(), grid[0][k]->y()));
  for(int i = 1; i <= N; i++)              // outer boundary: left → right
    blockPoly.push_back(SPoint2(grid[i][nbLayers_]->x(), grid[i][nbLayers_]->y()));
  for(int k = nbLayers_ - 1; k >= 0; k--) // right side: outer → profile
    blockPoly.push_back(SPoint2(grid[N][k]->x(), grid[N][k]->y()));
  for(int i = N - 1; i >= 1; i--)         // profile: right → left (closes polygon)
    blockPoly.push_back(SPoint2(grid[i][0]->x(), grid[i][0]->y()));

  auto inPoly = [&](double px, double py) -> bool {
    bool inside = false;
    int np = (int)blockPoly.size();
    for(int a = 0, b = np - 1; a < np; b = a++) {
      double xa = blockPoly[a].x(), ya = blockPoly[a].y();
      double xb = blockPoly[b].x(), yb = blockPoly[b].y();
      if(((ya > py) != (yb > py)) &&
         (px < (xb - xa) * (py - ya) / (yb - ya) + xa))
        inside = !inside;
    }
    return inside;
  };

  // --- Delete triangles inside the quad block; track their vertices ---
  std::set<MVertex *> deletedVSet;
  {
    std::vector<MTriangle *> keep;
    for(auto *tri : gf->triangles) {
      double cx = (tri->getVertex(0)->x() + tri->getVertex(1)->x() +
                   tri->getVertex(2)->x()) / 3.0;
      double cy = (tri->getVertex(0)->y() + tri->getVertex(1)->y() +
                   tri->getVertex(2)->y()) / 3.0;
      if(inPoly(cx, cy)) {
        for(int j = 0; j < 3; j++) deletedVSet.insert(tri->getVertex(j));
        delete tri;
      } else {
        keep.push_back(tri);
      }
    }
    gf->triangles = keep;
  }

  // --- Delete pre-existing quads inside the BC block (BL fan quads in the BC zone) ---
  // The BL algorithm generates fan quads at the arc_bl endpoint that span into the
  // BC zone. Only delete quads whose centroid is inside blockPoly; fan quads that
  // straddle the left boundary (centroid in the BL zone) are kept — they fill the
  // visual column at the BL/BC junction.
  {
    std::vector<MQuadrangle *> keepQ;
    for(auto *q : gf->quadrangles) {
      double cx = (q->getVertex(0)->x() + q->getVertex(1)->x() +
                   q->getVertex(2)->x() + q->getVertex(3)->x()) / 4.0;
      double cy = (q->getVertex(0)->y() + q->getVertex(1)->y() +
                   q->getVertex(2)->y() + q->getVertex(3)->y()) / 4.0;
      if(inPoly(cx, cy)) {
        for(int j = 0; j < 4; j++) deletedVSet.insert(q->getVertex(j));
        delete q;
      } else {
        keepQ.push_back(q);
      }
    }
    gf->quadrangles = keepQ;
  }

  // --- Merge grid[0][k] with adjacent BL outer vertices (conforming interface) ---
  // The BL last column's right-side outer vertices are at the same positions as
  // grid[0][1..nbLayers_]. Reusing those vertex objects makes the BL/BC interface
  // conforming. BL outer vertices may be on a BL-outer GEdge (dim=1), not the GFace,
  // so we search ALL quad vertices (no entity filter) with a loose tolerance.
  {
    const double tol2 = h1_ * h1_ * 0.01;  // (0.1 * h1)^2 — BL/BC positions differ by ~0.0003
    for(int k = 1; k <= nbLayers_; k++) {
      MVertex *gv = grid[0][k];
      double bestD2 = 1e30; MVertex *bestV = nullptr;
      for(auto *q : gf->quadrangles) {
        for(int j = 0; j < 4; j++) {
          MVertex *v = q->getVertex(j);
          if(v->onWhat() && v->onWhat()->dim() == 0) continue;  // skip model vertices
          double dx = v->x() - gv->x(), dy = v->y() - gv->y();
          double d2 = dx*dx + dy*dy;
          if(d2 < bestD2) { bestD2 = d2; bestV = v; }
        }
      }
      if(bestD2 < tol2) grid[0][k] = bestV;
    }
  }

  // --- Expand deletion to include fringe triangles ---
  // Fringe triangles straddle the BC block boundary: their centroid is outside
  // blockPoly but they share at least one vertex with the originally deleted zone.
  // Deleting them gives a clean, gap-free cavity for constrained retriangulation.
  const std::set<MVertex *> deletedVSet0 = deletedVSet;  // snapshot before expansion
  {
    std::vector<MTriangle *> keep;
    for(auto *tri : gf->triangles) {
      bool isFringe = false;
      for(int j = 0; j < 3 && !isFringe; j++)
        if(deletedVSet0.count(tri->getVertex(j))) isFringe = true;
      if(isFringe) {
        for(int j = 0; j < 3; j++) deletedVSet.insert(tri->getVertex(j));
        delete tri;
      }
      else keep.push_back(tri);
    }
    gf->triangles = keep;
  }

  // --- Collect outer cavity boundary edges from the surviving triangulation ---
  // A boundary edge of the surviving mesh that has at least one vertex in the
  // deleted zone is a cavity outer boundary edge (it borders the cavity from outside).
  std::vector<MEdge> constraints;
  {
    std::map<MEdge, int, Less_Edge> edgeCnt;
    for(auto *tri : gf->triangles)
      for(int j = 0; j < 3; j++)
        edgeCnt[tri->getEdge(j)]++;
    for(auto it = edgeCnt.begin(); it != edgeCnt.end(); ++it) {
      if(it->second == 1 &&
         (deletedVSet.count(it->first.getVertex(0)) ||
          deletedVSet.count(it->first.getVertex(1))))
        constraints.push_back(it->first);
    }
  }
  // BC outer boundary edges (top row of the BC block) — must be in the triangulation.
  for(int i = 0; i < N; i++)
    constraints.push_back(MEdge(grid[i][nbLayers_], grid[i + 1][nbLayers_]));
  // BC left and right column edges — prevent the triangulation from reaching inside
  // the BC block through the side columns.
  for(int k = 0; k < nbLayers_; k++) {
    constraints.push_back(MEdge(grid[0][k],     grid[0][k + 1]));
    constraints.push_back(MEdge(grid[N][k],     grid[N][k + 1]));
  }

  // --- Collect cavity vertices for re-triangulation ---
  // Include all freed vertices (deletedVSet) except those still held by surviving
  // BL/BC quads, plus the full BC outer boundary and column anchor nodes.
  std::set<MVertex *> survivingQuadVerts;
  for(auto *q : gf->quadrangles)
    for(int j = 0; j < 4; j++) survivingQuadVerts.insert(q->getVertex(j));

  std::set<MVertex *> cavityVSet;
  for(MVertex *v : deletedVSet)
    if(!survivingQuadVerts.count(v)) cavityVSet.insert(v);
  // BC outer boundary and side column nodes are constraint anchors; include them
  // regardless of whether they appeared in deletedVSet.
  for(int i = 0; i <= N; i++)
    cavityVSet.insert(grid[i][nbLayers_]);
  for(int k = 0; k <= nbLayers_; k++) {
    cavityVSet.insert(grid[0][k]);
    cavityVSet.insert(grid[N][k]);
  }

  std::vector<MVertex *> cavityVerts(cavityVSet.begin(), cavityVSet.end());

  // --- Re-triangulate the cavity using constrained Delaunay ---
  // delaunayMeshIn2D builds a fresh Delaunay triangulation of cavityVerts and
  // recovers all constraint edges by diagonal swaps.
  std::vector<MTriangle *> newTris;
  delaunayMeshIn2D(cavityVerts, newTris, /*removeBox=*/true, &constraints);

  // Discard any new triangle whose centroid falls inside the BC block (those cells
  // will be filled by BC quads) and add the rest to the face.
  int nNewTri = 0;
  for(auto *tri : newTris) {
    double cx = (tri->getVertex(0)->x() + tri->getVertex(1)->x() +
                 tri->getVertex(2)->x()) / 3.0;
    double cy = (tri->getVertex(0)->y() + tri->getVertex(1)->y() +
                 tri->getVertex(2)->y()) / 3.0;
    if(inPoly(cx, cy)) { delete tri; }
    else { gf->triangles.push_back(tri); nNewTri++; }
  }

  // --- Insert BC structured quads ---
  for(int i = 0; i < N; i++)
    for(int k = 0; k < nbLayers_; k++)
      gf->quadrangles.push_back(new MQuadrangle(
        grid[i][k],         grid[i + 1][k],
        grid[i + 1][k + 1], grid[i][k + 1]));

  // --- Final vertex purge: remove orphaned vertices ---
  {
    std::set<MVertex *> keep;
    for(auto *t : gf->triangles)
      for(int j = 0; j < 3; j++) keep.insert(t->getVertex(j));
    for(auto *q : gf->quadrangles)
      for(int j = 0; j < 4; j++) keep.insert(q->getVertex(j));
    std::vector<MVertex *> finalV;
    for(auto *v : gf->mesh_vertices) {
      if(keep.count(v)) finalV.push_back(v);
      else delete v;
    }
    gf->mesh_vertices = finalV;
  }

  Msg::Info(
    "BoundaryCorner: %d columns × %d layers = %d quads, "
    "cavity retriangulation: %d new triangles (face %d)",
    N, nbLayers_, N * nbLayers_, nNewTri, gf->tag());
}

double BoundaryCornerField::operator()(double x, double y, double z,
                                       GEntity *ge)
{
  computeParameters();
  double dx   = x - axisPoint_[0];
  double dy   = y - axisPoint_[1];
  double dist = std::sqrt(dx * dx + dy * dy);

  // Zone A: constant ColWidth arc-length over NbHeightControl columns,
  // innermost, adjacent to AxisPoint.
  double S_A = std::max(0, nbHeightControl_) * lBLeff_;

  // Zone B: geometric arc-length transition from ColWidth up to MaxColumnWidth
  // over NbLengthControl columns, outer part of the corner zone.
  int nl = std::max(1, nbLengthControl_);
  double eps_op;
  const double q_op = (w0max_ > 1e-100) ? lBLeff_ / w0max_ : 1.0;
  if(nl < 2 || q_op >= 1.0 - 1e-10)
    eps_op = 1.0;
  else
    eps_op = std::pow(q_op, 1.0 / (nl - 1));
  double S_B = (std::abs(eps_op - 1.0) < 1e-10)
      ? w0max_ * nl
      : w0max_ * (1.0 - std::pow(eps_op, nl)) / (1.0 - eps_op);

  if(dist < S_A)
    return lBLeff_;
  if(S_B > 0.0 && dist < S_A + S_B)
    return lBLeff_ + (w0max_ - lBLeff_) * (dist - S_A) / S_B;

  // Zone C: geometric growth from MaxColumnWidth up to FarSize, smoothing the
  // handoff to the background mesh size instead of returning unconstrained
  // (1e22) right at the edge of zone B.
  if(farSize_ > 0.0 && farSize_ > w0max_) {
    double distC = dist - (S_A + S_B);
    double size = w0max_ * std::pow(farRatio_, distC / w0max_);
    return std::min(farSize_, size);
  }
  return 1e22;
}

// ---------------------------------------------------------------------------
// BoundaryDoubleCornerField
// ---------------------------------------------------------------------------

std::string BoundaryDoubleCornerField::getDescription()
{
  return "Structured quad columns covering a full profile that meets the "
         "symmetry axis at 90 degrees at both ends (nose and tail). "
         "Compressed corner columns at each axis endpoint, constant-width "
         "columns in between.";
}

BoundaryDoubleCornerField::BoundaryDoubleCornerField()
  : h1_(1e-3), ratio_(1.15), w0max_(0.1), lBL_(-1.0), omega_(1.0),
    nbLayers_(30), nbLengthControl_(5), nbHeightControl_(5), skipAxisColumn_(0),
    farSize_(-1.0), farRatio_(1.2), lBLeff_(1e-3), hTotal_(0.0)
{
  nosePoint_[0] = nosePoint_[1] = 0.0;
  tailPoint_[0] = tailPoint_[1] = 0.0;

  options["CurvesList"] = new FieldOptionList(
    curvesList_, "Tags of the profile curves", &update_needed);
  options["NosePoint"] = new FieldOptionListDouble(
    nosePointList_, "Nose axis point [x_nose, 0.0]", &update_needed);
  options["TailPoint"] = new FieldOptionListDouble(
    tailPointList_, "Tail axis point [x_tail, 0.0]", &update_needed);
  options["Size"] = new FieldOptionDouble(
    h1_, "Height of the first BL row", &update_needed);
  options["Ratio"] = new FieldOptionDouble(
    ratio_, "Geometric growth ratio between successive BL rows", &update_needed);
  options["NbLayers"] = new FieldOptionInt(
    nbLayers_, "Number of BL rows", &update_needed);
  options["NbLengthControl"] = new FieldOptionInt(
    nbLengthControl_,
    "Number of columns, outer part of each corner zone, over which the column "
    "arc-length transitions geometrically from MaxColumnWidth down to ColWidth.  "
    "Row heights in this region stay BL-geometric (no blending).",
    &update_needed);
  options["NbHeightControl"] = new FieldOptionInt(
    nbHeightControl_,
    "Number of columns, inner part of each corner zone (adjacent to the axis "
    "point), over which the column arc-length is held constant at ColWidth.  "
    "Row heights stay BL-geometric all the way to the axis point.",
    &update_needed);
  options["MaxColumnWidth"] = new FieldOptionDouble(
    w0max_,
    "Column arc-length in the constant-width mid-section.",
    &update_needed);
  options["ColWidth"] = new FieldOptionDouble(
    lBL_,
    "Arc-length of the innermost BC column at each axis point.  "
    "-1 = use the height of the last (outermost) BL row.",
    &update_needed);
  options["Omega"] = new FieldOptionDouble(
    omega_, "Height scale factor for BC rows (default 1.0)", &update_needed);
  options["SkipAxisColumn"] = new FieldOptionInt(
    skipAxisColumn_,
    "Set to 1 to omit the structured quad columns exactly on y=0. "
    "Use for 3D wedge revolve meshes.",
    &update_needed);
  options["FarSize"] = new FieldOptionDouble(
    farSize_,
    "Target far-field size that the size field grows toward geometrically "
    "(zone C, beyond MaxColumnWidth) at FarRatio, smoothing the transition to "
    "the background mesh size instead of jumping straight to it.  "
    "<= 0 disables zone C (legacy behaviour: unconstrained beyond zone B).",
    &update_needed);
  options["FarRatio"] = new FieldOptionDouble(
    farRatio_, "Geometric growth ratio used in zone C (FarSize transition)",
    &update_needed);
}

void BoundaryDoubleCornerField::computeParameters()
{
  double hLastRow = h1_ * std::pow(ratio_, nbLayers_ - 1);
  lBLeff_ = (lBL_ > 0.0) ? lBL_ : hLastRow;
  if(lBLeff_ < 1e-100) lBLeff_ = 1e-100;

  {
    auto it = nosePointList_.begin();
    nosePoint_[0] = (it != nosePointList_.end()) ? *it++ : 0.0;
    nosePoint_[1] = (it != nosePointList_.end()) ? *it   : 0.0;
  }
  {
    auto it = tailPointList_.begin();
    tailPoint_[0] = (it != tailPointList_.end()) ? *it++ : 0.0;
    tailPoint_[1] = (it != tailPointList_.end()) ? *it   : 0.0;
  }

  if(std::abs(ratio_ - 1.0) < 1e-10)
    hTotal_ = h1_ * nbLayers_;
  else
    hTotal_ = h1_ * (std::pow(ratio_, nbLayers_) - 1.0) / (ratio_ - 1.0);
}

double BoundaryDoubleCornerField::arcLengthToParam(GEdge *ge, double x, double y)
{
  double t = 0.0;
  const SPoint3 p(x, y, 0.0);
  ge->closestPoint(p, t);
  return t;
}

SPoint2 BoundaryDoubleCornerField::normalAtPoint(GEdge *ge, double t)
{
  SVector3 d = ge->firstDer(t);
  double len = d.norm();
  if(len < 1e-14) return SPoint2(0.0, 1.0);
  return SPoint2(-d.y() / len, d.x() / len);
}

void BoundaryDoubleCornerField::subdivideAxisEdge(
  GEdge *axisEdge, MVertex *axisVert,
  const std::vector<MVertex *> &axisColVerts,
  GFace *gf)
{
  // axisColVerts = grid[end][1..nbLayers_], where grid[end][nbLayers_] is outermost.
  // We replace axis MLines that fall in [xNose, xOuter] with the chain:
  //   vReconnect → axisColVerts[nbLayers_-1] → ... → axisColVerts[0] → axisVert
  double xNose  = axisVert->x();
  double xOuter = axisColVerts.back()->x();
  double xMin = std::min(xNose, xOuter);
  double xMax = std::max(xNose, xOuter);

  MVertex *vReconnect = nullptr;
  std::vector<int> toRemove;
  for(int li = 0; li < (int)axisEdge->lines.size(); li++) {
    MLine *ml   = axisEdge->lines[li];
    MVertex *va = ml->getVertex(0), *vb = ml->getVertex(1);
    bool aIn = (va->x() >= xMin - 1e-14 && va->x() <= xMax + 1e-14);
    bool bIn = (vb->x() >= xMin - 1e-14 && vb->x() <= xMax + 1e-14);
    if(aIn && bIn) {
      toRemove.push_back(li);
    } else if(aIn != bIn) {
      toRemove.push_back(li);
      vReconnect = aIn ? vb : va;
    }
  }

  if(vReconnect && !toRemove.empty()) {
    int insertPos = toRemove.front();
    for(int i = (int)toRemove.size() - 1; i >= 0; i--) {
      delete axisEdge->lines[toRemove[i]];
      axisEdge->lines.erase(axisEdge->lines.begin() + toRemove[i]);
    }

    int nb = (int)axisColVerts.size();  // == nbLayers_
    // Chain: vReconnect → axisColVerts[nb-1] → ... → axisColVerts[0] → axisVert
    axisEdge->lines.insert(axisEdge->lines.begin() + insertPos,
                           new MLine(vReconnect, axisColVerts[nb - 1]));
    for(int k = nb - 2; k >= 0; k--)
      axisEdge->lines.insert(axisEdge->lines.begin() + ++insertPos,
                             new MLine(axisColVerts[k + 1], axisColVerts[k]));
    axisEdge->lines.insert(axisEdge->lines.begin() + ++insertPos,
                           new MLine(axisColVerts[0], axisVert));

    auto &mv = axisEdge->mesh_vertices;
    auto rmBegin = std::remove_if(mv.begin(), mv.end(),
      [xMin, xMax](MVertex *v) {
        return v->x() >= xMin - 1e-14 && v->x() <= xMax + 1e-14;
      });
    // NOTE: these vertices are NOT deleted here, even though they are being
    // dropped from axisEdge's own bookkeeping. At this point gf still holds
    // its *initial* (pre-BoundaryCorner) triangulation, built before this
    // field ran, and that triangulation can still reference these exact
    // MVertex objects as triangle corners (the axis curve is part of gf's
    // boundary, so its 1D-mesh vertices are also gf's 2D boundary vertices).
    // Deleting them here frees memory that's about to be read again by the
    // surrounding meshGenerator() re-triangulation pass (modifyInitialMesh-
    // ForBoundaryCorners -> deMeshGFace -> meshGenerator), producing a
    // use-after-free crash with a corrupted MVertex vtable. The initial
    // triangulation (and these orphaned vertices with it) is discarded
    // wholesale a few lines later in the caller via deMeshGFace/gf->
    // deleteMesh(), so simply leaking them here (bounded: at most a couple
    // of vertices per axis corner, once per mesh generation) is safe.
    mv.erase(rmBegin, mv.end());
  } else {
    Msg::Warning("BoundaryDoubleCorner: failed to subdivide axis GEdge %d",
                 axisEdge->tag());
  }
}

bool BoundaryDoubleCornerField::buildForFace(
  GFace *gf,
  const std::vector<MQuadrangle *> &blQuads,
  const std::set<MVertex *> &blVerts,
  std::vector<MQuadrangle *> &bcQuads,
  std::set<MVertex *> &verts,
  std::vector<MLine *> &outerLines,
  std::map<MVertex *, std::vector<MVertex *>> &junctionMap,
  std::vector<std::pair<GEdge *, std::vector<MVertex *>>> &axisReclassifyOut)
{
  computeParameters();
  if(curvesList_.empty() || nbLayers_ < 1) return false;

  GEdge *ge = nullptr;
  for(int tag : curvesList_) {
    GEdge *e = gf->model()->getEdgeByTag(tag);
    if(!e) continue;
    for(auto *fe : gf->edges()) if(fe == e) { ge = e; break; }
    if(ge) break;
  }
  if(!ge) return false;

  // Build baseVerts: GVertex-start + mesh_vertices + GVertex-end
  std::vector<MVertex *> baseVerts;
  if(ge->getBeginVertex() && !ge->getBeginVertex()->mesh_vertices.empty())
    baseVerts.push_back(ge->getBeginVertex()->mesh_vertices[0]);
  for(auto *v : ge->mesh_vertices) baseVerts.push_back(v);
  if(ge->getEndVertex() && !ge->getEndVertex()->mesh_vertices.empty())
    baseVerts.push_back(ge->getEndVertex()->mesh_vertices[0]);
  if((int)baseVerts.size() < 2) {
    Msg::Error("BoundaryDoubleCorner: GEdge %d has no 1D mesh yet", ge->tag());
    return false;
  }

  // Deduplicate consecutive baseVerts at the same (x,y) position.  This can
  // occur when a LC-spaced intermediate vertex lands exactly on the next
  // scatter-point knot, giving two coincident entries in ge->mesh_vertices and
  // producing zero-area base faces in the quad columns.
  {
    const double tol2 = 1e-20;
    std::vector<MVertex *> dedup;
    dedup.reserve(baseVerts.size());
    dedup.push_back(baseVerts[0]);
    for(int i = 1; i < (int)baseVerts.size(); i++) {
      double dx = baseVerts[i]->x() - dedup.back()->x();
      double dy = baseVerts[i]->y() - dedup.back()->y();
      if(dx*dx + dy*dy > tol2) dedup.push_back(baseVerts[i]);
    }
    baseVerts = std::move(dedup);
  }
  if((int)baseVerts.size() < 2) {
    Msg::Error("BoundaryDoubleCorner: GEdge %d: all baseVerts coincide", ge->tag());
    return false;
  }

  // Orient baseVerts so nosePoint_ is at front, tailPoint_ at back.
  auto sqDist2 = [](MVertex *v, double x, double y) {
    double dx = v->x() - x, dy = v->y() - y;
    return dx*dx + dy*dy;
  };
  {
    MVertex *vF = baseVerts.front(), *vB = baseVerts.back();
    double dFN = sqDist2(vF, nosePoint_[0], nosePoint_[1]);
    double dBN = sqDist2(vB, nosePoint_[0], nosePoint_[1]);
    if(dBN < dFN)
      std::reverse(baseVerts.begin(), baseVerts.end());
  }

  // Detect BL-consumed vertices at the nose end (bcStart) and tail end (bcEnd).
  std::set<MVertex *> blOuterVerts;
  for(auto *q : blQuads)
    for(int j = 2; j <= 3; j++) blOuterVerts.insert(q->getVertex(j));

  int bcStart = 0;
  for(int i = 1; i < (int)baseVerts.size(); i++) {
    if(blOuterVerts.count(baseVerts[i])) bcStart = i;
    else break;
  }

  int bcEnd = 0;
  for(int i = (int)baseVerts.size() - 2; i >= 0; i--) {
    if(blOuterVerts.count(baseVerts[i])) bcEnd = (int)baseVerts.size() - 1 - i;
    else break;
  }

  // Find axis GVertices and GEdges at nose and tail.
  auto findAxisEdge = [&](MVertex *axisVert, double axisY) -> std::pair<GVertex *, GEdge *> {
    GVertex *axisGV = nullptr;
    for(auto *gv : {ge->getBeginVertex(), ge->getEndVertex()}) {
      if(!gv || gv->mesh_vertices.empty()) continue;
      if(gv->mesh_vertices[0] == axisVert) { axisGV = gv; break; }
    }
    GEdge *axisEdge = nullptr;
    if(axisGV) {
      for(auto *adj : axisGV->edges()) {
        if(adj == ge) continue;
        GVertex *other = (adj->getBeginVertex() == axisGV) ? adj->getEndVertex()
                                                            : adj->getBeginVertex();
        if(other && std::abs(other->y() - axisY) < 1e-10) { axisEdge = adj; break; }
      }
    }
    return {axisGV, axisEdge};
  };

  MVertex *noseVert = baseVerts.front();
  MVertex *tailVert = baseVerts.back();
  auto [noseGV, noseAxisEdge] = findAxisEdge(noseVert, nosePoint_[1]);
  auto [tailGV, tailAxisEdge] = findAxisEdge(tailVert, tailPoint_[1]);

  // N = number of inter-vertex segments to cover with quad columns.
  // When axis edges are found we include the endpoint columns; otherwise skip them
  // (as BCF does) and leave a small gap for Delaunay.
  bool includeNose = (noseAxisEdge && !skipAxisColumn_);
  bool includeTail = (tailAxisEdge && !skipAxisColumn_);

  // Useful baseVerts range: [bcStart .. M-1-bcEnd].
  // Number of segments in that range = M - 1 - bcStart - bcEnd.
  // Including both endpoints means all segments; excluding one endpoint removes
  // one column from that end (start or finish the walk one step inside).
  int M = (int)baseVerts.size();
  int N = M - 1 - bcStart - bcEnd;
  if(!includeNose) N--;   // skip noseVert column
  if(!includeTail) N--;   // skip tailVert column
  if(N < 1) {
    Msg::Warning("BoundaryDoubleCorner: no segments remain after BL fan for GEdge %d",
                 ge->tag());
    return false;
  }

  // First base-vert index: bcStart when includeNose, bcStart+1 when not.
  int startIdx = includeNose ? bcStart : bcStart + 1;

  auto baseVertForCol = [&](int i) -> MVertex * {
    return baseVerts[startIdx + i];
  };

  // Origin-independent reference for "outward" direction: compare nose/tail x
  // to the profile's own mean x rather than to x=0, so a geometry shifted away
  // from the origin (e.g. body starting at x=1) still gets the correct sign.
  double xMeanProfile = 0.0;
  for(auto *v : baseVerts) xMeanProfile += v->x();
  xMeanProfile /= (double)baseVerts.size();

  // Determine normal orientation using nose point (must point outward from axis).
  bool flipNormal = false;
  {
    double xN = nosePoint_[0] - xMeanProfile;
    if(std::abs(xN) > 1e-10) {
      double tN = arcLengthToParam(ge, noseVert->x(), noseVert->y());
      SPoint2 nN = normalAtPoint(ge, tN);
      flipNormal = (nN.x() * xN < -1e-10);
    }
  }

  std::vector<std::vector<MVertex *>> grid(N + 1,
    std::vector<MVertex *>(nbLayers_ + 1, nullptr));

  for(int i = 0; i <= N; i++) {
    MVertex *bv = baseVertForCol(i);
    grid[i][0] = bv;

    double ti = arcLengthToParam(ge, bv->x(), bv->y());
    SPoint2 ni = normalAtPoint(ge, ti);
    if(flipNormal) ni = SPoint2(-ni.x(), -ni.y());
    double bx = bv->x(), by = bv->y();

    // Force axis columns onto y=0 with exact horizontal normal.
    bool isNoseCol = (includeNose && i == 0);
    bool isTailCol = (includeTail && i == N);
    if(isNoseCol) {
      double signX = (nosePoint_[0] >= xMeanProfile) ? 1.0 : -1.0;
      ni = SPoint2(signX, 0.0);
      by = nosePoint_[1];
    } else if(isTailCol) {
      double signX = (tailPoint_[0] >= xMeanProfile) ? 1.0 : -1.0;
      ni = SPoint2(signX, 0.0);
      by = tailPoint_[1];
    }

    for(int k = 1; k <= nbLayers_; k++) {
      double hk = (std::abs(ratio_ - 1.0) < 1e-10)
                      ? h1_ * omega_ * k
                      : h1_ * omega_ * (std::pow(ratio_, k) - 1.0) / (ratio_ - 1.0);
      grid[i][k] = new MVertex(bx + ni.x() * hk, by + ni.y() * hk, 0.0, gf);
    }
  }

  // Subdivide axis edges and record for later reclassification.
  if(includeNose) {
    std::vector<MVertex *> noseColVerts(nbLayers_);
    for(int k = 1; k <= nbLayers_; k++) noseColVerts[k - 1] = grid[0][k];
    subdivideAxisEdge(noseAxisEdge, noseVert, noseColVerts, gf);
    axisReclassifyOut.push_back({noseAxisEdge, std::move(noseColVerts)});
  }
  if(includeTail) {
    std::vector<MVertex *> tailColVerts(nbLayers_);
    for(int k = 1; k <= nbLayers_; k++) tailColVerts[k - 1] = grid[N][k];
    subdivideAxisEdge(tailAxisEdge, tailVert, tailColVerts, gf);
    axisReclassifyOut.push_back({tailAxisEdge, std::move(tailColVerts)});
  }

  // Collect interior BC vertices.
  for(int i = 0; i <= N; i++)
    for(int k = 1; k <= nbLayers_; k++)
      if(grid[i][k] && !blVerts.count(grid[i][k]))
        verts.insert(grid[i][k]);

  // Outer boundary MLines (k = nbLayers_ row).
  for(int i = 0; i < N; i++)
    outerLines.push_back(new MLine(grid[i][nbLayers_], grid[i + 1][nbLayers_]));

  // Build quads.
  for(int i = 0; i < N; i++)
    for(int k = 0; k < nbLayers_; k++)
      bcQuads.push_back(new MQuadrangle(
        grid[i][k],         grid[i + 1][k],
        grid[i + 1][k + 1], grid[i][k + 1]));

  Msg::Warning("BoundaryDoubleCorner (pre-mesh): %d columns x %d layers = %d quads"
               " (face %d), noseAxisEdge=%d tailAxisEdge=%d",
               N, nbLayers_, N * nbLayers_, gf->tag(),
               noseAxisEdge ? noseAxisEdge->tag() : -1,
               tailAxisEdge ? tailAxisEdge->tag() : -1);
  return true;
}

double BoundaryDoubleCornerField::operator()(double x, double y, double z,
                                             GEntity *ge)
{
  computeParameters();

  auto zoneDist = [&](double px, double py) {
    double dx = x - px, dy = y - py;
    return std::sqrt(dx*dx + dy*dy);
  };
  double dist = std::min(zoneDist(nosePoint_[0], nosePoint_[1]),
                         zoneDist(tailPoint_[0], tailPoint_[1]));

  // Zone A: constant ColWidth arc-length over NbHeightControl columns,
  // innermost, adjacent to the axis point.
  double S_A = std::max(0, nbHeightControl_) * lBLeff_;

  // Zone B: geometric arc-length transition from ColWidth up to MaxColumnWidth
  // over NbLengthControl columns, outer part of the corner zone.
  int nl = std::max(1, nbLengthControl_);
  double eps_op;
  const double q_op = (w0max_ > 1e-100) ? lBLeff_ / w0max_ : 1.0;
  if(nl < 2 || q_op >= 1.0 - 1e-10)
    eps_op = 1.0;
  else
    eps_op = std::pow(q_op, 1.0 / (nl - 1));
  double S_B = (std::abs(eps_op - 1.0) < 1e-10)
      ? w0max_ * nl
      : w0max_ * (1.0 - std::pow(eps_op, nl)) / (1.0 - eps_op);

  if(dist < S_A)
    return lBLeff_;
  if(S_B > 0.0 && dist < S_A + S_B)
    return lBLeff_ + (w0max_ - lBLeff_) * (dist - S_A) / S_B;

  if(farSize_ > 0.0 && farSize_ > w0max_) {
    double distC = dist - (S_A + S_B);
    double size = w0max_ * std::pow(farRatio_, distC / w0max_);
    return std::min(farSize_, size);
  }
  return 1e22;
}

