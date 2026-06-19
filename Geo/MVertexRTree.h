// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef MVERTEX_RTREE
#define MVERTEX_RTREE

#include <vector>
#include <cfloat>
#include <cmath>
#include "GmshMessage.h"
#include "MVertex.h"
#include "MElement.h"
#include "Context.h"
#include "rtree.h"

// Accumulates the smallest squared edge length found in a set of mesh
// elements (lines, triangles, quadrangles, ...) into minLen2. Used to size
// the matching tolerance of an MVertexRTree off the *local* mesh feature
// size of the entity actually being extruded, rather than the global model
// bounding box (CTX::instance()->lc) -- see localVertexRTreeTolerance below.
template <typename ElemPtr>
inline void accumulateMinEdgeLength2(const std::vector<ElemPtr> &elems,
                                     double &minLen2)
{
  for(std::size_t i = 0; i < elems.size(); i++) {
    MElement *e = elems[i];
    int n = (int)e->getNumVertices();
    for(int j = 0; j < n; j++) {
      MVertex *v0 = e->getVertex(j);
      MVertex *v1 = e->getVertex((j + 1) % n);
      double dx = v1->x() - v0->x();
      double dy = v1->y() - v0->y();
      double dz = v1->z() - v0->z();
      double len2 = dx * dx + dy * dy + dz * dz;
      if(len2 > 0. && len2 < minLen2) minLen2 = len2;
    }
  }
}

// Absolute position tolerance for matching vertices of an extrusion. An
// extruded entity's own mesh can contain elements many orders of magnitude
// smaller than the full model extent (e.g. thin boundary-layer quads inside
// a much larger far-field domain). Using geom.tolerance*globalLc in that
// case can exceed the local mesh spacing and merge vertices that are
// actually distinct (e.g. two different boundary-layer rows), producing
// degenerate elements. Anchor the tolerance on whichever is smaller: the
// global bounding-box diagonal, or the smallest real element edge found via
// accumulateMinEdgeLength2.
inline double localVertexRTreeTolerance(double minLen2)
{
  double lc = CTX::instance()->lc;
  if(minLen2 < DBL_MAX) {
    double minLen = std::sqrt(minLen2);
    if(minLen > 0. && minLen < lc) lc = minLen;
  }
  return CTX::instance()->geom.tolerance * lc;
}

// Stores MVertex pointers in an R-Tree so we can query unique vertices by their
// coordinates, up to a prescribed tolerance.
class MVertexRTree {
private:
  RTree<MVertex *, double, 3, double> *_rtree;
  double _tol;
  static bool rtree_callback(MVertex *v, void *ctx)
  {
    MVertex **out = static_cast<MVertex **>(ctx);
    *out = v;
    return false; // we're done searching
  }

public:
  MVertexRTree(double tolerance = 1.e-8)
  {
    _rtree = new RTree<MVertex *, double, 3, double>();
    _tol = tolerance;
  }
  ~MVertexRTree()
  {
    _rtree->RemoveAll();
    delete _rtree;
  }
  MVertex *insert(MVertex *v, bool warnIfExists = false,
                  std::set<MVertex *> *duplicates = 0)
  {
    MVertex *out;
    double _min[3] = {v->x() - _tol, v->y() - _tol, v->z() - _tol};
    double _max[3] = {v->x() + _tol, v->y() + _tol, v->z() + _tol};
    if(!_rtree->Search(_min, _max, rtree_callback, &out)) {
      _rtree->Insert(_min, _max, v);
      return 0;
    }
    else {
      if(duplicates) {
        duplicates->insert(out);
        duplicates->insert(v);
      }
      if(warnIfExists) {
        Msg::Warning("Vertex %d (%.16g, %.16g, %.16g) already exists in the "
                     "mesh with tolerance %g: Vertex %d (%.16g, %.16g, %.16g)",
                     v->getNum(), v->x(), v->y(), v->z(), _tol, out->getNum(),
                     out->x(), out->y(), out->z());
      }
      return out;
    }
  }
  int insert(std::vector<MVertex *> &v, bool warnIfExists = false,
             std::set<MVertex *> *duplicates = 0)
  {
    int num = 0;
    for(std::size_t i = 0; i < v.size(); i++)
      num += (insert(v[i], warnIfExists, duplicates) ? 1 : 0);
    return num; // number of vertices not inserted
  }
  MVertex *find(double x, double y, double z)
  {
    MVertex *out;
    double _min[3] = {x - _tol, y - _tol, z - _tol};
    double _max[3] = {x + _tol, y + _tol, z + _tol};
    if(_rtree->Search(_min, _max, rtree_callback, &out)) return out;
    return 0;
  }
  std::size_t size() { return _rtree->Count(); }
};

#endif
