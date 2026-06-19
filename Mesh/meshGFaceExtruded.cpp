// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include <set>
#include "GmshConfig.h"
#include "GModel.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MQuadrangle.h"
#include "ExtrudeParams.h"
#include "MVertexRTree.h"
#include "Context.h"
#include "GmshMessage.h"

#if defined(HAVE_QUADTRI)
#include "QuadTriExtruded2D.h"
#endif

static void addTriangle(MVertex *v1, MVertex *v2, MVertex *v3, GFace *to)
{
  to->triangles.push_back(new MTriangle(v1, v2, v3));
}

static void addQuadrangle(MVertex *v1, MVertex *v2, MVertex *v3, MVertex *v4,
                          GFace *to)
{
  to->quadrangles.push_back(new MQuadrangle(v1, v2, v3, v4));
}

static void
createQuaTri(std::vector<MVertex *> &v, GFace *to,
             std::set<std::pair<MVertex *, MVertex *> > *constrainedEdges,
             MLine *source, int tri_quad_flag)
{
  ExtrudeParams *ep = to->meshAttributes.extrude;
  if(v[0] == v[1] || v[1] == v[3])
    addTriangle(v[0], v[3], v[2], to);
  else if(v[0] == v[2] || v[2] == v[3])
    addTriangle(v[0], v[1], v[3], to);
  else if(v[0] == v[3] || v[1] == v[2])
    Msg::Error("Uncoherent extruded quadrangle in surface %d", to->tag());
  else {
    // Trevor Strickler added the tri_quad_flag stuff here.
    if((ep->mesh.Recombine && tri_quad_flag != 2) || tri_quad_flag == 1) {
      if(!constrainedEdges)
        addQuadrangle(v[0], v[1], v[3], v[2], to);
      else {
        std::pair<MVertex *, MVertex *> p1(std::min(v[1], v[2]),
                                           std::max(v[1], v[2]));
        std::pair<MVertex *, MVertex *> p2(std::min(v[0], v[3]),
                                           std::max(v[0], v[3]));
        if(constrainedEdges->count(p1)) {
          addTriangle(v[2], v[1], v[0], to);
          addTriangle(v[2], v[3], v[1], to);
        }
        else if(constrainedEdges->count(p2)) {
          addTriangle(v[2], v[3], v[0], to);
          addTriangle(v[0], v[3], v[1], to);
        }
        else
          addQuadrangle(v[0], v[1], v[3], v[2], to);
      }
    }
    else if(!constrainedEdges) {
      addTriangle(v[0], v[1], v[3], to);
      addTriangle(v[0], v[3], v[2], to);
    }
    else {
      std::pair<MVertex *, MVertex *> p(std::min(v[1], v[2]),
                                        std::max(v[1], v[2]));
      if(constrainedEdges->count(p)) {
        addTriangle(v[2], v[1], v[0], to);
        addTriangle(v[2], v[3], v[1], to);
      }
      else {
        addTriangle(v[2], v[3], v[0], to);
        addTriangle(v[0], v[3], v[1], to);
      }
    }
  }
}

static void
extrudeMesh(GEdge *from, GFace *to, MVertexRTree &pos,
            std::set<std::pair<MVertex *, MVertex *> > *constrainedEdges)
{
  ExtrudeParams *ep = to->meshAttributes.extrude;

  // create vertices (if the edges are constrained, they already exist)
  if(!constrainedEdges) {
    for(std::size_t i = 0; i < from->mesh_vertices.size(); i++) {
      std::vector<MVertex *> extruded_vertices;
      MVertex *v = from->mesh_vertices[i];
      MEdgeVertex *mv = (MEdgeVertex *)v;
      mv->bl_data = new MVertexBoundaryLayerData();
      for(int j = 0; j < ep->mesh.NbLayer; j++) {
        for(int k = 0; k < ep->mesh.NbElmLayer[j]; k++) {
          double x = v->x(), y = v->y(), z = v->z();
          ep->Extrude(j, k + 1, x, y, z);
          if(j != ep->mesh.NbLayer - 1 || k != ep->mesh.NbElmLayer[j] - 1) {
            MVertex *newv = 0;
            if(to->geomType() != GEntity::DiscreteSurface &&
               to->geomType() != GEntity::BoundaryLayerSurface) {
              SPoint3 xyz(x, y, z);
              SPoint2 uv = to->parFromPoint(xyz);
              newv = new MFaceVertex(x, y, z, to, uv[0], uv[1]);
            }
            else {
              newv = new MVertex(x, y, z, to);
            }
            to->mesh_vertices.push_back(newv);
            pos.insert(newv);
            extruded_vertices.push_back(newv);
          }
        }
      }
      mv->bl_data->addChildrenFamily(extruded_vertices);
    }
  }

  int tri_quad_flag = 0;

#if defined(HAVE_QUADTRI)
  // figure out whether to recombine this surface or not in the event of
  // quadToTri region neighbors (if QuadToTri, tri_quad_flag is an int flag that
  // lets createQuadTri() override the surface's intrinsic ep->mesh.Recombine
  // flag.  tri_quad_flag values: 0 = no override, 1 = mesh with quads, 2 = mesh
  // with triangles.)
  bool detectQuadToTriLateral = false;
  bool quadToTri_valid =
    IsValidQuadToTriLateral(to, &tri_quad_flag, &detectQuadToTriLateral);
  if(detectQuadToTriLateral && !quadToTri_valid)
    Msg::Error(
      "In MeshGFaceExtrudedSurface::extrudeMesh(), Mesh of QuadToTri Lateral "
      "surface %d likely has errors.",
      to->tag());
#endif

  // create elements (note that it would be faster to access the *interior*
  // nodes by direct indexing, but it's just simpler to query everything by
  // position)
  for(std::size_t i = 0; i < from->lines.size(); i++) {
    MVertex *v0 = from->lines[i]->getVertex(0);
    MVertex *v1 = from->lines[i]->getVertex(1);
    for(int j = 0; j < ep->mesh.NbLayer; j++) {
      for(int k = 0; k < ep->mesh.NbElmLayer[j]; k++) {
        std::vector<MVertex *> verts;
        double x[4] = {v0->x(), v1->x(), v0->x(), v1->x()};
        double y[4] = {v0->y(), v1->y(), v0->y(), v1->y()};
        double z[4] = {v0->z(), v1->z(), v0->z(), v1->z()};
        for(int p = 0; p < 2; p++) {
          ep->Extrude(j, k, x[p], y[p], z[p]);
          ep->Extrude(j, k + 1, x[p + 2], y[p + 2], z[p + 2]);
        }
        for(int p = 0; p < 4; p++) {
          MVertex *tmp = pos.find(x[p], y[p], z[p]);
          if(!tmp) {
            Msg::Error("Could not find extruded vertex (%.16g, %.16g, %.16g) "
                       "in surface %d",
                       x[p], y[p], z[p], to->tag());
            return;
          }
          verts.push_back(tmp);
        }
        createQuaTri(verts, to, constrainedEdges, from->lines[i],
                     tri_quad_flag);
      }
    }
  }
}

static void copyMesh(GFace *from, GFace *to, MVertexRTree &pos)
{
  ExtrudeParams *ep = to->meshAttributes.extrude;

  // interior vertices
  std::vector<MVertex *> mesh_vertices = from->mesh_vertices;

  // add all embedded vertices
  std::vector<MVertex *> embedded = from->getEmbeddedMeshVertices();
  mesh_vertices.insert(mesh_vertices.end(), embedded.begin(), embedded.end());

  // create extruded vertices
  for(std::size_t i = 0; i < mesh_vertices.size(); i++) {
    MVertex *v = mesh_vertices[i];
    double x = v->x(), y = v->y(), z = v->z();
    ep->Extrude(ep->mesh.NbLayer - 1, ep->mesh.NbElmLayer[ep->mesh.NbLayer - 1],
                x, y, z);
    MVertex *newv = 0;
    if(to->geomType() != GEntity::DiscreteSurface &&
       to->geomType() != GEntity::BoundaryLayerSurface) {
      SPoint3 xyz(x, y, z);
      SPoint2 uv = to->parFromPoint(xyz);
      newv = new MFaceVertex(x, y, z, to, uv[0], uv[1]);
    }
    else {
      newv = new MVertex(x, y, z, to);
    }
    to->mesh_vertices.push_back(newv);
    to->correspondingVertices[newv] = v;
    pos.insert(newv);
  }

  std::vector<double> tfo;
  //ep->GetAffineTransform(tfo); // TODO: check transform
  to->GEntity::setMeshMaster(from, tfo, false);

#if defined(HAVE_QUADTRI)
  // if performing QuadToTri mesh, cannot simply copy the mesh from the source.
  // The vertices and triangles can be copied directly though.  First, of
  // course, do some checks and make sure this is a valid QuadToTri top surface
  // before engaging in QuadToTri meshing.
  int quadToTri = NO_QUADTRI;
  bool detectQuadToTriTop = false;
  int quadToTri_valid =
    IsValidQuadToTriTop(to, &quadToTri, &detectQuadToTriTop);
  bool is_toroidal = quadToTri_valid >= 2 ? true : false;
  bool is_noaddverts = quadToTri_valid == 3 ? true : false;
  if(detectQuadToTriTop && !quadToTri_valid && !is_toroidal) {
    Msg::Error("In MeshGFaceExtrudedSurface::copyMesh(), Mesh of QuadToTri top "
               "surface %d likely has errors.",
               to->tag());
  }

  // if this is toroidal No New Vertices QuadToTri, then replace the root
  // dependency face's boundary quads with triangles for better meshing.
  if(is_toroidal && is_noaddverts) {
    GFace *root = findRootSourceFaceForFace(from);
    if(root == from) {
      ReplaceBndQuadsInFace(root);
      Msg::Warning(
        "To facilitate QuadToTri interface on surface %d, source "
        "surface %d was re-meshed with all triangles on boundary. "
        "To avoid this, use QuadTriAddVerts instead of QuadTriNoNewVerts",
        to->tag(), root->tag());
    }
  }
#endif

  // create triangle elements
  for(std::size_t i = 0; i < from->triangles.size(); i++) {
    std::vector<MVertex *> verts;
    for(int j = 0; j < 3; j++) {
      MVertex *v = from->triangles[i]->getVertex(j);
      double x = v->x(), y = v->y(), z = v->z();
      ep->Extrude(ep->mesh.NbLayer - 1,
                  ep->mesh.NbElmLayer[ep->mesh.NbLayer - 1], x, y, z);
      MVertex *tmp = pos.find(x, y, z);
      if(!tmp) {
        Msg::Error(
          "Could not find extruded vertex (%.16g, %.16g, %.16g) in surface %d",
          x, y, z, to->tag());
        return;
      }
      verts.push_back(tmp);
    }
    addTriangle(verts[0], verts[1], verts[2], to);
  }

#if defined(HAVE_QUADTRI)
  // Add triangles for divided quads for QuadTri. If quadtotri and not part of a
  // toroidal extrusion, mesh the top surface accordingly
  if(detectQuadToTriTop && !is_toroidal) {
    if(!MeshQuadToTriTopSurface(from, to, pos))
      Msg::Error("In MeshExtrudedSurface()::copyMesh(), mesh of QuadToTri top "
                 "surface %d failed.",
                 to->tag());
    return;
  }
#endif

  // create quadrangle elements if NOT QuadToTri and NOT toroidal
  for(std::size_t i = 0; i < from->quadrangles.size(); i++) {
    std::vector<MVertex *> verts;
    for(int j = 0; j < 4; j++) {
      MVertex *v = from->quadrangles[i]->getVertex(j);
      double x = v->x(), y = v->y(), z = v->z();
      ep->Extrude(ep->mesh.NbLayer - 1,
                  ep->mesh.NbElmLayer[ep->mesh.NbLayer - 1], x, y, z);
      MVertex *tmp = pos.find(x, y, z);
      if(!tmp) {
        Msg::Error(
          "Could not find extruded vertex (%.16g, %.16g, %.16g) in surface %d",
          x, y, z, to->tag());
        return;
      }
      verts.push_back(tmp);
    }
    addQuadrangle(verts[0], verts[1], verts[2], verts[3], to);
  }
}

int MeshExtrudedSurface(
  GFace *gf, std::set<std::pair<MVertex *, MVertex *> > *constrainedEdges)
{
  ExtrudeParams *ep = gf->meshAttributes.extrude;

  if(!ep || !ep->mesh.ExtrudeMesh) return 0;

  Msg::Info("Meshing surface %d (extruded)", gf->tag());

  // resolve the source entity first, so the rtree's matching tolerance can
  // be sized off its own mesh's local feature size (smallest element edge)
  // rather than the global model bounding box -- see localVertexRTreeTolerance
  // in MVertexRTree.h for why this matters (thin boundary-layer elements
  // inside a much larger far-field domain).
  GEdge *fromEdge = 0;
  GFace *fromFace = 0;
  if(ep->geo.Mode == EXTRUDED_ENTITY) {
    fromEdge = gf->model()->getEdgeByTag(std::abs(ep->geo.Source));
    if(!fromEdge) {
      Msg::Error("Unknown source curve %d for extrusion", ep->geo.Source);
      return 0;
    }
  }
  else {
    fromFace = gf->model()->getFaceByTag(std::abs(ep->geo.Source));
    if(!fromFace) {
      Msg::Error("Unknown source surface %d for extrusion", ep->geo.Source);
      return 0;
    }
    else if(fromFace->geomType() != GEntity::DiscreteSurface &&
            fromFace->meshStatistics.status != GFace::DONE) {
      // cannot mesh the face yet (the source face is not meshed):
      // will do it later
      return 1;
    }
  }

  double minLen2 = DBL_MAX;
  if(fromEdge) accumulateMinEdgeLength2(fromEdge->lines, minLen2);
  if(fromFace) {
    accumulateMinEdgeLength2(fromFace->triangles, minLen2);
    accumulateMinEdgeLength2(fromFace->quadrangles, minLen2);
  }
  MVertexRTree pos(localVertexRTreeTolerance(minLen2));

  // insert all the vertices on the boundary of the face gf
  std::vector<GEdge *> const &edges = gf->edges();
  std::vector<GEdge *>::const_iterator it = edges.begin();
  while(it != edges.end()) {
    pos.insert((*it)->mesh_vertices);
    if((*it)->getBeginVertex())
      pos.insert((*it)->getBeginVertex()->mesh_vertices);
    if((*it)->getEndVertex())
      pos.insert((*it)->getEndVertex()->mesh_vertices);
    ++it;
  }

  // if the edges of the mesh are constrained, the vertices already
  // exist on the face--so we add them to the set
  if(constrainedEdges) {
    pos.insert(gf->mesh_vertices);
  }

  if(fromEdge)
    extrudeMesh(fromEdge, gf, pos, constrainedEdges);
  else
    copyMesh(fromFace, gf, pos);

  gf->meshStatistics.status = GFace::DONE;
  return 1;
}
