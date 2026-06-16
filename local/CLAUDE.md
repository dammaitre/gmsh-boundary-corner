# BOUNDARY CORNER — IMPLEMENTATION REFERENCE

## Context

2D axisymmetric CFD (OpenFOAM wedge). Fluid is **outside** the body. `y=0` is the symmetry axis. A body profile **P** meets the axis at 90° at the **corner C**. `BoundaryCornerField` covers the **entire** profile with structured quad columns: a geometrically compressed section near **C** and a constant-width section for the rest.

## Key geometry

- **C** = AxisPoint: where P meets y=0 at 90°. Last quad column sits on the axis.
- **Compressed zone**: `NbCornerColumns` columns near C, arc-length decreasing from `MaxColumnWidth` to `ColWidth/Size`.
- **Constant zone**: remaining profile covered with constant `MaxColumnWidth` columns.
- **Layers** run outward along the surface normal (BL-style geometric growth) for all columns.
- Last column (at C) lies on y=0 because the normal at C is horizontal (+x).

## Options (Field.h / Field.cpp constructor)

| Option | Member | Semantics |
|---|---|---|
| `CurvesList` | `curvesList_` | GEdge tags of the full profile |
| `AxisPoint` | `axisPointList_` | [x_C, 0.0] — corner |
| `Size` | `h1_` | First BL layer normal height |
| `Ratio` | `ratio_` | BL layer geometric ratio |
| `NbLayers` | `nbLayers_` | Number of BL rows (k direction) |
| `NbCornerColumns` | `nbCornerColumns_` | Number of compressed columns near AxisPoint |
| `MaxColumnWidth` | `w0max_` | Column arc-length in the constant zone; max in the compressed zone |
| `ColWidth` | `lBL_` | Arc-length of the innermost column at AxisPoint. -1 → use Size |
| `Omega` | `omega_` | Layer height scale factor (default 1.0) |

**Parameter semantics (important):**
- `ColWidth`/`Size` is the **small** corner value; `MaxColumnWidth` is the **large** transition value.
- Compression ratio: `eps = (lBLeff / w0max_)^(1/(NbCornerColumns-1))`.
- `lBLeff = lBL_ if lBL_ > 0 else h1_`.

## Class layout (Mesh/Field.h)

```cpp
class BoundaryCornerField : public Field {
  // Options (stored as members, bound via FieldOption*):
  std::list<int>    curvesList_;
  std::list<double> axisPointList_;
  double h1_, ratio_, lBL_, omega_, w0max_;
  int    nbLayers_, nbCornerColumns_;

  // Computed:
  double axisPoint_[2];
  double lBLeff_;   // resolved ColWidth (lBL_ > 0 ? lBL_ : h1_)
  double hTotal_;   // total BL height

  void    computeParameters();
  double  arcLengthToParam(GEdge*, double x, double y);
  SPoint2 normalAtPoint(GEdge*, double t);
};
```

## Implementation (Mesh/Field.cpp ~L2987)

### `computeParameters()`
- Resolves `lBLeff_` and `d1eff` (no mutation of stored options).
- Guards: `lBLeff_ >= 1e-100`, `ratio==1` → linear `hTotal_`.
- `eps_ = pow(lBLeff_/d1eff, 1/(N-1))`.

### `arcLengthToParam(ge, x, y)` — critical
```cpp
const SPoint3 p(x, y, 0.0);  // const is mandatory
ge->closestPoint(p, t);       // uses virtual overload: GPoint(const SPoint3&, double&)
```
`const` on `p` is required to resolve the **virtual** `closestPoint(const SPoint3&, double&)` overload (t = output param). Without `const`, the non-virtual `SPoint3 closestPoint(SPoint3&, double tolerance)` overload is called with tolerance=0 → `closestPointFinder(0)` → `oversample(pts, 0)` → `d/0 = INT_MAX` iterations → **OOM / SIGKILL**.

### `normalAtPoint(ge, t)`
Returns outward surface normal: rotate tangent 90° CCW → `(-dy/|d|, dx/|d|)`. Correct for CW-parameterized profiles (e.g. ellipse arc from top to nose).

### `buildCornerColumns(gm)` — legacy post-mesh entry (kept for reference, superseded by `buildForFace`)

Called **post-mesh** (from `Generator.cpp` after the 2D mesh loop). Cannot be pre-mesh: `meshGFace()` calls `GFace::deleteMesh()` internally, wiping any pre-injected elements.

**Steps:**
1. Resolve GEdge + GFace from `curvesList_`.
2. Compute `t_start`, `t_end` via `arcLengthToParam`.
3. Integrate arc length `S_total` (200-point midpoint rule on `|firstDer|`).
4. Scale: `w0 = S_total*(1-eps)/(1-eps^N)` — first column width; series sums exactly to `S_total`.
5. Step N-1 profile points along the curve: `dt = Li / |firstDer(t0)|`, clamped to `t_end`.  `prof[N]` forced to axisPoint.
6. For each profile point, compute outward normal → place `(N+1)*(nbLayers+1)` MVertex on GFace.
7. Insert `N*nbLayers` MQuadrangle.
8. Build CCW `blockPoly` from actual grid vertices (left side → outer top → right side → profile base).
9. Delete triangles whose centroid is inside `blockPoly` (ray-cast); collect `deletedVSet`.
10. Snap fringe vertices (in `deletedVSet ∩ survivingVSet`, not in `gridSet`) to nearest quad outer-boundary vertex.
11. Delete degenerate triangles (snap can collapse a triangle to a line).
12. Purge unreferenced `mesh_vertices`.

### `operator()(x, y, z)`
Background field size hint (used before mesh generation):
```
dist = distance to AxisPoint
return lBLeff_ + (d1eff - lBLeff_) * dist / (hTotal_*5)   for dist < zone
return 1e22                                                  otherwise
```
Small size near C, BL-matching size at S.

## Pre-mesh architecture (current — branch `snap-remove`)

The old post-mesh `buildCornerColumns` has been replaced by a **pre-mesh** approach that mirrors how `BoundaryLayerField` works: structured quads are injected *before* the far-field triangulation, so the triangle mesh naturally conforms to the quad outer boundary.

### Key change: `buildForFace` (Mesh/Field.cpp)

```cpp
bool BoundaryCornerField::buildForFace(
    GFace *gf,
    const std::vector<MQuadrangle *> &blQuads,
    const std::set<MVertex *>        &blVerts,
    std::vector<MQuadrangle *>       &bcQuads,
    std::set<MVertex *>              &verts,
    std::vector<MLine *>             &outerLines);
```

Called from `modifyInitialMeshForBoundaryCorners` (meshGFace.cpp), itself called inside the outer `meshGenerator` after `modifyInitialMeshForBoundaryLayers`.

**Steps:**
1. Resolve `ge` (the BC profile edge) from `curvesList_`.
2. Build `baseVerts` = GVertex-start + `ge->mesh_vertices` + GVertex-end, oriented S→C.
3. Detect `bcStart`: how many leading `baseVerts[1..k]` are already consumed by BL fan quads (they appear as `j=2/3` vertices of `blQuads`). Start BC columns from `baseVerts[bcStart]`.
4. Set `N = baseVerts.size() - 2 - bcStart`: exclude the last column at the axis point (see Axis Column below).
5. Build `grid[i][k]`: base at `baseVerts[bcStart+i]`, outer rows at `base + k*h_k * normal`.
6. Fill `bcQuads` (N×nbLayers_ quads), `verts` (BC inner+outer vertices), `outerLines` (outer-row MLines).

**Axis column included (N = size-1 instead of size-2):**
The last column at `theta=0°` is now included as structured quads. Its outer vertices land on `y=0` at `(x_C + h_k, 0)`. To make the XOR bedges close correctly, the axis GEdge's 1D mesh (`axisEdge->lines`) is subdivided to replace the MLines that would overlap the new column with a chain: `vReconnect → grid[N][nbLayers_] → ... → grid[N][1] → axisVert`. The XOR then cancels all shared edges between the BC quads and the new axis chain, leaving a clean closed bedge loop.

**Key classification rule:** `grid[N][k]` (axis-column outer vertices) are classified on `gf` (dim=2), NOT on `axisEdge`. This is required so `reparamMeshVertexOnFace` can compute correct (x,y) parametric coordinates for the inner BDS Delaunay mesher. Classifying them on `axisEdge` causes `reparamOnFace(t=0)` to return the farfield endpoint instead of the actual vertex position, breaking edge recovery.

### `modifyInitialMeshForBoundaryCorners` (Mesh/meshGFace.cpp)

Mirrors `modifyInitialMeshForBoundaryLayers`:
1. Calls `buildForFace` for every `BoundaryCornerField` registered in `FieldManager`.
2. Builds `bedges` by XOR: domain edges ⊕ BL quads ⊕ BC quads → the boundary of the remaining-to-triangulate region.
3. Creates `discreteEdge ne(444445)` with `ne.lines` = bedges MLines.
4. **Protects BL outer vertices from `deMeshGFace`** (see bug fix below).
5. Calls `deMeshGFace(gf)` then `meshGenerator(gf, 0, 0, true, false, &hop)` to triangulate the remaining region.
6. **Strips protected vertices from `gf->mesh_vertices`** so the outer mesher's `verts`/`bcVerts` re-insertion (lines 1907/1912) doesn't double-free them.

### Integration points

- `Mesh/meshGFace.cpp`: `modifyInitialMeshForBoundaryCorners` called at ~L1830 (inside outer `meshGenerator`, after `modifyInitialMeshForBoundaryLayers`).
- `Mesh/meshGFace.cpp`: outer `meshGenerator` inserts `bcQuads` and `bcVerts` at ~L1911-1913.
- `Mesh/Field.cpp` ~L3295: `BoundaryCorner` registered in `FieldManager::map_type_name`.
- `Context.h`, `Options.cpp`, `DefaultOptions.h`: `Mesh.BoundaryCornerField` option (set to the BC field id).

## Build

```bash
cd build && make -j$(nproc) shared   # Python API (used by test_bc.py)
cd build && make -j$(nproc) gmsh     # binary (used for --gui)
```
mmg3d disabled: pre-existing linker bug in bundled v4.0, unrelated to this work.

## Test

```bash
python3 test_bc.py          # headless — prints triangle/quad counts
python3 test_bc.py --gui    # opens result in gmsh GUI
```

Expected output (H1=0.012, RATIO=1.20, N_LAY=6, N_COLS=8, W0_MAX=0.06):
- 354 quadrangles (full-profile BC structured quads, including axis column)
- 2173 triangles (far-field Delaunay)

## Bug fixes (branch `snap-remove`)

### 1. Use-after-free of BL outer vertices (`meshGFace.cpp`)

**Root cause:** After the BL inner mesher runs, `MFaceVertex` objects (BL outer row) are in `gf->mesh_vertices`. `modifyInitialMeshForBoundaryCorners` calls `deMeshGFace(gf)` → `GFace::deleteMesh()` → frees every `gf->mesh_vertices[i]`. But `ne.lines` (already built) holds raw pointers to those freed vertices. The 3rd `meshGenerator` then dereferences them → crash/UB.

**Fix:** Before `deMeshGFace`, collect all `ne.lines` vertices with `onWhat()==gf` into `protected_verts` and remove them from `gf->mesh_vertices`. After the inner mesher, `_deleteUnusedVertices` re-adds them (they appear on the triangulation boundary); strip them again so the outer mesher's `verts`/`bcVerts` insert doesn't create duplicates and double-free.

```cpp
// Protect BL outer vertices from deMeshGFace
std::set<MVertex *> protected_verts;
for(auto *l : ne.lines)
    for(int j = 0; j < 2; j++)
        if(l->getVertex(j)->onWhat() == gf) protected_verts.insert(l->getVertex(j));
auto &mv = gf->mesh_vertices;
mv.erase(std::remove_if(mv.begin(), mv.end(),
    [&](MVertex *v){ return protected_verts.count(v); }), mv.end());

deMeshGFace kil; kil(gf);
meshGenerator(gf, 0, 0, true, false, &hop);

// Strip again to avoid duplicate re-insertion
mv.erase(std::remove_if(mv.begin(), mv.end(),
    [&](MVertex *v){ return protected_verts.count(v); }), mv.end());
```

### 2. Axis column as structured quads (`Field.cpp`)

**Root cause (original):** At the nose (theta=0°) the BC outer column falls exactly on y=0. Collinear vertices between `p_nose` and the nearest axis mesh vertex prevented constrained-Delaunay edge recovery.

**Fix (current):** The axis GEdge's 1D mesh is subdivided to split around the new column's outer vertices (`grid[N][k]`). The XOR bedge mechanism then cancels all shared edges, giving a clean closed boundary for the inner Delaunay. `grid[N][k]` are classified on `gf` (not on `axisEdge`) so `reparamMeshVertexOnFace` returns their correct (x,y) position rather than the farfield endpoint. Old axis interior vertices in the replaced range are deleted from `axisEdge->mesh_vertices`.

## Known limitations

- `buildForFace` only processes `curvesList_.front()` (first curve). Multi-curve support not implemented.
- Axis column inclusion requires a valid `axisEdge` (GEdge adjacent to axisGV lying on y=0). If not found, falls back to old behaviour: N = size-2, last column skipped.
