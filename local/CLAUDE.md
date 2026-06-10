# BOUNDARY CORNER — IMPLEMENTATION REFERENCE

## Context

2D axisymmetric CFD (OpenFOAM wedge). Fluid is **outside** the body. `y=0` is the symmetry axis. A body profile **P** meets the axis at 90° at the **corner C**. `BoundaryCornerField` inserts a structured quad fan in the region between a StartPoint **S** on P and **C**, replacing the unstructured triangles that would otherwise span that corner.

## Key geometry

- **S** = StartPoint: on the profile, where BL ends and BC begins. Adjacent BL tangential cell size = l_BL.
- **C** = AxisPoint: where P meets y=0 at 90°. Last quad column sits on the axis.
- **Columns** run along the arc S→C (arc-length compression, decreasing widths).
- **Layers** run outward along the surface normal (BL-style geometric growth).
- Last column (at C) lies on y=0 because the normal at C is horizontal (+x).

## Options (Field.h / Field.cpp constructor)

| Option | Member | Semantics |
|---|---|---|
| `CurvesList` | `curvesList_` | GEdge tags of the profile |
| `AxisPoint` | `axisPointList_` | [x_C, 0.0] — corner |
| `StartPoint` | `startPointList_` | [x_S, y_S] — BC zone start |
| `Size` | `h1_` | First BL layer normal height |
| `Ratio` | `ratio_` | BL layer geometric ratio |
| `NbLayers` | `nbLayers_` | Number of BL rows (k direction) |
| `NbCornerColumns` | `nbCornerColumns_` | Number of quad columns (i direction) |
| `Delta1` | `delta1_` | **First** column arc-length at S (= l_BL, large). -1 → use ColWidth |
| `ColWidth` | `lBL_` | **Last** column arc-length at C (corner cell, small). -1 → use Size |
| `Omega` | `omega_` | Layer height scale factor (default 1.0) |

**Parameter semantics (important):**
- `Delta1` is the **large** BL-matching value at S; `Size`/`ColWidth` is the **small** corner value at C.
- Compression ratio: `eps = (lBLeff / d1eff)^(1/(N-1))` — requires `Delta1 > ColWidth/Size` for `eps < 1`.
- `lBLeff = lBL_ if lBL_ > 0 else h1_`; `d1eff = delta1_ if delta1_ > 0 else lBLeff`.

## Class layout (Mesh/Field.h)

```cpp
class BoundaryCornerField : public Field {
  // Options (stored as members, bound via FieldOption*):
  std::list<int>    curvesList_;
  std::list<double> axisPointList_, startPointList_;
  double h1_, ratio_, delta1_, lBL_, omega_;
  int    nbLayers_, nbCornerColumns_;

  // Computed:
  double axisPoint_[2], startPoint_[2];
  double eps_;      // tangential compression ratio (<1 = compress toward C)
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

### `buildCornerColumns(gm)` — main mesh injection

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

## Integration points

- **`Mesh/Generator.cpp`**: after the `for(GFace*)` 2D mesh loop, iterates `FieldManager` for `BoundaryCornerField` instances and calls `buildCornerColumns(gm)`.
- **`Mesh/Field.cpp` ~L3295**: `BoundaryCorner` registered in `FieldManager::map_type_name`.
- **`Context.h`, `Context.cpp`**: `CTX::mesh::boundaryCornerField` int option.
- **`Options.cpp`, `Options.h`, `DefaultOptions.h`**: `Mesh.BoundaryCornerField` Python-accessible option.
- **`test_bc.py`**: test script — half-ellipse nose, H1=0.012, RATIO=1.20, N_LAY=6, N_COL=7, DELTA1=0.08. Expects 42 quads.

## Build

```bash
cd build && cmake .. -DENABLE_MMG3D=OFF && make -j$(nproc)
```
mmg3d disabled: pre-existing linker bug in bundled v4.0, unrelated to this work.

## Known issues / open work

- **Non-conforming interface**: quad outer boundary nodes do not generally coincide with surrounding triangle nodes. The fringe snap is a heuristic. Proper fix: enforce quad boundary edges as Delaunay constraints, or pre-declare embedded edges before `meshGFace()`.
- `buildCornerColumns` only processes `curvesList_.front()` (first curve). Multi-curve support not implemented.
- **interface gap**: There is currently a column gap in interfaces between BLs & BCs.
