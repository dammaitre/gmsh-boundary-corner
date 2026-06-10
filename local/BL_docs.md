# Boundary Layer Fields in Gmsh — `Mesh/Field.h` and `Mesh/Field.cpp`

## 1. Class Overview

**Class:** `BoundaryLayerField` (inherits `Field`)  
**Header:** `Mesh/Field.h:178–236` | **Implementation:** `Mesh/Field.cpp:2626–2979`

---

## 2. Parameters and Data Members

| Parameter | Type | Purpose |
|---|---|---|
| `hwall_n` | double | Mesh size normal to the wall |
| `hfar` | double | Mesh size far from the wall |
| `ratio` | double | Growth ratio between successive layers |
| `thickness` | double | Max BL thickness |
| `tgt_aniso_ratio` (`AnisoMax`) | double | Threshold aniso ratio before creating a fan |
| `iRecombine` (`Quads`) | int | Generate quads in the BL |
| `iIntersect` (`IntersectMetrics`) | int | Intersect metrics from all attractors |
| `edges_id` | `list<int>` | Curve tags that receive a BL |
| `nodes_id` | `list<int>` | Point tags where the BL ends |
| `fan_nodes_id` | `list<int>` | Points where a fan (wedge) is inserted |
| `hwall_n_nodes` | `list<double>` | Per-node override of `hwall_n` |
| `excluded_faces_id` | `list<int>` | Surfaces to exclude from BL |
| `_att_fields` | `list<DistanceField*>` | Cached KD-tree distance attractors |
| `current_distance` | double | Distance at the last evaluated point |
| `_closest_point` | `SPoint3` | Closest attractor point found |

---

## 3. Attractor Initialization (Lazy, On-Demand)

The first time the field is evaluated (or after `removeAttractors()` resets `update_needed = true`), attractors are built:

- One `DistanceField(0, tag, 100000)` per **node** in `nodes_id` — 100k sample points
- One `DistanceField(1, tag, 300000)` per **curve** in `edges_id` — 300k sample points

These wrap KD-tree structures and are stored in `_att_fields` for reuse.

---

## 4. Scalar Size Computation — `operator()(x,y,z,ge)` (Field.cpp:2781)

1. Lazily initialize attractors (if `update_needed`).
2. Query every `DistanceField` in `_att_fields` → take the minimum `dist`.
3. If `dist > thickness * ratio` → return `1e22` (field has no influence here).
4. Return mesh size using the exponential growth law:

```
lc = dist * (ratio - 1) + hwall_n
return min(hfar, lc)
```

At `dist=0` this gives `hwall_n`; each successive layer is `ratio` times thicker than the previous.

---

## 5. Anisotropic Metric Computation — `operator()(x,y,z,SMetric3&,ge)` (Field.cpp:2938)

This is the main path used by the anisotropic mesher.

### Normal and tangential sizes

```
lc_n = dist * (ratio - 1) + hwall_n    // normal to wall — fine
lc_t = min(lc_n * AnisoMax, hfar)      // tangential — coarse
```

### Metric construction by attractor dimension

- **Vertex (dim=0):** direction toward vertex → `buildMetricTangentToCurve(t1, lc_n, lc_n)`
- **Curve (dim=1):** tangent + curvature of the edge at the closest point. Curvature adapts `lc_t` via:

  ```
  1/d² = 0.5/lc_t² * (1 + sqrt(1 + 4*crv²*lc_t⁴ / (lc_n²*β²)))
  ```

  → `buildMetricTangentToCurve(tangent, sqrt(d²), lc_n)`

- **Face (dim=2):** principal curvatures at the closest surface point, applied to both principal directions; aniso ratio clamped to `tgt_aniso_ratio` → `buildMetricTangentToSurface(...)`

### Metric intersection

- If `iIntersect=1`: intersect metrics from **all** attractors (conservative — takes the finest).
- Otherwise: keep only the closest attractor's metric.

---

## 6. Setup Methods — Scoping the Field to 1D/2D Context

### `setupFor1d(int iE)` — Field.cpp:2683

Called before meshing curve `iE`. Rebuilds `nodes_id` to contain only the BL-endpoint nodes belonging to that curve, so the 1D size field acts only along this curve.

### `setupFor2d(int iF)` — Field.cpp:2714

Called before meshing surface `iF`:

1. Skips if `iF` is in `excluded_faces_id`.
2. Walks all edges of the face and filters to curves that are **not** extruded entities and whose adjacent non-extruded face count == 1 (genuine BL walls, not shared interior faces).
3. Rebuilds `edges_id` and `nodes_id` for this surface scope.
4. Calls `removeAttractors()` to invalidate the KD-tree cache.

This scoping mechanism reconfigures the same field object before each surface is meshed, so distance computations only consider walls relevant to that surface.

---

## 7. What `BoundaryLayerField` Does NOT Do

**`Field.cpp` / `Field.h` contain no BL geometry creation.** The field is purely a **size/metric prescription**:

- Only **reads** from `GModel` (via `getEdgeByTag`, `getFaceByTag`, `getVertexByTag`) to retrieve tangent/curvature data.
- Does not call `GModel::add*()`, does not touch `ExtrudeParams`, does not create new curves, surfaces, or volumes.

The actual **creation of extruded BL geometry** (structured layers of quads/hexes, fan elements at corners, extrusion topology) is handled elsewhere — specifically in `Mesh/BoundaryLayers.cpp` and the `Mesh/Generator*.cpp` pipeline, which call `setupFor2d`, evaluate this field's metric, and then drive the extrusion.

---

## 8. FieldManager Integration

`FieldManager` (declared in `Field.h:132–171`) maintains:

- `_boundaryLayer_fields`: `vector<int>` of registered BL field IDs
- `addBoundaryLayerFieldId(int id)` / `addBoundaryLayerFieldId(vector<int>&)`
- `getBoundaryLayerField(int i)` / `getNumBoundaryLayerFields()`

The type is registered at startup:

```cpp
map_type_name["BoundaryLayer"] = new FieldFactoryT<BoundaryLayerField>();
```

---

## 9. Data Flow Summary

```
User sets: edges_id, nodes_id, hwall_n, ratio, thickness
       |
       v
Field registered in FieldManager._boundaryLayer_fields
       |
       v
Mesher calls setupFor2d(iF)  -->  scopes field to surface iF
       |
       v
Mesher evaluates field(x,y,z) or field.metric(x,y,z) at each mesh point
       |
       v
Distance to wall curves/points computed via KD-tree DistanceField attractors
       |
       v
Exponential size law --> scalar lc   OR   anisotropic SMetric3 (with curvature)
       |
       v
Meshing algorithm uses sizes to place vertices in structured BL layers
       |
       v
BL geometry extrusion handled in Mesh/BoundaryLayers.cpp (outside Field.cpp)
```
