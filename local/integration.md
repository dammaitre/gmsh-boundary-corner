# BoundaryDoubleCorner — Integration into AEROD/core/phase_mesh.py

## What changes and what stays

| Layer | Old | New |
|---|---|---|
| Body curves | 3 segments (nose cap, main, tail cap) | 1 single spline nose→top→tail |
| Axis layout | 4 lines with A_nose / A_tail breakpoints | 2 lines, no breakpoints |
| BL field | `BoundaryLayer` on main curve only | `BoundaryDoubleCorner` on full spline |
| Mesh generation | `generate(2)` → `recombine()` → `extrude` | `extrude` geometry first → `generate(3)` |
| Physical groups | Bounding-box scan of 3D surfaces | `curve_to_lat` mapping from extruded geometry |
| `convert_step` | unchanged | unchanged |
| `patch_fix_step` | fixes `axis` → `symmetryPlane` | rename `defaultFaces` → `symAxis`, fix `symmetry` type; drop `axis` fix |
| `validate_mesh` | unchanged | unchanged |

---

## Parameter mapping

| `build_mesh` param | BDC field option | Notes |
|---|---|---|
| `h1` | `Size` | first layer height — identical semantics |
| `n_bl` | `NbLayers` | layer count — identical semantics |
| `r_bl` | `Ratio` | growth ratio — identical semantics |
| `lc_wall` | `MaxColumnWidth` | tangential cell size in constant zone; must be > `h1` |
| — | `NbCornerColumns` | fan columns near each axis tip; **new param**, 15–25 for large Re |
| — | `ColWidth` | innermost column arc-length; defaults to `Size` when omitted |

There is no `Thickness` option — BDC infers the BL height from `Size`, `Ratio`, `NbLayers`.

---

## `build_mesh` — what to rewrite

### 1. Single spline, simple axis

Replace the three-curve split (nose cap / main / tail cap) and the four axis segments with:

```python
# profile_ptags is the full list: nose (y=0) → shoulder → tail (y=0)
CURVE_BODY = gmsh.model.geo.addSpline(profile_ptags)

x_nose = float(pts[0,  0])
x_tail = float(pts[-1, 0])

P_outlet_axis = gmsh.model.geo.addPoint(x_max, 0.0,   0.0, lc_far)
P_inlet_axis  = gmsh.model.geo.addPoint(x_min, 0.0,   0.0, lc_far)
P_inlet_top   = gmsh.model.geo.addPoint(x_min, y_max, 0.0, lc_far)
P_far_top     = gmsh.model.geo.addPoint(x_max, y_max, 0.0, lc_far)

ptag_nose = profile_ptags[0]
ptag_tail = profile_ptags[-1]

L_AX_R  = gmsh.model.geo.addLine(P_outlet_axis, ptag_nose)   # right axis
L_RIGHT = gmsh.model.geo.addLine(P_outlet_axis, P_far_top)
L_TOP   = gmsh.model.geo.addLine(P_far_top,     P_inlet_top)
L_LEFT  = gmsh.model.geo.addLine(P_inlet_top,   P_inlet_axis)
L_AX_L  = gmsh.model.geo.addLine(P_inlet_axis,  ptag_tail)   # left axis

loop = gmsh.model.geo.addCurveLoop(
    [-L_AX_R, L_RIGHT, L_TOP, L_LEFT, L_AX_L, -CURVE_BODY])
surf = gmsh.model.geo.addPlaneSurface([loop])
```

### 2. Extrude geometry before meshing

```python
DZ = 1.0
out   = gmsh.model.geo.extrude(
    [(2, surf)], 0, 0, DZ, numElements=[1], recombine=True)
s_top = out[0][1]   # back face
vol   = out[1][1]

gmsh.model.geo.synchronize()
```

### 3. Build the curve → lateral surface map

```python
orig_curves = {abs(c) for _, c in gmsh.model.getBoundary([(2, surf)], oriented=True)}

curve_to_lat = {}
for _, s in out[2:]:
    for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
        ctag = abs(c)
        if ctag in orig_curves:
            curve_to_lat[ctag] = s
            break
```

### 4. BoundaryDoubleCorner field

```python
N_COLS    = 20        # ~15-25; more → smoother fan near tip
COL_WIDTH = h1        # or a fixed fraction of lc_wall

bdc = gmsh.model.mesh.field.add("BoundaryDoubleCorner")
gmsh.model.mesh.field.setNumbers(bdc, "CurvesList",      [CURVE_BODY])
gmsh.model.mesh.field.setNumbers(bdc, "NosePoint",       [x_nose, 0.0])
gmsh.model.mesh.field.setNumbers(bdc, "TailPoint",       [x_tail, 0.0])
gmsh.model.mesh.field.setNumber (bdc, "Size",            h1)
gmsh.model.mesh.field.setNumber (bdc, "Ratio",           r_bl)
gmsh.model.mesh.field.setNumber (bdc, "NbLayers",        n_bl)
gmsh.model.mesh.field.setNumber (bdc, "NbCornerColumns", N_COLS)
gmsh.model.mesh.field.setNumber (bdc, "MaxColumnWidth",  lc_wall)
gmsh.model.mesh.field.setNumber (bdc, "ColWidth",        COL_WIDTH)

gmsh.model.mesh.field.setAsBackgroundMesh(bdc)
gmsh.option.setNumber("Mesh.BoundaryCornerField",     bdc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_far)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h1 * 0.5)
```

You can still add a `Box` wake refinement field and wrap everything in a `Min` field —
just set **that** `Min` field as background mesh and keep `Mesh.BoundaryCornerField = bdc`.

```python
# Optional: wake refinement
fid_wake = gmsh.model.mesh.field.add("Box")
# ... configure ...

fid_min = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(fid_min, "FieldsList", [bdc, fid_wake])
gmsh.model.mesh.field.setAsBackgroundMesh(fid_min)
gmsh.option.setNumber("Mesh.BoundaryCornerField", bdc)  # still the BDC tag
```

### 5. Generate mesh and assign physical groups

```python
gmsh.model.mesh.generate(3)   # no generate(2)/recombine() before this

# Physical groups via curve_to_lat
gmsh.model.addPhysicalGroup(3, [vol],                                         name="fluid")
gmsh.model.addPhysicalGroup(2, [surf],                                        name="frontAndBack")
gmsh.model.addPhysicalGroup(2, [s_top],                                       name="frontAndBack")
gmsh.model.addPhysicalGroup(2, [curve_to_lat[CURVE_BODY]],                    name="body")
gmsh.model.addPhysicalGroup(2, [curve_to_lat[L_RIGHT]],                       name="outlet")
gmsh.model.addPhysicalGroup(2, [curve_to_lat[L_TOP]],                         name="farfield")
gmsh.model.addPhysicalGroup(2, [curve_to_lat[L_LEFT]],                        name="inlet")
gmsh.model.addPhysicalGroup(2, [curve_to_lat[L_AX_R], curve_to_lat[L_AX_L]], name="axis")
```

Note: `frontAndBack` collects both the original 2D face (`surf`) and the extruded top face
(`s_top`) — pass both tags in a single `addPhysicalGroup` call:
```python
gmsh.model.addPhysicalGroup(2, [surf, s_top], name="frontAndBack")
```

---

## `patch_fix_step` — what to change

gmshToFoam puts structured axis-column faces in `defaultFaces`. Rename it and fix types:

```python
def patch_fix_step(case_dir):
    boundary_file = Path(case_dir) / "constant" / "polyMesh" / "boundary"
    txt = boundary_file.read_text()

    # Rename gmshToFoam's default symmetry plane
    txt = re.sub(r'\bdefaultFaces\b', 'symAxis', txt)

    # Set correct types
    for patch, ptype in [
        ("frontAndBack", "empty"),
        ("axis",         "symmetry"),
        ("symAxis",      "symmetry"),
        ("body",         "wall"),
    ]:
        txt = re.sub(
            rf'({re.escape(patch)}\s*{{[^}}]*?type\s+)\w+\s*;',
            rf'\g<1>{ptype};', txt, flags=re.DOTALL)
        txt = re.sub(
            rf'({re.escape(patch)}\s*{{[^}}]*?physicalType\s+)\w+\s*;',
            rf'\g<1>{ptype};', txt, flags=re.DOTALL)

    boundary_file.write_text(txt)
```

---

## `controlDict` — writePrecision

Set `writePrecision 15` to avoid false non-orthogonality spikes from float truncation.
This is a solved problem (see nonortho_fix memory entry); just ensure the skeleton writer
in `convert_step` includes it:

```
writePrecision  15;
```

---

## What to delete from `build_mesh`

- `r_cap_threshold`, `nose_split_idx`, `tail_split_idx`, `x_nose_inner`, `x_tail_inner`,
  `r_split_nose`, `r_split_tail`, `lc_cap` — cap logic entirely removed.
- `CURVE_NOSE_CAP`, `CURVE_BODY_MAIN`, `CURVE_TAIL_CAP` → replaced by `CURVE_BODY`.
- `A_nose`, `A_tail` points and the four near/far axis lines.
- `N_body_nose`, `N_body_tail`, `setTransfiniteCurve` calls.
- `fid_bl` (`BoundaryLayer` field) and `fid_dist` / `fid_thr` (`Distance`/`Threshold`).
- `generate(2)`, `recombine()`, `optimize("Laplace2D")`, `optimize("Relocate2D")`.
- The bounding-box physical group loop (replaced by `curve_to_lat`).
- `_fix_inverted_hexes_in_msh` — pre-mesh structured quads prevent inversion.

---

## Reference

Working end-to-end example: `testing/bigdoublecorner.py`.
Field options: `local/CLAUDE.md` § Options.
Non-ortho precision fix: memory entry `nonortho_fix`.
