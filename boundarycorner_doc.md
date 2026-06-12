# BoundaryCorner — Agent Setup Guide

## What it solves

Standard Gmsh boundary-layer generation fails when a body surface meets the symmetry axis at 90° (e.g. a nose or tail stagnation point).  At that corner, averaging the wall-normal with the axis-normal produces degenerate or inverted cells.

`BoundaryCornerField` replaces that zone with a **structured fan of quads** that sweeps the 90° rotation of the outward normal from the body shoulder down to the axis tip, producing clean polar cells with no degeneracy.

## Geometry assumptions

- **x** = axial direction (symmetry axis = y = 0).
- **y** = radial direction (y ≥ 0, upper half-plane only).
- The body profile is a curve from some point at y > 0 (shoulder) to the **AxisPoint** at y = 0 (nose or tail), where the tangent is vertical — perpendicular to x.
- The fluid domain is the region **outside** the body (the field only operates in the half-plane above the axis).

## Wedge meshing does NOT work

Revolving the 2D mesh into a 3D wedge (gmsh extrude + rotate) crashes or produces invalid cells at the axis.  **Do not attempt it.**

The correct approach for axisymmetric CFD:

1. Generate a 2D planar mesh (this fork).
2. Extrude it 1 cell deep along z (`DZ = 1.0`) to get a flat volume.
3. Set `front`/`back` patches to `empty` in OpenFOAM.
4. Handle the 2πr axisymmetry factor in the solver / post-processing.

## API — BoundaryCornerField options

| Option | Type | Description |
|---|---|---|
| `CurvesList` | int list | Tags of the GEdge(s) forming the body profile in the corner zone |
| `AxisPoint` | float list | `[x_a, 0.0]` — the tip point where the curve meets the axis |
| `Size` | float | Height of the first boundary-layer row (`h1`) |
| `Ratio` | float | BL growth ratio (e.g. 1.2) |
| `NbLayers` | int | Number of BL rows in the normal direction |
| `NbCornerColumns` | int | Number of fan columns in the tangential direction near the corner |
| `MaxColumnWidth` | float | Maximum tangential column width at the StartPoint (typically = far-field `lc` on the body) |

`Omega` (default 1.0) and `ColWidth` are additional optional parameters.

## Required global option

```python
gmsh.option.setNumber("Mesh.BoundaryCornerField", <field_tag>)
```

This tells the mesher which field drives the structured fan insertion.  Pass the tag of **one** of the BoundaryCorner fields (if you have both nose and tail, pass either — both are active via the Min field).

## Minimal example (nose only)

```python
bc = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc, "CurvesList",      [arc_front])
gmsh.model.mesh.field.setNumbers(bc, "AxisPoint",       [A, 0.0])
gmsh.model.mesh.field.setNumber (bc, "Size",            0.012)
gmsh.model.mesh.field.setNumber (bc, "Ratio",           1.20)
gmsh.model.mesh.field.setNumber (bc, "NbLayers",        6)
gmsh.model.mesh.field.setNumber (bc, "NbCornerColumns", 8)
gmsh.model.mesh.field.setNumber (bc, "MaxColumnWidth",  0.06)

gmsh.model.mesh.field.setAsBackgroundMesh(bc)
gmsh.option.setNumber("Mesh.BoundaryCornerField", bc)
```

## Two-corner example (nose + tail)

```python
def add_bc(curves, axis_pt):
    f = gmsh.model.mesh.field.add("BoundaryCorner")
    gmsh.model.mesh.field.setNumbers(f, "CurvesList",      curves)
    gmsh.model.mesh.field.setNumbers(f, "AxisPoint",       axis_pt)
    gmsh.model.mesh.field.setNumber (f, "Size",            H1)
    gmsh.model.mesh.field.setNumber (f, "Ratio",           RATIO)
    gmsh.model.mesh.field.setNumber (f, "NbLayers",        N_LAY)
    gmsh.model.mesh.field.setNumber (f, "NbCornerColumns", N_COLS)
    gmsh.model.mesh.field.setNumber (f, "MaxColumnWidth",  W0_MAX)
    return f

bc_front = add_bc([arc_front], [A,  0.0])
bc_back  = add_bc([arc_back],  [-A, 0.0])

min_f = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_f, "FieldsList", [bc_front, bc_back])
gmsh.model.mesh.field.setAsBackgroundMesh(min_f)

gmsh.option.setNumber("Mesh.BoundaryCornerField", bc_front)  # either tag
```

## OpenFOAM integration

After meshing, extrude 1 cell along z and export as MSH 2.2:

```python
out = gmsh.model.geo.extrude([(2, sf)], 0, 0, DZ, numElements=[1], recombine=True)
# ...
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("mesh.msh")
```

Convert and fix patches:

```bash
gmshToFoam mesh.msh
```

The BoundaryCorner axis-column cells at y = 0 are placed in `defaultFaces` by gmshToFoam.  Rename it to a distinct symmetry patch (e.g. `symAxis`) and set types:

| Patch | OpenFOAM type |
|---|---|
| `front`, `back` | `empty` |
| `body` | `wall` |
| `symmetry`, `symAxis` | `symmetry` |
| `inlet`, `outlet`, `top` | `patch` |

See `testing/2dplanar.py → fix_of_boundary()` for a regex-based approach.

## Full working example

`testing/2dplanar.py` — half-ellipse body (A=2, B=1.2) with nose and tail corners, k-ω SST, SIMPLE solver.

```bash
# source OpenFOAM first, then:
python testing/2dplanar.py          # headless
python testing/2dplanar.py --gui    # open mesh in gmsh GUI first
```

Output lands in `testing/ellipse2dplanar_of/`.
