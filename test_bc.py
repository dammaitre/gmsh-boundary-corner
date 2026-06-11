"""
test_bc.py — BoundaryCorner covering the full half-ellipse nose profile.

Geometry (external axisymmetric CFD, upper-half 2D cross-section):
  - Ellipse semi-axes a=2.0 (axial), b=1.2 (radial)
  - Body surface = quarter ellipse arc from (0, b) down to the nose (a, 0)
  - Fluid domain = region outside the ellipse, bounded by a rectangular far-field
  - Symmetry axis = y=0 (horizontal)
  - AxisPoint = (a, 0) — nose tip, where wall meets axis at 90°

The BoundaryCorner field covers the entire profile arc:
  - NbCornerColumns columns near the nose use geometric arc-length compression
    (decreasing column width toward AxisPoint, from MaxColumnWidth to ColWidth/Size)
  - The remaining profile is covered with constant-width columns of MaxColumnWidth

Run:  python test_bc.py           # headless
      python test_bc.py --gui     # open result in gmsh GUI
"""

import sys, os, math, subprocess

_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry parameters ───────────────────────────────────────────────────────
A = 2.0     # ellipse semi-major axis (axial, x direction)
B = 1.2     # ellipse semi-minor axis (radial, y direction)
X_FAR = 5.0
Y_FAR = 3.5

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.25
LC_BODY = 0.06
LC_NOSE = 0.03

# ── BoundaryCorner parameters ─────────────────────────────────────────────────
H1       = 0.012   # first normal layer height
RATIO    = 1.20    # geometric ratio, normal direction
N_LAY    = 6       # number of normal layers
N_COLS   = 8       # compressed columns near AxisPoint
W0_MAX   = LC_BODY  # column arc-length in the constant-width section

gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("test_bc_ellipse")

# ── Points ────────────────────────────────────────────────────────────────────
p_origin = gmsh.model.geo.addPoint(0,    0,    0, LC_BODY)
p_nose   = gmsh.model.geo.addPoint(A,    0,    0, LC_NOSE)   # AxisPoint
p_top    = gmsh.model.geo.addPoint(0,    B,    0, LC_BODY)   # y-axis / ellipse top
p_ax_r   = gmsh.model.geo.addPoint(X_FAR, 0,  0, LC_FAR)
p_far_tr = gmsh.model.geo.addPoint(X_FAR, Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint(0,    Y_FAR,  0, LC_FAR)

# ── Curves ────────────────────────────────────────────────────────────────────
# Full quarter-ellipse profile from p_top to p_nose (CW)
arc_profile = gmsh.model.geo.addEllipseArc(p_top, p_origin, p_nose, p_nose)

l_l   = gmsh.model.geo.addLine(p_top,    p_far_tl)
l_top = gmsh.model.geo.addLine(p_far_tl, p_far_tr)
l_r   = gmsh.model.geo.addLine(p_far_tr, p_ax_r)
l_ax  = gmsh.model.geo.addLine(p_ax_r,   p_nose)

# Fluid domain: CCW loop; arc is CW so reversed with minus sign
cl = gmsh.model.geo.addCurveLoop([-arc_profile, l_l, l_top, l_r, l_ax])
sf = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize()

# ── BoundaryCorner field on full profile (p_top → p_nose) ────────────────────
bc = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc, "CurvesList",       [arc_profile])
gmsh.model.mesh.field.setNumbers(bc, "AxisPoint",        [A, 0.0])
gmsh.model.mesh.field.setNumber (bc, "Size",             H1)
gmsh.model.mesh.field.setNumber (bc, "Ratio",            RATIO)
gmsh.model.mesh.field.setNumber (bc, "NbLayers",         N_LAY)
gmsh.model.mesh.field.setNumber (bc, "NbCornerColumns",  N_COLS)
gmsh.model.mesh.field.setNumber (bc, "MaxColumnWidth",   W0_MAX)

gmsh.model.mesh.field.setAsBackgroundMesh(bc)

gmsh.option.setNumber("Mesh.BoundaryCornerField",        bc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax",    LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin",    H1 * 0.5)

# ── Generate ──────────────────────────────────────────────────────────────────
gmsh.model.mesh.generate(2)

# ── Results ───────────────────────────────────────────────────────────────────
elem_types, elem_tags, _ = gmsh.model.mesh.getElements(2)
print("\n=== Results ===")
type_names = {2: "triangle", 3: "quadrangle"}
n_quads = n_tris = 0
for t, tags in zip(elem_types, elem_tags):
    if not len(tags):
        continue
    name = type_names.get(t, f"type{t}")
    print(f"  {name}: {len(tags)}")
    if t == 3: n_quads = len(tags)
    if t == 2: n_tris  = len(tags)

print(f"\n  NbCornerColumns={N_COLS},  NbLayers={N_LAY},  W0_MAX={W0_MAX}")
print(f"  total quads : {n_quads}")

out = "/tmp/test_bc_ellipse.msh"
gmsh.write(out)
print(f"\nwritten: {out}")

if gui:
    dev_bin = os.path.join(_root, "build", "gmsh")
    import shutil
    gmsh_bin = dev_bin if shutil.which(dev_bin) else "gmsh"
    cache = os.path.join(_root, "build", "CMakeCache.txt")
    has_fltk = False
    if os.path.exists(cache):
        with open(cache) as _f:
            for _line in _f:
                if "FLTK_BASE_LIBRARY_RELEASE:FILEPATH=" in _line and "NOTFOUND" not in _line:
                    has_fltk = True
                    break
    if not has_fltk:
        gmsh_bin = shutil.which("gmsh") or dev_bin
    subprocess.run([gmsh_bin, out])

gmsh.finalize()
