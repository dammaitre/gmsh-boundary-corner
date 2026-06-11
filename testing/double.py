"""
test_2bc.py — BoundaryCorner on a full half-ellipse (two axis intersections).

Geometry (external axisymmetric CFD, upper-half 2D cross-section):
  - Ellipse semi-axes a=2.0 (axial), b=1.2 (radial)
  - Body surface = FULL upper half-ellipse from tail (-a,0) through top (0,b) to nose (a,0)
  - Two axis intersection points: nose (a,0) and tail (-a,0)
  - Fluid domain = region outside the ellipse, bounded by a rectangular far-field

Strategy — two separate arcs, one BoundaryCorner field each:
  arc_front: (0,b) → (a,0)    covered by bc_front  (AxisPoint = nose)
  arc_back:  (-a,0) → (0,b)   covered by bc_back   (AxisPoint = tail)

modifyInitialMeshForBoundaryCorners loops over ALL BoundaryCornerField instances, so
both fields are processed automatically in a single mesh generation call.

Run:  python test_2bc.py           # headless
      python test_2bc.py --gui     # open result in gmsh GUI
"""

import sys, os, math, subprocess

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry parameters ───────────────────────────────────────────────────────
A = 2.0     # ellipse semi-major axis (axial, x direction)
B = 1.2     # ellipse semi-minor axis (radial, y direction)
X_FAR = 6.0
Y_FAR = 3.5

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.30
LC_BODY = 0.06
LC_NOSE = 0.03

# ── BoundaryCorner parameters ─────────────────────────────────────────────────
H1       = 0.012   # first normal layer height
RATIO    = 1.20    # geometric ratio, normal direction
N_LAY    = 6       # number of normal layers
N_COLS   = 8       # compressed columns near each AxisPoint
W0_MAX   = LC_BODY  # column arc-length in the constant-width section

gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("test_2bc_ellipse")

# ── Points ────────────────────────────────────────────────────────────────────
p_origin  = gmsh.model.geo.addPoint(  0,      0,     0, LC_BODY)
p_nose    = gmsh.model.geo.addPoint(  A,      0,     0, LC_NOSE)   # right AxisPoint
p_tail    = gmsh.model.geo.addPoint( -A,      0,     0, LC_NOSE)   # left AxisPoint
p_top     = gmsh.model.geo.addPoint(  0,      B,     0, LC_BODY)   # top of ellipse
p_ax_r    = gmsh.model.geo.addPoint(  X_FAR,  0,     0, LC_FAR)
p_ax_l    = gmsh.model.geo.addPoint( -X_FAR,  0,     0, LC_FAR)
p_far_tr  = gmsh.model.geo.addPoint(  X_FAR,  Y_FAR, 0, LC_FAR)
p_far_tl  = gmsh.model.geo.addPoint( -X_FAR,  Y_FAR, 0, LC_FAR)

# ── Curves ────────────────────────────────────────────────────────────────────
# Front quarter-ellipse: (0,B) → (A,0)  [CW, same as single-corner test]
arc_front = gmsh.model.geo.addEllipseArc(p_top,  p_origin, p_nose, p_nose)
# Back quarter-ellipse:  (-A,0) → (0,B)  [CCW]
arc_back  = gmsh.model.geo.addEllipseArc(p_tail, p_origin, p_nose, p_top)

# Axis edges defined to ARRIVE AT the axis intersection point.
# This matches the buildForFace convention: axisEdge->getEndVertex()==axisGV.
l_ax_r  = gmsh.model.geo.addLine(p_ax_r,   p_nose)   # downstream axis → nose
l_right = gmsh.model.geo.addLine(p_ax_r,   p_far_tr)
l_top   = gmsh.model.geo.addLine(p_far_tr, p_far_tl)
l_left  = gmsh.model.geo.addLine(p_far_tl, p_ax_l)
l_ax_l  = gmsh.model.geo.addLine(p_ax_l,   p_tail)   # upstream axis → tail

# Fluid domain loop (CCW):
#   p_nose → p_ax_r → p_far_tr → p_far_tl → p_ax_l → p_tail → p_top → p_nose
# -l_ax_r : traverse l_ax_r (p_ax_r→p_nose) reversed, i.e. p_nose→p_ax_r
# +l_ax_l : traverse l_ax_l (p_ax_l→p_tail) forward,  i.e. p_ax_l→p_tail
cl = gmsh.model.geo.addCurveLoop(
    [-l_ax_r, l_right, l_top, l_left, l_ax_l, arc_back, arc_front])
sf = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize()

# ── BoundaryCorner field — FRONT (nose, right axis intersection) ──────────────
bc_front = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc_front, "CurvesList",       [arc_front])
gmsh.model.mesh.field.setNumbers(bc_front, "AxisPoint",        [A, 0.0])
gmsh.model.mesh.field.setNumber (bc_front, "Size",             H1)
gmsh.model.mesh.field.setNumber (bc_front, "Ratio",            RATIO)
gmsh.model.mesh.field.setNumber (bc_front, "NbLayers",         N_LAY)
gmsh.model.mesh.field.setNumber (bc_front, "NbCornerColumns",  N_COLS)
gmsh.model.mesh.field.setNumber (bc_front, "MaxColumnWidth",   W0_MAX)

# ── BoundaryCorner field — BACK (tail, left axis intersection) ────────────────
bc_back = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc_back, "CurvesList",       [arc_back])
gmsh.model.mesh.field.setNumbers(bc_back, "AxisPoint",        [-A, 0.0])
gmsh.model.mesh.field.setNumber (bc_back, "Size",             H1)
gmsh.model.mesh.field.setNumber (bc_back, "Ratio",            RATIO)
gmsh.model.mesh.field.setNumber (bc_back, "NbLayers",         N_LAY)
gmsh.model.mesh.field.setNumber (bc_back, "NbCornerColumns",  N_COLS)
gmsh.model.mesh.field.setNumber (bc_back, "MaxColumnWidth",   W0_MAX)

# Min field combines both BC size functions for background mesh guidance.
# Each BC field returns a small size near its AxisPoint and 1e22 elsewhere.
min_f = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_f, "FieldsList", [bc_front, bc_back])
gmsh.model.mesh.field.setAsBackgroundMesh(min_f)

# modifyInitialMeshForBoundaryCorners iterates over ALL BoundaryCornerField
# instances; setting BoundaryCornerField to any one of them is sufficient to
# signal that the BC pre-mesh step should run (the field itself is found by type).
gmsh.option.setNumber("Mesh.BoundaryCornerField",     bc_front)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

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

print(f"\n  NbCornerColumns={N_COLS} × 2,  NbLayers={N_LAY},  W0_MAX={W0_MAX}")
print(f"  total quads : {n_quads}")

out = "/tmp/test_2bc_ellipse.msh"
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
