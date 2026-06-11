"""
test_bc.py — BoundaryCorner on a half-ellipse nose profile,
             with an adjacent BoundaryLayer field on the upper arc.

Geometry (external axisymmetric CFD, upper-half 2D cross-section):
  - Ellipse semi-axes a=2.0 (axial), b=1.2 (radial)
  - Body surface = quarter ellipse arc from (0, b) down to the nose (a, 0)
  - Fluid domain = region outside the ellipse, bounded by a rectangular far-field
  - Symmetry axis = y=0 (horizontal)
  - AxisPoint = (a, 0) — nose tip, where wall meets axis at 90°
  - StartPoint = point on the arc at theta=55°, junction between BL and BC zones

Arc is split at StartPoint:
  arc_bl  :  p_top  → p_start  (theta 90°→55°) — BoundaryLayer field
  arc_bc  :  p_start → p_nose  (theta 55°→ 0°) — BoundaryCorner field

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
THETA_START = 55 * math.pi / 180     # profile angle where BC zone starts
X_START = A * math.cos(THETA_START)  # ≈ 1.147
Y_START = B * math.sin(THETA_START)  # ≈ 0.983

H1     = 0.012   # first BL layer height (normal direction)
RATIO  = 1.20    # geometric ratio, normal direction
N_LAY  = 6       # number of BL layers
W0_MAX = LC_BODY  # max BC column arc-length (at S); N columns derived automatically

gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("test_bc_ellipse")

# ── Points ────────────────────────────────────────────────────────────────────
# p_origin is the ellipse centre — used only as reference for addEllipseArc
p_origin = gmsh.model.geo.addPoint(0,       0,      0, LC_BODY)
p_start  = gmsh.model.geo.addPoint(X_START, Y_START,0, LC_BODY)  # BL/BC junction
p_nose   = gmsh.model.geo.addPoint(A,       0,      0, LC_NOSE)  # AxisPoint
p_top    = gmsh.model.geo.addPoint(0,       B,      0, LC_BODY)  # y-axis / ellipse top
p_ax_r   = gmsh.model.geo.addPoint(X_FAR,   0,      0, LC_FAR)   # far-right on axis
p_far_tr = gmsh.model.geo.addPoint(X_FAR,   Y_FAR,  0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint(0,       Y_FAR,  0, LC_FAR)

# ── Curves ────────────────────────────────────────────────────────────────────
# Ellipse arcs: both CW (addEllipseArc: start, centre, major-axis-pt, end)
arc_bl = gmsh.model.geo.addEllipseArc(p_top,   p_origin, p_nose, p_start)  # upper (BL)
arc_bc = gmsh.model.geo.addEllipseArc(p_start, p_origin, p_nose, p_nose)   # lower (BC)

l_l   = gmsh.model.geo.addLine(p_top,    p_far_tl)   # left boundary (upward)
l_top = gmsh.model.geo.addLine(p_far_tl, p_far_tr)   # top boundary (rightward)
l_r   = gmsh.model.geo.addLine(p_far_tr, p_ax_r)     # right boundary (downward)
l_ax  = gmsh.model.geo.addLine(p_ax_r,   p_nose)     # axis y=0 (leftward, to nose)

# Fluid domain: CCW loop — both arcs are CW so reversed with minus signs
# Path: p_nose →[-arc_bc]→ p_start →[-arc_bl]→ p_top →[l_l]→ ... →[l_ax]→ p_nose
cl = gmsh.model.geo.addCurveLoop([-arc_bc, -arc_bl, l_l, l_top, l_r, l_ax])
sf = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize()

# Total BL thickness (geometric series: H1 + H1*r + ... + H1*r^(N-1))
H_TOTAL = H1 * (RATIO**N_LAY - 1.0) / (RATIO - 1.0)

# ── BoundaryLayer field on upper arc (p_top → p_start) ───────────────────────
bl = gmsh.model.mesh.field.add("BoundaryLayer")
gmsh.model.mesh.field.setNumbers(bl, "EdgesList", [arc_bl])
gmsh.model.mesh.field.setNumber (bl, "hwall_n",   H1)
gmsh.model.mesh.field.setNumber (bl, "ratio",     RATIO)
gmsh.model.mesh.field.setNumber (bl, "thickness", H_TOTAL)
gmsh.model.mesh.field.setNumber (bl, "Quads",     1)

# ── BoundaryCorner field on lower arc (p_start → p_nose) ─────────────────────
bc = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc, "CurvesList",      [arc_bc])
gmsh.model.mesh.field.setNumbers(bc, "AxisPoint",       [A,       0.0    ])
gmsh.model.mesh.field.setNumbers(bc, "StartPoint",      [X_START, Y_START])
gmsh.model.mesh.field.setNumber (bc, "Size",            H1)
gmsh.model.mesh.field.setNumber (bc, "Ratio",           RATIO)
gmsh.model.mesh.field.setNumber (bc, "NbLayers",        N_LAY)
gmsh.model.mesh.field.setNumber (bc, "MaxColumnWidth",  W0_MAX)

# setAsBoundaryLayer registers bl in _boundaryLayer_fields, which causes
# buildAdditionalPoints2D → modifyInitialMeshForBoundaryLayers to run during
# meshGFace and generate structured BL quads in the upper arc.
# setAsBackgroundMesh on bc provides size guidance for the rest of the domain.
gmsh.model.mesh.field.setAsBoundaryLayer(bl)
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

print(f"\n  BL layers : {N_LAY},  W0_MAX={W0_MAX}  (N columns derived by field)")
print(f"  total quads : {n_quads}  (BL quads counted separately not yet available here)")

out = "/tmp/test_bc_ellipse.msh"
gmsh.write(out)
print(f"\nwritten: {out}")

if gui:
    # Prefer the dev build; fall back to the system gmsh if the dev build
    # has no FLTK GUI support (libfltk1.3-dev not installed).
    dev_bin = os.path.join(_root, "build", "gmsh")
    import shutil
    gmsh_bin = dev_bin if shutil.which(dev_bin) else "gmsh"
    # Quick check: if the dev binary has no FLTK, it exits immediately.
    # Detect by checking whether FLTK is listed in its cmake features.
    import subprocess as _sp
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
