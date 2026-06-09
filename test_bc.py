"""
test_bc.py — BoundaryCorner on a half-ellipse nose profile.

Geometry (external axisymmetric CFD, upper-half 2D cross-section):
  - Ellipse semi-axes a=2.0 (axial), b=1.2 (radial)
  - Body surface = quarter ellipse arc from (0, b) down to the nose (a, 0)
  - Fluid domain = region outside the ellipse, bounded by a rectangular far-field
  - Symmetry axis = y=0 (horizontal)
  - AxisPoint = (a, 0) — nose tip, where wall meets axis at 90°
  - StartPoint = point on the arc at theta=55°, beginning of BC zone

Run:  python test_bc.py           # headless
      python test_bc.py --gui     # open result in gmsh GUI
"""

import sys, os, math, subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "api"))
os.environ["GMSH_LIB"] = os.path.join(os.path.dirname(__file__), "build", "libgmsh.so")

import gmsh

# ── Geometry parameters ───────────────────────────────────────────────────────
A = 2.0     # ellipse semi-major axis (axial, x direction)
B = 1.2     # ellipse semi-minor axis (radial, y direction)
X_FAR = 5.0
Y_FAR = 3.5

# ── BoundaryCorner parameters ─────────────────────────────────────────────────
THETA_START = 55 * math.pi / 180     # profile angle where BC zone starts
X_START = A * math.cos(THETA_START)  # ≈ 1.147
Y_START = B * math.sin(THETA_START)  # ≈ 0.983

H1     = 0.012   # first BL layer height (normal direction)
RATIO  = 1.20    # geometric ratio, normal direction
N_LAY  = 6       # number of BL layers
N_COL  = 7       # number of BC columns along the arc
DELTA1 = 0.08    # arc length of last column at nose (< first column → compression)

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.25
LC_BODY = 0.06
LC_NOSE = 0.03

gui = "--gui" in sys.argv
gmsh.initialize(["gmsh"] if gui else ["gmsh", "-nopopup"])
gmsh.model.add("test_bc_ellipse")

# ── Points ────────────────────────────────────────────────────────────────────
# p_origin is the ellipse centre — used only as reference for addEllipseArc
p_origin = gmsh.model.geo.addPoint(0,     0,     0, LC_BODY)
p_nose   = gmsh.model.geo.addPoint(A,     0,     0, LC_NOSE)  # AxisPoint
p_top    = gmsh.model.geo.addPoint(0,     B,     0, LC_BODY)  # y-axis / ellipse top
p_ax_r   = gmsh.model.geo.addPoint(X_FAR, 0,    0, LC_FAR)   # far-right on axis
p_far_tr = gmsh.model.geo.addPoint(X_FAR, Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint(0,    Y_FAR,  0, LC_FAR)

# ── Curves ────────────────────────────────────────────────────────────────────
# Ellipse arc: clockwise from p_top (0,B) to p_nose (A,0)
# (addEllipseArc: start, centre, major-axis point, end)
arc   = gmsh.model.geo.addEllipseArc(p_top, p_origin, p_nose, p_nose)

l_l   = gmsh.model.geo.addLine(p_top,    p_far_tl)   # left boundary (upward)
l_top = gmsh.model.geo.addLine(p_far_tl, p_far_tr)   # top boundary (rightward)
l_r   = gmsh.model.geo.addLine(p_far_tr, p_ax_r)     # right boundary (downward)
l_ax  = gmsh.model.geo.addLine(p_ax_r,   p_nose)     # axis y=0 (leftward, to nose)

# Fluid domain: CCW loop  (arc is CW, so we reverse it with −arc)
# Path: p_nose →[-arc]→ p_top →[l_l]→ p_far_tl →[l_top]→ p_far_tr →[l_r]→ p_ax_r →[l_ax]→ p_nose
cl = gmsh.model.geo.addCurveLoop([-arc, l_l, l_top, l_r, l_ax])
sf = gmsh.model.geo.addPlaneSurface([cl])
gmsh.model.geo.synchronize()

# ── BoundaryCorner field ───────────────────────────────────────────────────────
bc = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc, "CurvesList",      [arc])
gmsh.model.mesh.field.setNumbers(bc, "AxisPoint",       [A,       0.0    ])
gmsh.model.mesh.field.setNumbers(bc, "StartPoint",      [X_START, Y_START])
gmsh.model.mesh.field.setNumber (bc, "Size",            H1)
gmsh.model.mesh.field.setNumber (bc, "Ratio",           RATIO)
gmsh.model.mesh.field.setNumber (bc, "NbLayers",        N_LAY)
gmsh.model.mesh.field.setNumber (bc, "NbCornerColumns", N_COL)
gmsh.model.mesh.field.setNumber (bc, "Delta1",          DELTA1)
gmsh.model.mesh.field.setAsBackgroundMesh(bc)

gmsh.option.setNumber("Mesh.BoundaryCornerField",        bc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax",    LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin",    H1)

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

expected = N_COL * N_LAY
print(f"\n  expected quads : {expected}  (N_COL={N_COL} × N_LAY={N_LAY})")
status = "✓" if n_quads == expected else "✗"
print(f"  {status} quad count : {n_quads}")

out = "/tmp/test_bc_ellipse.msh"
gmsh.write(out)
print(f"\nwritten: {out}")

if gui:
    subprocess.run(["/usr/bin/gmsh", out])

gmsh.finalize()
