"""
testing_bc_bug.py — Reproduce the negative-volume / open-cell bug we see in
phase_mesh.py's BoundaryCorner implementation, using a simple half-ellipse body
(same geometry as ../gmsh/testing/2dplanar.py) so we can diff against the working
reference.

The reference (2dplanar.py) uses:
  - EllipseArc body curves
  - CW curve loop  ([-l_ax_r, l_right, l_top, l_left, l_ax_l, arc_back, arc_front])
  - Axis lines that go from far-field to body tip (p_ax_r→p_nose, p_ax_l→p_tail)

Our pipeline (phase_mesh.py) uses:
  - Spline body curves split into nose/main/tail
  - CCW curve loop
  - Axis lines that go from body tip to far-field (ptag_nose→P_outlet_axis etc.)

We reproduce the pipeline approach here with the ellipse geometry so we can run
checkMesh and compare with the working reference.

Usage:
    source .venv/bin/activate
    python testing_bc_bug.py          # our pipeline approach → should show bug
    python testing_bc_bug.py --ref    # 2dplanar.py reference approach → should be clean
"""

import sys, os, re, subprocess, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("GMSH_LIB",
    "/mnt/DAMIEN/DEV/gmsh/build/libgmsh.so")

import gmsh

REF_MODE = "--ref" in sys.argv

# ── Geometry (same as 2dplanar.py) ────────────────────────────────────────────
A      = 2.0    # ellipse semi-major (axial)
B      = 1.2    # ellipse semi-minor (radial)
X_FAR  = 6.0
Y_FAR  = 3.5
DZ     = 1.0

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.30
LC_BODY = 0.06
LC_NOSE = 0.03

# ── BoundaryCorner parameters (same as reference) ─────────────────────────────
H1     = 0.012
RATIO  = 1.20
N_LAY  = 6
N_COLS = 8
W0_MAX = LC_BODY

# ── Output ────────────────────────────────────────────────────────────────────
label    = "ref" if REF_MODE else "ours"
CASE_DIR = os.path.join(os.path.dirname(__file__), f"bc_bug_{label}")

print(f"\n{'='*60}")
print(f"  Mode: {'REFERENCE (2dplanar.py approach)' if REF_MODE else 'OURS (phase_mesh.py approach)'}")
print(f"  Output: {CASE_DIR}")
print(f"{'='*60}\n")

if os.path.isdir(CASE_DIR):
    shutil.rmtree(CASE_DIR)
os.makedirs(CASE_DIR)

gmsh.initialize()
gmsh.option.setNumber("General.Verbosity",  1)
gmsh.option.setNumber("General.ExpertMode", 1)
gmsh.model.add("bc_bug")

# ── Build a 200-point cosine-spaced spline profile of the half-ellipse ────────
# x ∈ [-A, A],  r(x) = B * sqrt(1 - (x/A)^2)
# nose at (A, 0), tail at (-A, 0)  — same convention as 2dplanar.py
N_PTS = 200
s     = (1 - np.cos(np.linspace(0, np.pi, N_PTS))) / 2   # 0…1, cosine-spaced
xs    = A - 2*A*s                                          # A … -A
rs    = B * np.sqrt(np.maximum(1 - (xs/A)**2, 0))

pts_xy = np.column_stack([xs, rs])   # nose=(A,0) … tail=(-A,0)

# ── Profile split at the body shoulder (max radius) ───────────────────────────
top_idx = int(np.argmax(rs))   # index of max r — shoulder at (0, B)

x_nose = float(pts_xy[0,  0])   #  A
x_tail = float(pts_xy[-1, 0])   # -A

print(f"  top_idx={top_idx}  shoulder r={pts_xy[top_idx, 1]:.4f}")

# ── Add all profile points ─────────────────────────────────────────────────────
prof_ptags = [gmsh.model.geo.addPoint(float(x), float(r), 0.0, LC_BODY)
              for x, r in pts_xy]

ptag_nose = prof_ptags[0]
ptag_tail = prof_ptags[-1]

# ── Farfield points ────────────────────────────────────────────────────────────
if REF_MODE:
    # 2dplanar.py layout: axis points on the FAR side of body tips
    p_ax_r   = gmsh.model.geo.addPoint( X_FAR,  0,     0, LC_FAR)
    p_ax_l   = gmsh.model.geo.addPoint(-X_FAR,  0,     0, LC_FAR)
    p_far_tr = gmsh.model.geo.addPoint( X_FAR,  Y_FAR, 0, LC_FAR)
    p_far_tl = gmsh.model.geo.addPoint(-X_FAR,  Y_FAR, 0, LC_FAR)
else:
    # Our layout: axis points are inlet/outlet
    P_inlet_axis   = gmsh.model.geo.addPoint(-X_FAR,  0,     0, LC_FAR)
    P_inlet_top    = gmsh.model.geo.addPoint(-X_FAR,  Y_FAR, 0, LC_FAR)
    P_farfield_top = gmsh.model.geo.addPoint( X_FAR,  Y_FAR, 0, LC_FAR)
    P_outlet_axis  = gmsh.model.geo.addPoint( X_FAR,  0,     0, LC_FAR)

# ── Body spline curves: front half (nose→shoulder) and back half (shoulder→tail)
CURVE_FRONT = gmsh.model.geo.addSpline(prof_ptags[0:top_idx + 1])
CURVE_BACK  = gmsh.model.geo.addSpline(prof_ptags[top_idx:])

# ── Axis and farfield curves ───────────────────────────────────────────────────
if REF_MODE:
    # 2dplanar.py: axis lines go FROM far TO body tip
    l_ax_r  = gmsh.model.geo.addLine(p_ax_r,   ptag_nose)   # right far → nose
    l_ax_l  = gmsh.model.geo.addLine(p_ax_l,   ptag_tail)   # left far  → tail
    l_right = gmsh.model.geo.addLine(p_ax_r,   p_far_tr)
    l_top   = gmsh.model.geo.addLine(p_far_tr, p_far_tl)
    l_left  = gmsh.model.geo.addLine(p_far_tl, p_ax_l)
else:
    # Both axis edges go far→C (axisVert at END, matching reference convention)
    CURVE_AXIS_FWD = gmsh.model.geo.addLine(P_inlet_axis,  ptag_tail)   # (-6,0)→(-2,0)
    CURVE_AXIS_AFT = gmsh.model.geo.addLine(P_outlet_axis, ptag_nose)   # (6,0)→(2,0)
    CURVE_INLET    = gmsh.model.geo.addLine(P_inlet_axis,   P_inlet_top)
    CURVE_FARFIELD = gmsh.model.geo.addLine(P_inlet_top,    P_farfield_top)
    CURVE_OUTLET   = gmsh.model.geo.addLine(P_farfield_top, P_outlet_axis)

# ── Curve loop ─────────────────────────────────────────────────────────────────
if REF_MODE:
    # CW loop matching 2dplanar.py: body traversal goes tail→back→front→nose (reversed)
    cl = gmsh.model.geo.addCurveLoop([
        -l_ax_r,        # ptag_nose → p_ax_r
        l_right,
        l_top,
        l_left,
        l_ax_l,         # p_ax_l → ptag_tail
        -CURVE_BACK,    # ptag_tail → shoulder (reversed: going UP)
        -CURVE_FRONT,   # shoulder → ptag_nose (reversed: going DOWN to nose)
    ])
else:
    # CCW loop: left-axis → tail → (body reversed) → nose → right-axis → farfield
    cl = gmsh.model.geo.addCurveLoop([
        CURVE_AXIS_FWD,
        -CURVE_BACK,
        -CURVE_FRONT,
        -CURVE_AXIS_AFT,
        -CURVE_OUTLET,
        -CURVE_FARFIELD,
        -CURVE_INLET,
    ])

sf  = gmsh.model.geo.addPlaneSurface([cl])

# ── Extrude geometry before meshing ───────────────────────────────────────────
out   = gmsh.model.geo.extrude([(2, sf)], 0, 0, DZ, numElements=[1], recombine=True)
s_top = out[0][1]
vol   = out[1][1]

gmsh.model.geo.synchronize()

# ── Map curves → lateral surfaces ─────────────────────────────────────────────
orig_curves = {abs(c) for _, c in gmsh.model.getBoundary([(2, sf)], oriented=True)}
curve_to_lat = {}
for _, s in out[2:]:
    for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
        ctag = abs(c)
        if ctag in orig_curves and ctag not in curve_to_lat:
            curve_to_lat[ctag] = s
            break

# ── BoundaryCorner fields ──────────────────────────────────────────────────────
def add_bc_field(curves, axis_pt):
    f = gmsh.model.mesh.field.add("BoundaryCorner")
    gmsh.model.mesh.field.setNumbers(f, "CurvesList",      curves)
    gmsh.model.mesh.field.setNumbers(f, "AxisPoint",       axis_pt)
    gmsh.model.mesh.field.setNumber (f, "Size",            H1)
    gmsh.model.mesh.field.setNumber (f, "Ratio",           RATIO)
    gmsh.model.mesh.field.setNumber (f, "NbLayers",        N_LAY)
    gmsh.model.mesh.field.setNumber (f, "NbCornerColumns", N_COLS)
    gmsh.model.mesh.field.setNumber (f, "MaxColumnWidth",  W0_MAX)
    return f

bc_front = add_bc_field([CURVE_FRONT], [x_nose, 0.0])
bc_back  = add_bc_field([CURVE_BACK],  [x_tail, 0.0])

min_f = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_f, "FieldsList", [bc_front, bc_back])
gmsh.model.mesh.field.setAsBackgroundMesh(min_f)

gmsh.option.setNumber("Mesh.BoundaryCornerField",     bc_front)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

# ── Generate 3D mesh ───────────────────────────────────────────────────────────
print("Generating 3D mesh...")
gmsh.model.mesh.generate(3)
gmsh.model.mesh.setOrder(1)
try:
    _, etags, _ = gmsh.model.mesh.getElements(3)
    n3d = sum(len(t) for t in etags)
    print(f"  3D elements generated: {n3d}")
except Exception as e:
    print(f"  [WARN] Could not count 3D elements: {e}")

# ── Physical groups ────────────────────────────────────────────────────────────
def pg(dim, tags, tag, name):
    gmsh.model.addPhysicalGroup(dim, tags, tag)
    gmsh.model.setPhysicalName(dim, tag, name)

body_lats = [curve_to_lat[c] for c in [CURVE_FRONT, CURVE_BACK]
             if c in curve_to_lat]

if REF_MODE:
    axis_lats = [curve_to_lat[c] for c in [l_ax_r, l_ax_l] if c in curve_to_lat]
    inlet_lats  = [curve_to_lat[l_left]]  if l_left  in curve_to_lat else []
    top_lats    = [curve_to_lat[l_top]]   if l_top   in curve_to_lat else []
    outlet_lats = [curve_to_lat[l_right]] if l_right in curve_to_lat else []
else:
    axis_lats   = [curve_to_lat[c] for c in [CURVE_AXIS_FWD, CURVE_AXIS_AFT]
                   if c in curve_to_lat]
    inlet_lats  = [curve_to_lat[CURVE_INLET]]    if CURVE_INLET    in curve_to_lat else []
    top_lats    = [curve_to_lat[CURVE_FARFIELD]]  if CURVE_FARFIELD in curve_to_lat else []
    outlet_lats = [curve_to_lat[CURVE_OUTLET]]    if CURVE_OUTLET   in curve_to_lat else []

pg(3, [vol],           1, "fluid")
pg(2, [sf, s_top],     2, "frontAndBack")
pg(2, body_lats,       3, "body")
pg(2, axis_lats,       4, "axis")
if inlet_lats:  pg(2, inlet_lats,  5, "inlet")
if top_lats:    pg(2, top_lats,    6, "farfield")
if outlet_lats: pg(2, outlet_lats, 7, "outlet")

# ── Write mesh ─────────────────────────────────────────────────────────────────
msh_path = os.path.join(CASE_DIR, "mesh.msh")
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(msh_path)
gmsh.finalize()
print(f"Mesh written: {msh_path}  ({os.path.getsize(msh_path)/1e3:.0f} kB)")

# ── gmshToFoam ────────────────────────────────────────────────────────────────
sys_dir = os.path.join(CASE_DIR, "system")
os.makedirs(sys_dir, exist_ok=True)
with open(os.path.join(sys_dir, "controlDict"), "w") as f:
    f.write("FoamFile{version 2.0;format ascii;class dictionary;object controlDict;}\n"
            "application simpleFoam;\nstartFrom startTime;\nstartTime 0;\n"
            "stopAt endTime;\nendTime 1;\ndeltaT 1;\nwriteControl timeStep;\n"
            "writeInterval 1;\n")

r = subprocess.run(["gmshToFoam", "-case", CASE_DIR, msh_path],
                   capture_output=True, text=True)
if r.returncode != 0:
    print("gmshToFoam FAILED:\n", r.stderr[-500:])
    sys.exit(1)
print("gmshToFoam OK")

# ── Patch: rename defaultFaces → axisCap, fix types ───────────────────────────
bnd = os.path.join(CASE_DIR, "constant", "polyMesh", "boundary")
txt = open(bnd).read()

patches = re.findall(r'^\s+(\w+)\s*$', txt, re.MULTILINE)
print("Patches from gmshToFoam:", patches)

txt = re.sub(r'\bdefaultFaces\b', 'axisCap', txt)
for name, ptype in [("frontAndBack","empty"), ("axis","symmetry"),
                     ("axisCap","symmetry"), ("body","wall")]:
    txt = re.sub(
        rf'({re.escape(name)}\s*\{{[^}}]*?type\s+)\w+\s*;',
        rf'\g<1>{ptype};', txt, flags=re.DOTALL)
    txt = re.sub(
        rf'({re.escape(name)}\s*\{{[^}}]*?physicalType\s+)\w+\s*;',
        rf'\g<1>{ptype};', txt, flags=re.DOTALL)

open(bnd, "w").write(txt)
print("Boundary patched.")

# ── checkMesh ─────────────────────────────────────────────────────────────────
print("\n── checkMesh " + "─"*50)
r = subprocess.run(["checkMesh", "-case", CASE_DIR], capture_output=True, text=True)
out = r.stdout + r.stderr
for line in out.splitlines():
    if any(k in line for k in ["***", "Failed", "OK", "negative", "open cell",
                                "unused", "oriented", "non-ortho", "skew", "Mesh non"]):
        print(line)
