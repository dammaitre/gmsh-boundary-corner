"""
doublecorner.py — BoundaryDoubleCorner on a full half-ellipse as a SINGLE spline.

Geometry (external axisymmetric CFD, upper-half 2D cross-section):
  - Ellipse semi-axes a=2.0 (axial), b=1.2 (radial)
  - Body surface = single spline from nose (a,0) through top (0,b) to tail (-a,0)
  - Fluid domain = region outside the ellipse, bounded by a rectangular far-field

One BoundaryDoubleCornerField covers the entire spline:
  - NbCornerColumns compressed columns near the nose (right axis intersection)
  - Constant-width columns in the mid-section
  - NbCornerColumns compressed columns near the tail (left axis intersection)

Run:  python doublecorner.py           # headless — mesh + gmshToFoam + checkMesh
      python doublecorner.py --gui     # open result in gmsh GUI after meshing
"""

import sys, os, math, re, subprocess, shutil
import numpy as np

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry parameters ───────────────────────────────────────────────────────
A = 2.0     # ellipse semi-major axis (axial, x direction)
B = 1.2     # ellipse semi-minor axis (radial, y direction)
X_FAR = 6.0
Y_FAR = 3.5

N_SCATTER = 200   # interior scatter points on the spline

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.30
LC_BODY = 0.06
LC_NOSE = 0.03

# ── BoundaryDoubleCorner parameters ──────────────────────────────────────────
H1     = 0.012
RATIO  = 1.20
N_LAY  = 6
N_COLS = 8
W0_MAX = LC_BODY

# ── Extrusion ─────────────────────────────────────────────────────────────────
DZ = 1.0

# ── Output ────────────────────────────────────────────────────────────────────
CASE_DIR   = os.path.join(_root, "testing", "doublecorner_of")
OF_BASHRC  = "/opt/openfoam13/etc/bashrc"

# ─────────────────────────────────────────────────────────────────────────────
gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("doublecorner_ellipse")

# ── Points ────────────────────────────────────────────────────────────────────
p_nose   = gmsh.model.geo.addPoint(  A,      0,     0, LC_NOSE)
p_tail   = gmsh.model.geo.addPoint( -A,      0,     0, LC_NOSE)
p_ax_r   = gmsh.model.geo.addPoint(  X_FAR,  0,     0, LC_FAR)
p_ax_l   = gmsh.model.geo.addPoint( -X_FAR,  0,     0, LC_FAR)
p_far_tr = gmsh.model.geo.addPoint(  X_FAR,  Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint( -X_FAR,  Y_FAR, 0, LC_FAR)

# ── Single spline: nose → top → tail  (t from 0 to pi, 1-cos distribution) ──
js = np.arange(1, N_SCATTER + 1)
s  = (1.0 - np.cos(np.pi * js / (N_SCATTER + 1))) / 2.0
ts = np.pi * s
interior = []
for t in ts:
    interior.append(gmsh.model.geo.addPoint(A * math.cos(t), B * math.sin(t), 0, LC_BODY))

arc_profile = gmsh.model.geo.addSpline([p_nose] + interior + [p_tail])

# ── Axis and far-field lines ──────────────────────────────────────────────────
l_ax_r  = gmsh.model.geo.addLine(p_ax_r,   p_nose)
l_right = gmsh.model.geo.addLine(p_ax_r,   p_far_tr)
l_top   = gmsh.model.geo.addLine(p_far_tr, p_far_tl)
l_left  = gmsh.model.geo.addLine(p_far_tl, p_ax_l)
l_ax_l  = gmsh.model.geo.addLine(p_ax_l,   p_tail)

cl = gmsh.model.geo.addCurveLoop(
    [-l_ax_r, l_right, l_top, l_left, l_ax_l, -arc_profile])
sf = gmsh.model.geo.addPlaneSurface([cl])

# ── 3D extrusion (geometry only, before meshing) ──────────────────────────────
out = gmsh.model.geo.extrude(
    [(2, sf)], 0, 0, DZ,
    numElements=[1], recombine=True
)
s_top = out[0][1]
vol   = out[1][1]

gmsh.model.geo.synchronize()

# ── Map boundary curves → lateral surfaces ────────────────────────────────────
orig_curves = {abs(c) for _, c in gmsh.model.getBoundary([(2, sf)], oriented=True)}

curve_to_lat = {}
for _, s in out[2:]:
    for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
        ctag = abs(c)
        if ctag in orig_curves:
            curve_to_lat[ctag] = s
            break

# ── BoundaryDoubleCorner field ────────────────────────────────────────────────
bdc = gmsh.model.mesh.field.add("BoundaryDoubleCorner")
gmsh.model.mesh.field.setNumbers(bdc, "CurvesList",      [arc_profile])
gmsh.model.mesh.field.setNumbers(bdc, "NosePoint",       [A,  0.0])
gmsh.model.mesh.field.setNumbers(bdc, "TailPoint",       [-A, 0.0])
gmsh.model.mesh.field.setNumber (bdc, "Size",            H1)
gmsh.model.mesh.field.setNumber (bdc, "Ratio",           RATIO)
gmsh.model.mesh.field.setNumber (bdc, "NbLayers",        N_LAY)
gmsh.model.mesh.field.setNumber (bdc, "NbCornerColumns", N_COLS)
gmsh.model.mesh.field.setNumber (bdc, "MaxColumnWidth",  W0_MAX)

gmsh.model.mesh.field.setAsBackgroundMesh(bdc)
gmsh.option.setNumber("Mesh.BoundaryCornerField",     bdc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

# ── Generate ──────────────────────────────────────────────────────────────────
print(f"\nBuilding spline with {N_SCATTER} interior scatter points.")
gmsh.model.mesh.generate(3)

# ── Results ───────────────────────────────────────────────────────────────────
elem_types, elem_tags, _ = gmsh.model.mesh.getElements(2)
print("\n=== 2D element counts ===")
type_names = {2: "triangle", 3: "quadrangle"}
n_quads = n_tris = 0
for t, tags in zip(elem_types, elem_tags):
    if not len(tags):
        continue
    name = type_names.get(t, f"type{t}")
    print(f"  {name}: {len(tags)}")
    if t == 3: n_quads = len(tags)
    if t == 2: n_tris  = len(tags)

print(f"\n  NbCornerColumns={N_COLS} (each end),  NbLayers={N_LAY},  W0_MAX={W0_MAX}")
print(f"  total quads : {n_quads}")

# ── Physical groups ────────────────────────────────────────────────────────────
def pg(dim, tags, name):
    t = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, t, name)

pg(3, [vol],                                                    "fluid")
pg(2, [sf],                                                     "front")
pg(2, [s_top],                                                  "back")
pg(2, [curve_to_lat[arc_profile]],                              "body")
pg(2, [curve_to_lat[l_right]],                                  "outlet")
pg(2, [curve_to_lat[l_top]],                                    "top")
pg(2, [curve_to_lat[l_left]],                                   "inlet")
pg(2, [curve_to_lat[l_ax_r], curve_to_lat[l_ax_l]],            "symmetry")

# ── Write mesh ────────────────────────────────────────────────────────────────
if os.path.isdir(CASE_DIR):
    shutil.rmtree(CASE_DIR)
os.makedirs(CASE_DIR)
msh_path = os.path.join(CASE_DIR, "ellipse.msh")
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(msh_path)
print(f"\nMesh written: {msh_path}")

# ── Optional GUI preview ──────────────────────────────────────────────────────
if gui:
    dev_bin  = os.path.join(_root, "build", "gmsh")
    gmsh_bin = dev_bin if shutil.which(dev_bin) else shutil.which("gmsh") or dev_bin
    cache    = os.path.join(_root, "build", "CMakeCache.txt")
    has_fltk = any(
        "FLTK_BASE_LIBRARY_RELEASE:FILEPATH=" in line and "NOTFOUND" not in line
        for line in open(cache)
    ) if os.path.exists(cache) else False
    if not has_fltk:
        gmsh_bin = shutil.which("gmsh") or dev_bin
    subprocess.run([gmsh_bin, msh_path])

gmsh.finalize()

# ── Helper: run an OpenFOAM command with the OF environment sourced ───────────
def of_run(cmd, **kwargs):
    return subprocess.run(
        ["bash", "-c", f"source {OF_BASHRC} && {cmd}"],
        **kwargs
    )

# ── Minimal case skeleton (gmshToFoam requires system/controlDict) ─────────────
sys_dir = os.path.join(CASE_DIR, "system")
os.makedirs(sys_dir, exist_ok=True)
with open(os.path.join(sys_dir, "controlDict"), "w") as f:
    f.write(
        "FoamFile\n{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       dictionary;\n"
        "    location    \"system\";\n"
        "    object      controlDict;\n"
        "}\n\n"
        "application     foamRun;\n"
        "startFrom       startTime;\n"
        "startTime       0;\n"
        "stopAt          endTime;\n"
        "endTime         1;\n"
        "deltaT          1;\n"
        "writeControl    timeStep;\n"
        "writeInterval   1;\n"
    )

# ── gmshToFoam ────────────────────────────────────────────────────────────────
print("\n── gmshToFoam " + "─" * 60)
r = of_run(f"gmshToFoam {os.path.basename(msh_path)}", cwd=CASE_DIR)
if r.returncode != 0:
    print("\ngmshToFoam failed — check that OpenFOAM 13 is installed at", OF_BASHRC)
    sys.exit(r.returncode)

# ── Fix boundary patch types ──────────────────────────────────────────────────
_PATCH_TYPES = {
    "front":    "empty",
    "back":     "empty",
    "body":     "wall",
    "symmetry": "symmetry",
    "symAxis":  "symmetry",
    "inlet":    "patch",
    "outlet":   "patch",
    "top":      "patch",
}

def fix_of_boundary(boundary_path):
    with open(boundary_path) as f:
        txt = f.read()
    txt = re.sub(r'\bdefaultFaces\b', 'symAxis', txt)

    def _fix_block(m):
        name  = m.group(1)
        inner = m.group(2)
        ptype = _PATCH_TYPES.get(name, "patch")
        inner = re.sub(r'(type\s+)\w+(;)',         rf'\g<1>{ptype}\2',  inner)
        inner = re.sub(r'(physicalType\s+)\w+(;)', rf'\g<1>{ptype}\2',  inner)
        return f'{name}\n    {{\n{inner}    }}'

    txt = re.sub(
        r'(\w+)\s*\n\s*\{\s*\n(.*?)\n\s*\}',
        _fix_block, txt, flags=re.DOTALL
    )
    with open(boundary_path, "w") as f:
        f.write(txt)
    print(f"Patched: {boundary_path}")

fix_of_boundary(os.path.join(CASE_DIR, "constant", "polyMesh", "boundary"))

# ── checkMesh ─────────────────────────────────────────────────────────────────
print("\n── checkMesh " + "─" * 61)
of_run(f"checkMesh -case {CASE_DIR}")

print("\nDone — mesh in", CASE_DIR)
