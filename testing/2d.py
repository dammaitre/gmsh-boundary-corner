"""
2d.py — OpenFOAM 2D axisymmetric wedge mesh for the full half-ellipse domain.

Geometry identical to double.py.  After 2D meshing the surface is extruded
by a 5° rotation around the x-axis (symmetry axis) to produce a thin wedge
volume suitable for OpenFOAM wedge boundary conditions.

Post-processing:
  gmshToFoam  → reads the .msh file into an OpenFOAM case
  fix_boundary → sets correct patch types (wedge/axis/wall) in polyMesh/boundary
  checkMesh   → validates the resulting polyMesh

Run:  python 2d.py           # headless
      python 2d.py --gui     # open 3D mesh in gmsh GUI before OpenFOAM steps
"""

import sys, os, math, re, subprocess, shutil

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry parameters (same as double.py) ───────────────────────────────────
A      = 2.0    # ellipse semi-major axis (axial, x)
B      = 1.2    # ellipse semi-minor axis (radial, y)
X_FAR  = 6.0
Y_FAR  = 3.5

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.30
LC_BODY = 0.06
LC_NOSE = 0.03

# ── BoundaryCorner parameters ─────────────────────────────────────────────────
H1     = 0.012
RATIO  = 1.20
N_LAY  = 6
N_COLS = 8
W0_MAX = LC_BODY

# ── Wedge angle (OpenFOAM 2D axisymmetric standard) ──────────────────────────
WEDGE_ANGLE = math.radians(5.0)

# ── Output ────────────────────────────────────────────────────────────────────
CASE_DIR = "/tmp/ellipse2d_of"

# ─────────────────────────────────────────────────────────────────────────────
gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("ellipse2d")

# ── Points ────────────────────────────────────────────────────────────────────
p_origin = gmsh.model.geo.addPoint(  0,      0,     0, LC_BODY)
p_nose   = gmsh.model.geo.addPoint(  A,      0,     0, LC_NOSE)   # right axis point
p_tail   = gmsh.model.geo.addPoint( -A,      0,     0, LC_NOSE)   # left  axis point
p_top    = gmsh.model.geo.addPoint(  0,      B,     0, LC_BODY)
p_ax_r   = gmsh.model.geo.addPoint(  X_FAR,  0,     0, LC_FAR)
p_ax_l   = gmsh.model.geo.addPoint( -X_FAR,  0,     0, LC_FAR)
p_far_tr = gmsh.model.geo.addPoint(  X_FAR,  Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint( -X_FAR,  Y_FAR, 0, LC_FAR)

# ── Curves ────────────────────────────────────────────────────────────────────
arc_front = gmsh.model.geo.addEllipseArc(p_top,  p_origin, p_nose, p_nose)  # (0,B)→(A,0)
arc_back  = gmsh.model.geo.addEllipseArc(p_tail, p_origin, p_nose, p_top)   # (-A,0)→(0,B)

l_ax_r  = gmsh.model.geo.addLine(p_ax_r,   p_nose)    # downstream axis → nose
l_right = gmsh.model.geo.addLine(p_ax_r,   p_far_tr)
l_top   = gmsh.model.geo.addLine(p_far_tr, p_far_tl)
l_left  = gmsh.model.geo.addLine(p_far_tl, p_ax_l)
l_ax_l  = gmsh.model.geo.addLine(p_ax_l,   p_tail)    # upstream axis → tail

cl = gmsh.model.geo.addCurveLoop(
    [-l_ax_r, l_right, l_top, l_left, l_ax_l, arc_back, arc_front])
sf = gmsh.model.geo.addPlaneSurface([cl])

# ── 3D wedge extrusion (geometry, before meshing) ────────────────────────────
# Rotate the 2D surface around the x-axis by WEDGE_ANGLE.
#   numElements=[1] : one circumferential cell layer
#   recombine=True  : quads → hexes, tris → prisms
# revolve output: [0]=(2,s_top), [1]=(3,vol), [2..N]= lateral surfaces.
# Axis curves (y=0) are degenerate under rotation and are SKIPPED by revolve
# — they produce no lateral surface.
out   = gmsh.model.geo.revolve(
    [(2, sf)], 0, 0, 0, 1, 0, 0, WEDGE_ANGLE,
    numElements=[1], recombine=True
)
s_top = out[0][1]   # rotated copy of sf (wedge_back face)
vol   = out[1][1]   # 3D volume

gmsh.model.geo.synchronize()

# ── Map each original boundary curve → its lateral surface ───────────────────
# The revolve output order does NOT match getBoundary order, so we identify
# each lateral surface by checking which original curve appears in its boundary.
orig_curves = {abs(c) for _, c in gmsh.model.getBoundary([(2, sf)], oriented=True)}

curve_to_lat = {}
for _, s in out[2:]:
    bnd = gmsh.model.getBoundary([(2, s)], oriented=False)
    for _, c in bnd:
        ctag = abs(c)
        if ctag in orig_curves:
            curve_to_lat[ctag] = s
            break

# ── BoundaryCorner field — front (nose) ───────────────────────────────────────
bc_front = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc_front, "CurvesList",      [arc_front])
gmsh.model.mesh.field.setNumbers(bc_front, "AxisPoint",       [A, 0.0])
gmsh.model.mesh.field.setNumber (bc_front, "Size",            H1)
gmsh.model.mesh.field.setNumber (bc_front, "Ratio",           RATIO)
gmsh.model.mesh.field.setNumber (bc_front, "NbLayers",        N_LAY)
gmsh.model.mesh.field.setNumber (bc_front, "NbCornerColumns", N_COLS)
gmsh.model.mesh.field.setNumber (bc_front, "MaxColumnWidth",  W0_MAX)
gmsh.model.mesh.field.setNumber (bc_front, "SkipAxisColumn",  1)   # avoid y=0 cells under revolve

# ── BoundaryCorner field — back (tail) ────────────────────────────────────────
bc_back = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc_back, "CurvesList",      [arc_back])
gmsh.model.mesh.field.setNumbers(bc_back, "AxisPoint",       [-A, 0.0])
gmsh.model.mesh.field.setNumber (bc_back, "Size",            H1)
gmsh.model.mesh.field.setNumber (bc_back, "Ratio",           RATIO)
gmsh.model.mesh.field.setNumber (bc_back, "NbLayers",        N_LAY)
gmsh.model.mesh.field.setNumber (bc_back, "NbCornerColumns", N_COLS)
gmsh.model.mesh.field.setNumber (bc_back, "MaxColumnWidth",  W0_MAX)
gmsh.model.mesh.field.setNumber (bc_back, "SkipAxisColumn",  1)   # avoid y=0 cells under revolve

min_f = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_f, "FieldsList", [bc_front, bc_back])
gmsh.model.mesh.field.setAsBackgroundMesh(min_f)

gmsh.option.setNumber("Mesh.BoundaryCornerField",     bc_front)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

# ── 3D mesh (sweeps 2D mesh of sf into the wedge volume) ─────────────────────
gmsh.model.mesh.generate(3)

# ── Physical groups (OpenFOAM patch names) ────────────────────────────────────
def pg(dim, tags, name):
    t = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, t, name)

pg(3, [vol],                                              "fluid")
pg(2, [sf],                                               "wedge_front")  # θ = 0
pg(2, [s_top],                                            "wedge_back")   # θ = WEDGE_ANGLE
pg(2, [curve_to_lat[arc_front], curve_to_lat[arc_back]], "body")
pg(2, [curve_to_lat[l_right]],                            "outlet")
pg(2, [curve_to_lat[l_top]],                              "top")
pg(2, [curve_to_lat[l_left]],                             "inlet")

# ── Write mesh ────────────────────────────────────────────────────────────────
os.makedirs(CASE_DIR, exist_ok=True)
msh_path = os.path.join(CASE_DIR, "ellipse.msh")
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # gmshToFoam reads v2
gmsh.write(msh_path)
print(f"\nMesh written: {msh_path}")

# ── Optional GUI preview ──────────────────────────────────────────────────────
if gui:
    dev_bin  = os.path.join(_root, "build", "gmsh")
    gmsh_bin = dev_bin if shutil.which(dev_bin) else shutil.which("gmsh") or dev_bin
    cache    = os.path.join(_root, "build", "CMakeCache.txt")
    has_fltk = False
    if os.path.exists(cache):
        with open(cache) as _f:
            for _line in _f:
                if "FLTK_BASE_LIBRARY_RELEASE:FILEPATH=" in _line and "NOTFOUND" not in _line:
                    has_fltk = True
                    break
    if not has_fltk:
        gmsh_bin = shutil.which("gmsh") or dev_bin
    subprocess.run([gmsh_bin, msh_path])

gmsh.finalize()

# ── Minimal OpenFOAM case skeleton ────────────────────────────────────────────
sys_dir = os.path.join(CASE_DIR, "system")
os.makedirs(sys_dir, exist_ok=True)

ctrl_path = os.path.join(sys_dir, "controlDict")
if not os.path.exists(ctrl_path):
    with open(ctrl_path, "w") as f:
        f.write(
            "FoamFile\n{\n"
            "    version     2.0;\n"
            "    format      ascii;\n"
            "    class       dictionary;\n"
            "    location    \"system\";\n"
            "    object      controlDict;\n"
            "}\n\n"
            "application     icoFoam;\n"
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
r = subprocess.run(
    ["gmshToFoam", os.path.basename(msh_path)],
    cwd=CASE_DIR
)
if r.returncode != 0:
    print("\ngmshToFoam failed — is OpenFOAM sourced in this shell?")
    sys.exit(r.returncode)

# ── Fix boundary patch types ──────────────────────────────────────────────────
# gmshToFoam writes every patch as type=patch.
# OpenFOAM axisymmetric requires:
#   wedge       → type wedge
#   defaultFaces → rename to axis, type axis   (boundary faces on y=0 that
#                  revolve skips because axis curves produce no gmsh surface)
#   body        → type wall
#
# Note: the 14 gmshToFoam "Could not match" warnings are the degenerate
# axis-column quads from BoundaryCornerField (the outermost BC column sits
# entirely on y=0 and produces zero-volume cells under revolve).  Those faces
# are discarded by gmshToFoam; the resulting 12 "unused points" in checkMesh
# are their vertices.  This is a known limitation of revolving a mesh that
# includes BoundaryCorner axis columns.

# Desired types per patch name (after renaming defaultFaces→axis)
_PATCH_TYPES = {
    "wedge_front": "wedge",
    "wedge_back":  "wedge",
    "axis":        "axis",
    "body":   "wall",
    "inlet":  "patch",
    "outlet": "patch",
    "top":    "patch",
}

def fix_of_boundary(boundary_path):
    with open(boundary_path) as f:
        txt = f.read()

    # Rename defaultFaces → axis
    txt = re.sub(r'\bdefaultFaces\b', 'axis', txt)

    # For each named patch block, set type and physicalType
    def _fix_block(m):
        name  = m.group(1)
        inner = m.group(2)
        ptype = _PATCH_TYPES.get(name, "patch")
        inner = re.sub(r'(type\s+)\w+(;)',         rf'\g<1>{ptype}\2',  inner)
        inner = re.sub(r'(physicalType\s+)\w+(;)', rf'\g<1>{ptype}\2',  inner)
        return f'{name}\n    {{\n{inner}    }}'

    txt = re.sub(
        r'(\w+)\s*\n\s*\{\s*\n(.*?)\n\s*\}',
        _fix_block,
        txt,
        flags=re.DOTALL
    )

    with open(boundary_path, "w") as f:
        f.write(txt)
    print(f"Patched: {boundary_path}")

boundary_path = os.path.join(CASE_DIR, "constant", "polyMesh", "boundary")
fix_of_boundary(boundary_path)

# ── transformPoints ───────────────────────────────────────────────────────────
# The mesh currently spans θ ∈ [0°, 5°].  OpenFOAM wedge requires the two
# wedge faces to be symmetric about a coordinate plane, i.e. ±2.5°.
# Rotating by -2.5° around the x-axis achieves this without touching the
# 2D geometry or the meshing step.
print("\n── transformPoints " + "─" * 55)
# OF13 Rx(angle) is broken; use rotate=(n1 n2) to rotate y-axis by -WEDGE_ANGLE/2
_c = math.cos(-WEDGE_ANGLE / 2)
_s = math.sin(-WEDGE_ANGLE / 2)
_rot = f"rotate=((0 1 0)(0 {_c:.10f} {_s:.10f}))"
subprocess.run(["transformPoints", "-case", CASE_DIR, _rot])

# ── checkMesh ─────────────────────────────────────────────────────────────────
print("\n── checkMesh " + "─" * 61)
subprocess.run(["checkMesh", "-case", CASE_DIR])
