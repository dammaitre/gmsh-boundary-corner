"""
2dscatter.py — same geometry as 2dplanar.py but the ellipse arcs are replaced
by discrete/scattered numpy point clouds fed into a gmsh Spline.

Expected issues vs. the analytic arc version:
  • The spline tangent at nose (A, 0) and tail (-A, 0) is not exactly horizontal
    because it depends on the slope between the last two sample points, not the
    true ellipse tangent.  So the outward normal there is not exactly ±x and the
    BoundaryCorner last-column vertices are not exactly on y = 0.
  • Increasing N_SCATTER reduces (but never eliminates) the discretisation error.

After meshing, gmshToFoam + checkMesh are run.  foamRun is skipped so the
experiment can be examined without a full OpenFOAM solve.

Usage:
    python 2dscatter.py           # headless
    python 2dscatter.py --gui     # open the mesh in the gmsh GUI
"""

import sys, os, math, re, subprocess, shutil
import numpy as np

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry ──────────────────────────────────────────────────────────────────
A      = 2.0    # ellipse semi-major axis (axial, x)
B      = 1.2    # ellipse semi-minor axis (radial, y)
X_FAR  = 6.0
Y_FAR  = 3.5

# ── Scatter sampling ──────────────────────────────────────────────────────────
# Number of *interior* sample points per arc (endpoints are the GVertex points).
# Low value → visible tangent error at axis; high value → error is tiny but still
# present (unlike the analytic arc which is exact).
N_SCATTER = 12

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.30
LC_BODY = 0.06
LC_NOSE = 0.03

# ── BoundaryCorner parameters ─────────────────────────────────────────────────
H1     = 0.012
RATIO  = 1.15
N_LAY  = 6
N_COLS = 8
W0_MAX = LC_BODY

# ── Extrusion ─────────────────────────────────────────────────────────────────
DZ = 1.0

# ── Output ────────────────────────────────────────────────────────────────────
CASE_DIR = os.path.join(_root, "testing", "ellipse2dscatter_of")
MAX_TIME = 5000

# ─────────────────────────────────────────────────────────────────────────────
gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("ellipse2dscatter")

# ── Fixed anchor points ───────────────────────────────────────────────────────
p_origin = gmsh.model.geo.addPoint(  0,      0,     0, LC_BODY)
p_nose   = gmsh.model.geo.addPoint(  A,      0,     0, LC_NOSE)
p_tail   = gmsh.model.geo.addPoint( -A,      0,     0, LC_NOSE)
p_top    = gmsh.model.geo.addPoint(  0,      B,     0, LC_BODY)
p_ax_r   = gmsh.model.geo.addPoint(  X_FAR,  0,     0, LC_FAR)
p_ax_l   = gmsh.model.geo.addPoint( -X_FAR,  0,     0, LC_FAR)
p_far_tr = gmsh.model.geo.addPoint(  X_FAR,  Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint( -X_FAR,  Y_FAR, 0, LC_FAR)

# ── Build scattered-point arcs ────────────────────────────────────────────────
# Ellipse parametric: x(t) = A*cos(t), y(t) = B*sin(t)
#   Front arc:  t in [pi/2 .. 0]  → (0,B) → (A,0)
#   Back arc:   t in [pi  .. pi/2] → (-A,0) → (0,B)
#
# Interior points only — endpoints reuse existing GVertex tags so the topology
# is correct.  The spline will NOT have an exact horizontal tangent at t=0,π
# because the slope is estimated from the last two sample points.

def make_arc_points(t_start, t_end, lc=LC_BODY):
    """Return list of interior gmsh point tags sampled uniformly in t."""
    ts = np.linspace(t_start, t_end, N_SCATTER + 2)[1:-1]  # drop endpoints
    tags = []
    for t in ts:
        x = A * math.cos(t)
        y = B * math.sin(t)
        tags.append(gmsh.model.geo.addPoint(x, y, 0, lc))
    return tags

# Front arc: (0,B)=p_top → interior → (A,0)=p_nose
front_interior = make_arc_points(math.pi / 2, 0.0)
arc_front = gmsh.model.geo.addSpline([p_top] + front_interior + [p_nose])

# Back arc: (-A,0)=p_tail → interior → (0,B)=p_top
back_interior = make_arc_points(math.pi, math.pi / 2)
arc_back = gmsh.model.geo.addSpline([p_tail] + back_interior + [p_top])

# ── Far-field straight lines ──────────────────────────────────────────────────
l_ax_r  = gmsh.model.geo.addLine(p_ax_r,   p_nose)
l_right = gmsh.model.geo.addLine(p_ax_r,   p_far_tr)
l_top   = gmsh.model.geo.addLine(p_far_tr, p_far_tl)
l_left  = gmsh.model.geo.addLine(p_far_tl, p_ax_l)
l_ax_l  = gmsh.model.geo.addLine(p_ax_l,   p_tail)

cl = gmsh.model.geo.addCurveLoop(
    [-l_ax_r, l_right, l_top, l_left, l_ax_l, arc_back, arc_front])
sf = gmsh.model.geo.addPlaneSurface([cl])

# ── 3D extrusion (geometry only, before meshing) ──────────────────────────────
out   = gmsh.model.geo.extrude(
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

# ── BoundaryCorner fields ──────────────────────────────────────────────────────
def add_bc_field(curves, axis_point):
    f = gmsh.model.mesh.field.add("BoundaryCorner")
    gmsh.model.mesh.field.setNumbers(f, "CurvesList",      curves)
    gmsh.model.mesh.field.setNumbers(f, "AxisPoint",       axis_point)
    gmsh.model.mesh.field.setNumber (f, "Size",            H1)
    gmsh.model.mesh.field.setNumber (f, "Ratio",           RATIO)
    gmsh.model.mesh.field.setNumber (f, "NbLayers",        N_LAY)
    gmsh.model.mesh.field.setNumber (f, "NbCornerColumns", N_COLS)
    gmsh.model.mesh.field.setNumber (f, "MaxColumnWidth",  W0_MAX)
    return f

bc_front = add_bc_field([arc_front], [A,  0.0])
bc_back  = add_bc_field([arc_back],  [-A, 0.0])

min_f = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_f, "FieldsList", [bc_front, bc_back])
gmsh.model.mesh.field.setAsBackgroundMesh(min_f)

gmsh.option.setNumber("Mesh.BoundaryCornerField",     bc_front)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

# ── Generate mesh ─────────────────────────────────────────────────────────────
print(f"\nBuilding spline arcs with {N_SCATTER} interior scatter points per arc.")
print("Note: spline tangent at axis endpoints is approximate → normal not exactly ±x\n")
gmsh.model.mesh.generate(3)

# ── Inspect axis-column outer vertex positions ────────────────────────────────
# The axis-column outer vertices (nose side) must sit exactly on y=0.
# Expected x positions: A + h1*omega*(ratio^k - 1)/(ratio - 1) for k=1..N_LAY.
print("\n── Axis-column outer vertex check (nose, y must be 0) ────────────────")
nodes, coords, _ = gmsh.model.mesh.getNodes()
coords = np.array(coords).reshape(-1, 3)

# Compute expected x positions for the axis-column outer vertices
x_expected = []
for k in range(1, N_LAY + 1):
    hk = H1 * (RATIO**k - 1.0) / (RATIO - 1.0)
    x_expected.append(A + hk)

print(f"  Expected x range: [{x_expected[0]:.4f} .. {x_expected[-1]:.4f}], y=0 for all")
print()

found = []
for i in range(len(nodes)):
    x, y = coords[i, 0], coords[i, 1]
    # Match against expected x positions within a tight tolerance
    for xe in x_expected:
        if abs(x - xe) < H1 * 0.1 and abs(y) < H1 * 0.5:
            found.append((x, y))
            break
found.sort()

ok = True
for x, y in found:
    flag = "" if abs(y) < 1e-10 else "  ← OFF AXIS"
    print(f"  x={x:.6f}  y={y:.10f}{flag}")
    if abs(y) > 1e-10:
        ok = False

if not found:
    print("  (no axis-column outer vertices matched — check tolerances)")
elif ok:
    print(f"\n  All {len(found)} axis-column outer vertices on y=0 ✓")

# ── Physical groups ────────────────────────────────────────────────────────────
def pg(dim, tags, name):
    t = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, t, name)

pg(3, [vol],                                              "fluid")
pg(2, [sf],                                               "front")
pg(2, [s_top],                                            "back")
pg(2, [curve_to_lat[arc_front], curve_to_lat[arc_back]], "body")
pg(2, [curve_to_lat[l_right]],                            "outlet")
pg(2, [curve_to_lat[l_top]],                              "top")
pg(2, [curve_to_lat[l_left]],                             "inlet")
pg(2, [curve_to_lat[l_ax_r], curve_to_lat[l_ax_l]],      "symmetry")

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

# ── OpenFOAM case skeleton ────────────────────────────────────────────────────
sys_dir   = os.path.join(CASE_DIR, "system")
zero_dir  = os.path.join(CASE_DIR, "0")
const_dir = os.path.join(CASE_DIR, "constant")
for d in (sys_dir, zero_dir, const_dir):
    os.makedirs(d)

def foam_header(cls, location, obj):
    return (
        "FoamFile\n{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        f"    class       {cls};\n"
        f"    location    \"{location}\";\n"
        f"    object      {obj};\n"
        "}\n\n"
    )

with open(os.path.join(sys_dir, "controlDict"), "w") as f:
    f.write(
        foam_header("dictionary", "system", "controlDict") +
        "application     foamRun;\n"
        "solver          incompressibleFluid;\n\n"
        "startFrom       startTime;\n"
        "startTime       0;\n"
        "stopAt          endTime;\n"
        f"endTime         {MAX_TIME};\n"
        "deltaT          1;\n\n"
        "writeControl    timeStep;\n"
        "writeInterval   100;\n\n"
        "runTimeModifiable yes;\n"
    )

with open(os.path.join(sys_dir, "fvSchemes"), "w") as f:
    f.write(
        foam_header("dictionary", "system", "fvSchemes") +
        "ddtSchemes  { default steadyState; }\n\n"
        "gradSchemes { default Gauss linear; }\n\n"
        "divSchemes\n{\n"
        "    default                              none;\n"
        "    div(phi,U)                           Gauss linearUpwind grad(U);\n"
        "    div(phi,k)                           Gauss limitedLinear 1;\n"
        "    div(phi,omega)                       Gauss limitedLinear 1;\n"
        "    div((nuEff*dev(T(grad(U)))))         Gauss linear;\n"
        "}\n\n"
        "laplacianSchemes { default Gauss linear corrected; }\n\n"
        "interpolationSchemes { default linear; }\n\n"
        "snGradSchemes { default corrected; }\n\n"
        "wallDist { method meshWave; }\n"
    )

with open(os.path.join(sys_dir, "fvSolution"), "w") as f:
    f.write(
        foam_header("dictionary", "system", "fvSolution") +
        "solvers\n{\n"
        "    p\n    {\n"
        "        solver          GAMG;\n"
        "        smoother        GaussSeidel;\n"
        "        tolerance       1e-6;\n"
        "        relTol          0.01;\n"
        "    }\n"
        "    U\n    {\n"
        "        solver          smoothSolver;\n"
        "        smoother        GaussSeidel;\n"
        "        tolerance       1e-8;\n"
        "        relTol          0.1;\n"
        "    }\n"
        "    \"(k|omega)\"\n    {\n"
        "        solver          smoothSolver;\n"
        "        smoother        GaussSeidel;\n"
        "        tolerance       1e-8;\n"
        "        relTol          0.1;\n"
        "    }\n"
        "}\n\n"
        "SIMPLE\n{\n"
        "    nNonOrthogonalCorrectors 1;\n"
        "    residualControl\n    {\n"
        "        p       1e-4;\n"
        "        U       1e-4;\n"
        "        k       1e-4;\n"
        "        omega   1e-4;\n"
        "    }\n"
        "}\n\n"
        "relaxationFactors\n{\n"
        "    fields      { p 0.3; }\n"
        "    equations   { U 0.7; k 0.7; omega 0.7; }\n"
        "}\n"
    )

with open(os.path.join(const_dir, "momentumTransport"), "w") as f:
    f.write(
        foam_header("dictionary", "constant", "momentumTransport") +
        "simulationType  RAS;\n\n"
        "RAS\n{\n"
        "    model           kOmegaSST;\n"
        "    turbulence      on;\n"
        "    printCoeffs     on;\n"
        "}\n"
    )

with open(os.path.join(const_dir, "transportProperties"), "w") as f:
    f.write(
        foam_header("dictionary", "constant", "transportProperties") +
        "viscosityModel  constant;\n"
        "nu              1.5e-5;\n"
    )

_NU        = 1.5e-5
_U_INF     = 44.44
_k_inf     = 1.5 * (_U_INF * 0.01) ** 2
_omega_inf = _k_inf / (10.0 * _NU)

_PATCHES = ("inlet", "outlet", "top", "body", "symmetry", "symAxis", "front", "back")

def _bc(patch, field):
    empty  = "{ type empty; }"
    sym    = "{ type symmetry; }"
    zeroG  = "{ type zeroGradient; }"
    if patch in ("front", "back"):
        return empty
    if patch in ("symmetry", "symAxis"):
        return sym
    if field == "U":
        if patch == "inlet":   return "{ type fixedValue; value uniform (44.44 0 0); }"
        if patch == "body":    return "{ type noSlip; }"
        return zeroG
    if field == "p":
        if patch == "outlet":  return "{ type fixedValue; value uniform 0; }"
        return zeroG
    if field == "k":
        if patch == "inlet":   return f"{{ type fixedValue; value uniform {_k_inf:.4g}; }}"
        if patch == "body":    return "{ type kqRWallFunction; value uniform 0; }"
        return zeroG
    if field == "omega":
        if patch == "inlet":   return f"{{ type fixedValue; value uniform {_omega_inf:.4g}; }}"
        if patch == "body":    return "{ type omegaWallFunction; value uniform 1; }"
        return zeroG
    if field == "nut":
        if patch == "body":    return "{ type nutkWallFunction; value uniform 0; }"
        return "{ type calculated; value uniform 0; }"
    return zeroG

def write_field(path, cls, dims, internal, field_name):
    lines = [foam_header(cls, "0", field_name),
             f"dimensions  {dims};\n",
             f"internalField {internal};\n\n",
             "boundaryField\n{\n"]
    for p in _PATCHES:
        lines.append(f"    {p:<12} {_bc(p, field_name)}\n")
    lines.append("}\n")
    with open(path, "w") as f:
        f.writelines(lines)

write_field(os.path.join(zero_dir, "U"),
            "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (44.44 0 0)", "U")
write_field(os.path.join(zero_dir, "p"),
            "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", "p")
write_field(os.path.join(zero_dir, "k"),
            "volScalarField", "[0 2 -2 0 0 0 0]", f"uniform {_k_inf:.4g}", "k")
write_field(os.path.join(zero_dir, "omega"),
            "volScalarField", "[0 0 -1 0 0 0 0]", f"uniform {_omega_inf:.4g}", "omega")
write_field(os.path.join(zero_dir, "nut"),
            "volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", "nut")

# ── gmshToFoam ────────────────────────────────────────────────────────────────
print("\n── gmshToFoam " + "─" * 60)
r = subprocess.run(["gmshToFoam", os.path.basename(msh_path)], cwd=CASE_DIR)
if r.returncode != 0:
    print("\ngmshToFoam failed — is OpenFOAM sourced in this shell?")
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
subprocess.run(["checkMesh", "-case", CASE_DIR])

print("\nDone — scatter mesh in", CASE_DIR)
print("(foamRun skipped — examine mesh quality from checkMesh output above)")
