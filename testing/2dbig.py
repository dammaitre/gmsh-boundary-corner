"""
2dbig.py — large-Reynolds BL stress test: 2×100m × 2×27m ellipse with a very
thin first layer (h1=1e-4, 37 layers, r=1.2).  Geometry is a spline-sampled
ellipse (same approach as 2dscatter.py).  foamRun is skipped; the test checks
that meshing and gmshToFoam succeed without geometry or topology failures.

Usage:
    python 2dbig.py           # headless
    python 2dbig.py --gui     # open the mesh in the gmsh GUI
"""

import sys, os, math, re, subprocess, shutil
import numpy as np

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry ──────────────────────────────────────────────────────────────────
A      = 100.0   # ellipse semi-major axis (axial, x)  [m]
B      =  27.0   # ellipse semi-minor axis (radial, y) [m]
X_FAR  = 300.0   # domain half-length along x
Y_FAR  = 150.0   # domain half-height along y

# ── Scatter sampling ──────────────────────────────────────────────────────────
# 1-cos distribution clusters points near both endpoints of each arc.
N_SCATTER = 200  # interior sample points per arc

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 10.0
LC_BODY =   0.50
#LC_NOSE =   5e-3
LC_NOSE = .01

# ── BoundaryCorner parameters ─────────────────────────────────────────────────
H1        = 1e-4    # first BL layer height
RATIO     = 1.20    # BL growth ratio
N_LAY     = 37      # number of BL layers
N_COLS    = 10      # fan columns near axis corner
#COL_WIDTH = 5e-3    # ColWidth = 50 * H1 (innermost column arc-length at AxisPoint)
COL_WIDTH = .01
W0_MAX    = LC_BODY # MaxColumnWidth (constant-zone column arc-length)

# ── Extrusion ─────────────────────────────────────────────────────────────────
DZ = 1.0

# ── Output ────────────────────────────────────────────────────────────────────
CASE_DIR = os.path.join(_root, "testing", "ellipse2dbig_of")
MAX_TIME = 5000

# ─────────────────────────────────────────────────────────────────────────────
gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("ellipse2dbig")

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
def make_arc_points(t_start, t_end, lc=LC_BODY):
    """Return list of interior gmsh point tags sampled with 1-cos distribution in t."""
    js = np.arange(1, N_SCATTER + 1)
    s  = (1.0 - np.cos(np.pi * js / (N_SCATTER + 1))) / 2.0
    ts = t_start + (t_end - t_start) * s
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
    gmsh.model.mesh.field.setNumber (f, "ColWidth",        COL_WIDTH)
    gmsh.model.mesh.field.setNumber (f, "Omega",           1.0)
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
print(f"Ellipse: {2*A}m × {2*B}m  |  BL: h1={H1:.0e}, {N_LAY} layers, r={RATIO}")
print(f"ColWidth={COL_WIDTH:.0e}, NbCornerColumns={N_COLS}, W0_MAX={W0_MAX}\n")
gmsh.model.mesh.generate(3)

# ── Inspect axis-column outer vertex positions ────────────────────────────────
print("\n── Axis-column outer vertex check (nose, y must be 0) ────────────────")
nodes, coords, _ = gmsh.model.mesh.getNodes()
coords = np.array(coords).reshape(-1, 3)

x_expected = []
for k in range(1, N_LAY + 1):
    hk = H1 * (RATIO**k - 1.0) / (RATIO - 1.0)
    x_expected.append(A + hk)

print(f"  Expected x range: [{x_expected[0]:.6f} .. {x_expected[-1]:.4f}], y=0 for all")
print()

found = []
for i in range(len(nodes)):
    x, y = coords[i, 0], coords[i, 1]
    for xe in x_expected:
        if abs(x - xe) < H1 * 0.1 and abs(y) < H1 * 0.5:
            found.append((x, y))
            break
found.sort()

ok = True
for x, y in found:
    flag = "" if abs(y) < 1e-10 else "  <- OFF AXIS"
    print(f"  x={x:.6f}  y={y:.10f}{flag}")
    if abs(y) > 1e-10:
        ok = False

if not found:
    print("  (no axis-column outer vertices matched — check tolerances)")
elif ok:
    print(f"\n  All {len(found)} axis-column outer vertices on y=0 OK")

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

print("\nDone — big-Re mesh in", CASE_DIR)
print("(foamRun skipped — examine mesh quality from checkMesh output above)")
