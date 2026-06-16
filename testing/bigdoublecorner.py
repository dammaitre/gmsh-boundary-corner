"""
bigdoublecorner.py — large-Reynolds BoundaryDoubleCorner stress test.

Same geometry and BL parameters as 2dbig.py (A=100m, B=27m, h1=1e-4, 37 layers,
r=1.2) but the body is a SINGLE spline (nose → top → tail) covered by one
BoundaryDoubleCorner field instead of two separate BoundaryCorner fields.

Usage:
    python bigdoublecorner.py           # headless, mesh + checkMesh only
    python bigdoublecorner.py --gui     # open the mesh in the gmsh GUI
    python bigdoublecorner.py --of      # also run foamRun after checkMesh
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
N_SCATTER = 200  # interior sample points on the spline

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 10.0
LC_BODY =  0.25
LC_NOSE =  0.01

# ── BoundaryDoubleCorner parameters ──────────────────────────────────────────
H1        = 1e-4    # first BL layer height
RATIO     = 1.20    # BL growth ratio
N_LAY     = 37      # number of BL layers
N_COLS    = 10      # compressed columns near each axis corner
COL_WIDTH = 0.01    # arc-length of innermost column at each corner
W0_MAX    = LC_BODY # MaxColumnWidth (constant-zone column arc-length)

# ── Extrusion ─────────────────────────────────────────────────────────────────
DZ = 1.0

# ── Output ────────────────────────────────────────────────────────────────────
CASE_DIR = os.path.join(_root, "testing", "bigdoublecorner_of")
MAX_TIME = 5000

# ─────────────────────────────────────────────────────────────────────────────
gui    = "--gui" in sys.argv
run_of = "--of"  in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("bigdoublecorner")

# ── Points ────────────────────────────────────────────────────────────────────
p_nose   = gmsh.model.geo.addPoint(  A,      0,     0, LC_BODY)
p_tail   = gmsh.model.geo.addPoint( -A,      0,     0, LC_BODY)
p_ax_r   = gmsh.model.geo.addPoint(  X_FAR,  0,     0, LC_FAR)
p_ax_l   = gmsh.model.geo.addPoint( -X_FAR,  0,     0, LC_FAR)
p_far_tr = gmsh.model.geo.addPoint(  X_FAR,  Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint( -X_FAR,  Y_FAR, 0, LC_FAR)

# ── Single spline: nose → top → tail  (1-cos distribution in t) ──────────────
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
gmsh.model.mesh.field.setNumber (bdc, "ColWidth",        COL_WIDTH)

gmsh.model.mesh.field.setAsBackgroundMesh(bdc)
gmsh.option.setNumber("Mesh.BoundaryCornerField",     bdc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

# ── Generate mesh ─────────────────────────────────────────────────────────────
print(f"\nBuilding spline with {N_SCATTER} interior scatter points.")
print(f"Ellipse: {2*A}m × {2*B}m  |  BL: h1={H1:.0e}, {N_LAY} layers, r={RATIO}")
print(f"ColWidth={COL_WIDTH:.0e}, NbCornerColumns={N_COLS}, W0_MAX={W0_MAX}\n")
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

# ── Axis-column outer vertex check ────────────────────────────────────────────
print("\n── Axis-column outer vertex check (nose and tail, y must be 0) ──────────")
nodes, coords, _ = gmsh.model.mesh.getNodes()
coords = np.array(coords).reshape(-1, 3)

x_nose_expected = []
x_tail_expected = []
for k in range(1, N_LAY + 1):
    hk = H1 * (RATIO**k - 1.0) / (RATIO - 1.0)
    x_nose_expected.append(A + hk)
    x_tail_expected.append(-A - hk)

found_nose, found_tail = [], []
for i in range(len(nodes)):
    x, y = coords[i, 0], coords[i, 1]
    if abs(y) < H1 * 0.5:
        for xe in x_nose_expected:
            if abs(x - xe) < H1 * 0.1:
                found_nose.append((x, y)); break
        for xe in x_tail_expected:
            if abs(x - xe) < H1 * 0.1:
                found_tail.append((x, y)); break

ok = True
for label, found in (("nose", found_nose), ("tail", found_tail)):
    found.sort()
    print(f"\n  {label}:")
    for x, y in found:
        flag = "" if abs(y) < 1e-10 else "  <- OFF AXIS"
        print(f"    x={x:.6f}  y={y:.10f}{flag}")
        if abs(y) > 1e-10:
            ok = False
    if not found:
        print(f"    (no axis-column outer vertices matched — check tolerances)")

if ok and (found_nose or found_tail):
    print(f"\n  All axis-column outer vertices on y=0 OK")

# ── Physical groups ────────────────────────────────────────────────────────────
def pg(dim, tags, name):
    t = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, t, name)

pg(3, [vol],                                              "fluid")
pg(2, [sf],                                               "front")
pg(2, [s_top],                                            "back")
pg(2, [curve_to_lat[arc_profile]],                        "body")
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
OF_BASHRC = "/opt/openfoam13/etc/bashrc"

def of_run(cmd, **kwargs):
    return subprocess.run(
        ["bash", "-c", f"source {OF_BASHRC} && {cmd}"],
        **kwargs
    )

sys_dir   = os.path.join(CASE_DIR, "system")
zero_dir  = os.path.join(CASE_DIR, "0")
const_dir = os.path.join(CASE_DIR, "constant")
for d in (sys_dir, zero_dir, const_dir):
    os.makedirs(d, exist_ok=True)

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
        "writeFormat     binary;\n\n"
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
        "    nNonOrthogonalCorrectors 2;\n"
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

# ── foamRun ───────────────────────────────────────────────────────────────────
if run_of:
    print("\n── potentialFoam " + "─" * 57)
    of_run(f"potentialFoam -case {CASE_DIR} -noFunctionObjects")

    print("\n── foamRun " + "─" * 63)
    r = of_run(f"foamRun -case {CASE_DIR}")
    if r.returncode != 0:
        print("\nfoamRun failed — check OpenFOAM sourcing and log above.")
        sys.exit(r.returncode)
    open(os.path.join(CASE_DIR, "case.foam"), "w").close()

print("\nDone — big-Re double-corner mesh in", CASE_DIR)
if not run_of:
    print("(foamRun skipped — pass --of to run the solver)")
