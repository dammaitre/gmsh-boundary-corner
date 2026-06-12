"""
2dplanar.py — OpenFOAM 2D planar mesh for a half-ellipse body at stagnation.

Geometry: half-ellipse (semi-major A along x, semi-minor B along y) inside a
rectangular far-field.  The mesh lives in the y ≥ 0 half-plane (x-axis = symmetry
axis).  The BoundaryCorner field handles the 90° intersection of the body surface
with the axis at nose and tail.

After 2D meshing the surface is extruded 1 cell deep (DZ) along z to produce a
flat volume for OpenFOAM 2D planar cases (empty front/back patches).

Axisymmetric note: revolving this mesh into a wedge does NOT work — use the 2D
planar approach and handle axisymmetry in the solver (2πr integration, etc.).

Usage:
    python 2dplanar.py           # headless
    python 2dplanar.py --gui     # open the mesh in the gmsh GUI before OpenFOAM steps
"""

import sys, os, math, re, subprocess, shutil

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry ──────────────────────────────────────────────────────────────────
A      = 2.0    # ellipse semi-major axis (axial, x)
B      = 1.2    # ellipse semi-minor axis (radial, y)
X_FAR  = 6.0
Y_FAR  = 3.5

# ── Mesh sizes ────────────────────────────────────────────────────────────────
LC_FAR  = 0.30
LC_BODY = 0.06
LC_NOSE = 0.03

# ── BoundaryCorner parameters ─────────────────────────────────────────────────
H1     = 0.012   # first BL layer height
RATIO  = 1.20    # BL growth ratio
N_LAY  = 6       # number of BL layers
N_COLS = 8       # fan columns near axis corner
W0_MAX = LC_BODY # max tangential column width at StartPoint

# ── Extrusion ─────────────────────────────────────────────────────────────────
DZ = 1.0

# ── Output ────────────────────────────────────────────────────────────────────
CASE_DIR = os.path.join(_root, "testing", "ellipse2dplanar_of")
MAX_TIME = 5000

# ─────────────────────────────────────────────────────────────────────────────
gui = "--gui" in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("ellipse2dplanar")

# ── Points ────────────────────────────────────────────────────────────────────
p_origin = gmsh.model.geo.addPoint(  0,      0,     0, LC_BODY)
p_nose   = gmsh.model.geo.addPoint(  A,      0,     0, LC_NOSE)
p_tail   = gmsh.model.geo.addPoint( -A,      0,     0, LC_NOSE)
p_top    = gmsh.model.geo.addPoint(  0,      B,     0, LC_BODY)
p_ax_r   = gmsh.model.geo.addPoint(  X_FAR,  0,     0, LC_FAR)
p_ax_l   = gmsh.model.geo.addPoint( -X_FAR,  0,     0, LC_FAR)
p_far_tr = gmsh.model.geo.addPoint(  X_FAR,  Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint( -X_FAR,  Y_FAR, 0, LC_FAR)

# ── Curves ────────────────────────────────────────────────────────────────────
# Front arc: shoulder (0, B) → nose (A, 0)
arc_front = gmsh.model.geo.addEllipseArc(p_top,  p_origin, p_nose, p_nose)
# Back arc:  tail (-A, 0) → shoulder (0, B)
arc_back  = gmsh.model.geo.addEllipseArc(p_tail, p_origin, p_nose, p_top)

l_ax_r  = gmsh.model.geo.addLine(p_ax_r,   p_nose)
l_right = gmsh.model.geo.addLine(p_ax_r,   p_far_tr)
l_top   = gmsh.model.geo.addLine(p_far_tr, p_far_tl)
l_left  = gmsh.model.geo.addLine(p_far_tl, p_ax_l)
l_ax_l  = gmsh.model.geo.addLine(p_ax_l,   p_tail)

cl = gmsh.model.geo.addCurveLoop(
    [-l_ax_r, l_right, l_top, l_left, l_ax_l, arc_back, arc_front])
sf = gmsh.model.geo.addPlaneSurface([cl])

# ── 3D extrusion (geometry only, before meshing) ──────────────────────────────
# One cell layer along z; recombine produces hexes from quads, prisms from tris.
# out[0] = (2, s_top) — the extruded copy of sf
# out[1] = (3, vol)   — the volume
# out[2:] = lateral surfaces, one per bounding curve of sf
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
gmsh.model.mesh.generate(3)

# ── Physical groups (OpenFOAM patch names) ────────────────────────────────────
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

# k-ω SST freestream values: 1% turbulence intensity, ν_t/ν = 10
_NU        = 1.5e-5
_U_INF     = 44.44
_k_inf     = 1.5 * (_U_INF * 0.01) ** 2
_omega_inf = _k_inf / (10.0 * _NU)

_PATCHES = ("inlet", "outlet", "top", "body", "symmetry", "symAxis", "front", "back")

def _bc(patch, field):
    """Return the boundary condition string for a given patch and field."""
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

    # BoundaryCorner axis-column cells at y=0 land in defaultFaces; rename to
    # symAxis so it gets a separate symmetry patch from the main symmetry strip.
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

# ── foamRun ───────────────────────────────────────────────────────────────────
print("\n── foamRun " + "─" * 63)
r = subprocess.run(["foamRun", "-case", CASE_DIR])
if r.returncode != 0:
    print("\nfoamRun failed — check OpenFOAM sourcing and log above.")
    sys.exit(r.returncode)

open(os.path.join(CASE_DIR, "case.foam"), "w").close()
print("\nDone — results in", CASE_DIR)
