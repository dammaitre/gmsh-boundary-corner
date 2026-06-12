"""
2dwedge.py — OpenFOAM 2D axisymmetric wedge mesh for the full half-ellipse domain.

Geometry identical to double.py.  After 2D meshing the surface is extruded
by a 5° rotation around the x-axis (symmetry axis) to produce a thin wedge
volume suitable for OpenFOAM wedge boundary conditions.

Post-processing:
  gmshToFoam  → reads the .msh file into an OpenFOAM case
  fix_boundary → sets correct patch types (wedge/axis/wall) in polyMesh/boundary
  checkMesh   → validates the resulting polyMesh

Run:  python 2dwedge.py           # headless
      python 2dwedge.py --gui     # open 3D mesh in gmsh GUI before OpenFOAM steps
"""

import sys, os, math, re, subprocess, shutil

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

maxTime = 5000

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
CASE_DIR = os.path.join(_root, "testing", "ellipse2d_of")

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

# ── BoundaryCorner field — back (tail) ────────────────────────────────────────
bc_back = gmsh.model.mesh.field.add("BoundaryCorner")
gmsh.model.mesh.field.setNumbers(bc_back, "CurvesList",      [arc_back])
gmsh.model.mesh.field.setNumbers(bc_back, "AxisPoint",       [-A, 0.0])
gmsh.model.mesh.field.setNumber (bc_back, "Size",            H1)
gmsh.model.mesh.field.setNumber (bc_back, "Ratio",           RATIO)
gmsh.model.mesh.field.setNumber (bc_back, "NbLayers",        N_LAY)
gmsh.model.mesh.field.setNumber (bc_back, "NbCornerColumns", N_COLS)
gmsh.model.mesh.field.setNumber (bc_back, "MaxColumnWidth",  W0_MAX)

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
if os.path.isdir(CASE_DIR):
    shutil.rmtree(CASE_DIR)
os.makedirs(CASE_DIR)
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

# ── OpenFOAM case skeleton ────────────────────────────────────────────────────
sys_dir = os.path.join(CASE_DIR, "system")
zero_dir = os.path.join(CASE_DIR, "0")
const_dir = os.path.join(CASE_DIR, "constant")
for d in (sys_dir, zero_dir, const_dir):
    os.makedirs(d)

def _foam_header(cls, location, obj):
    return (
        "FoamFile\n{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        f"    class       {cls};\n"
        f"    location    \"{location}\";\n"
        f"    object      {obj};\n"
        "}\n\n"
    )

# controlDict — foamRun / SIMPLE, 1000 iterations
ctrl_path = os.path.join(sys_dir, "controlDict")
with open(ctrl_path, "w") as f:
    f.write(
        _foam_header("dictionary", "system", "controlDict") +
        "application     foamRun;\n"
        "solver          incompressibleFluid;\n\n"
        "startFrom       startTime;\n"
        "startTime       0;\n"
        "stopAt          endTime;\n"
        f"endTime         {maxTime};\n"
        "deltaT          1;\n\n"
        "writeControl    timeStep;\n"
        "writeInterval   100;\n\n"
        "runTimeModifiable yes;\n"
    )

# k-ω SST freestream values — 1% turbulence intensity, ν_t/ν = 10
_NU       = 1.5e-5
_U_INF    = 44.44
_k_inf    = 1.5 * (_U_INF * 0.01) ** 2        # ≈ 0.296 m²/s²
_omega_inf = _k_inf / (10.0 * _NU)             # ≈ 1973 s⁻¹

# fvSchemes
with open(os.path.join(sys_dir, "fvSchemes"), "w") as f:
    f.write(
        _foam_header("dictionary", "system", "fvSchemes") +
        "ddtSchemes  { default steadyState; }\n\n"
        "gradSchemes { default Gauss linear; }\n\n"
        "divSchemes\n{\n"
        "    default             none;\n"
        "    div(phi,U)          Gauss linearUpwind grad(U);\n"
        "    div(phi,k)          Gauss limitedLinear 1;\n"
        "    div(phi,omega)      Gauss limitedLinear 1;\n"
        "    div((nuEff*dev(T(grad(U))))) Gauss linear;\n"
        "}\n\n"
        "laplacianSchemes { default Gauss linear corrected; }\n\n"
        "interpolationSchemes { default linear; }\n\n"
        "snGradSchemes { default corrected; }\n\n"
        "wallDist { method meshWave; }\n"
    )

# fvSolution — SIMPLE + k-ω SST
with open(os.path.join(sys_dir, "fvSolution"), "w") as f:
    f.write(
        _foam_header("dictionary", "system", "fvSolution") +
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

# momentumTransport — k-ω SST
with open(os.path.join(const_dir, "momentumTransport"), "w") as f:
    f.write(
        _foam_header("dictionary", "constant", "momentumTransport") +
        "simulationType  RAS;\n\n"
        "RAS\n{\n"
        "    model           kOmegaSST;\n"
        "    turbulence      on;\n"
        "    printCoeffs     on;\n"
        "}\n"
    )

# transportProperties — air at ~15 m/s (Re ~ 2e6 on 2m body)
with open(os.path.join(const_dir, "transportProperties"), "w") as f:
    f.write(
        _foam_header("dictionary", "constant", "transportProperties") +
        "viscosityModel  constant;\n"
        "nu              1.5e-5;\n"
    )

# 0/U — uniform inflow along +x, body no-slip, wedge faces
_wedge_U  = "    type    wedge;\n"
_axis_U   = "    type    empty;\n"
with open(os.path.join(zero_dir, "U"), "w") as f:
    f.write(
        _foam_header("volVectorField", "0", "U") +
        "dimensions  [0 1 -1 0 0 0 0];\n"
        "internalField uniform (44.44 0 0);\n\n"
        "boundaryField\n{\n"
        "    inlet       { type fixedValue; value uniform (44.44 0 0); }\n"
        "    outlet      { type zeroGradient; }\n"
        "    top         { type zeroGradient; }\n"
        "    body        { type noSlip; }\n"
        "    wedge_front { type wedge; }\n"
        "    wedge_back  { type wedge; }\n"
        "    axis        { type empty; }\n"
        "}\n"
    )

# 0/p — kinematic pressure
with open(os.path.join(zero_dir, "p"), "w") as f:
    f.write(
        _foam_header("volScalarField", "0", "p") +
        "dimensions  [0 2 -2 0 0 0 0];\n"
        "internalField uniform 0;\n\n"
        "boundaryField\n{\n"
        "    inlet       { type zeroGradient; }\n"
        "    outlet      { type fixedValue; value uniform 0; }\n"
        "    top         { type zeroGradient; }\n"
        "    body        { type zeroGradient; }\n"
        "    wedge_front { type wedge; }\n"
        "    wedge_back  { type wedge; }\n"
        "    axis        { type empty; }\n"
        "}\n"
    )

# 0/k
with open(os.path.join(zero_dir, "k"), "w") as f:
    f.write(
        _foam_header("volScalarField", "0", "k") +
        "dimensions  [0 2 -2 0 0 0 0];\n"
        f"internalField uniform {_k_inf:.4g};\n\n"
        "boundaryField\n{\n"
        f"    inlet       {{ type fixedValue; value uniform {_k_inf:.4g}; }}\n"
        "    outlet      { type zeroGradient; }\n"
        "    top         { type zeroGradient; }\n"
        "    body        { type kqRWallFunction; value uniform 0; }\n"
        "    wedge_front { type wedge; }\n"
        "    wedge_back  { type wedge; }\n"
        "    axis        { type empty; }\n"
        "}\n"
    )

# 0/omega
with open(os.path.join(zero_dir, "omega"), "w") as f:
    f.write(
        _foam_header("volScalarField", "0", "omega") +
        "dimensions  [0 0 -1 0 0 0 0];\n"
        f"internalField uniform {_omega_inf:.4g};\n\n"
        "boundaryField\n{\n"
        f"    inlet       {{ type fixedValue; value uniform {_omega_inf:.4g}; }}\n"
        "    outlet      { type zeroGradient; }\n"
        "    top         { type zeroGradient; }\n"
        "    body        { type omegaWallFunction; value uniform 1; }\n"
        "    wedge_front { type wedge; }\n"
        "    wedge_back  { type wedge; }\n"
        "    axis        { type empty; }\n"
        "}\n"
    )

# 0/nut
with open(os.path.join(zero_dir, "nut"), "w") as f:
    f.write(
        _foam_header("volScalarField", "0", "nut") +
        "dimensions  [0 2 -1 0 0 0 0];\n"
        "internalField uniform 0;\n\n"
        "boundaryField\n{\n"
        "    inlet       { type calculated; value uniform 0; }\n"
        "    outlet      { type calculated; value uniform 0; }\n"
        "    top         { type calculated; value uniform 0; }\n"
        "    body        { type nutkWallFunction; value uniform 0; }\n"
        "    wedge_front { type wedge; }\n"
        "    wedge_back  { type wedge; }\n"
        "    axis        { type empty; }\n"
        "}\n"
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
# BoundaryCornerField generates MTriangle elements for the columns adjacent to
# the axis points.  Under revolve these triangles become prism cells, which
# are the correct topology for cells adjacent to the axis in OpenFOAM.

# Desired types per patch name (after renaming defaultFaces→axis)
_PATCH_TYPES = {
    "wedge_front": "wedge",
    "wedge_back":  "wedge",
    "axis":        "empty",
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

# ── foamRun (SIMPLE, 1000 iterations) ────────────────────────────────────────
print("\n── foamRun " + "─" * 63)
r = subprocess.run(["foamRun", "-case", CASE_DIR])
if r.returncode != 0:
    print("\nfoamRun failed — check OpenFOAM sourcing and log above.")
    sys.exit(r.returncode)

open(os.path.join(CASE_DIR, "case.foam"), "w").close()
print("\nDone — results in", CASE_DIR)
