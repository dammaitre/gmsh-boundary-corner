"""
standard3d.py — pipeline-scale BoundaryDoubleCorner case, full closed 3D
revolution mesh (N_WEDGES wedges covering the full 360 degrees around the
x-axis, the symmetry axis). Same body geometry and BL parameters as
standard.py. Unlike a plain OpenFOAM 2-face `wedge` case (a single thin
angular slice with `wedge` boundary conditions on its two flat faces), this
builds the COMPLETE solid of revolution: N_WEDGES volumes chained end to end
around the axis, closing back onto itself. There is no `wedge`/`axis` patch
at all — the mesh is topologically closed around x, so `body`/`inlet`/
`outlet`/`top` are the only boundary patches.

The body profile is the same single nose-to-tail spline as standard.py,
covered by one BoundaryDoubleCornerField. gmsh's geo-kernel `revolve` used to
refuse to extrude a single non-straight curve whose *both* endpoints lie on
the rotation axis (see the historical FIXME that was in Geo/Geo.cpp
ExtrudeCurve: "the resulting surface would have 2 bounding edges (the axis
and the curve); we cannot handle this case") — exactly this profile's shape
(nose and tail both on y=0). That FIXME has been fixed: ExtrudeCurve now
builds a genuine 2-curve ("digon") lateral surface — the source curve and
its rotated copy — for this case, distinguishing it from a curve that lies
entirely on the axis (still correctly produces no surface). See Geo/Geo.cpp
(ExtrudeCurve) and Geo/GeoInterpolation.cpp (InterpolateRuledSurface).

Building the full revolution as ONE `revolve` call over the full 2*pi angle
does NOT work with this kernel: a full-turn rotation maps every point back
onto itself (regardless of whether it's on the axis), so the "whole curve is
axis-invariant" check that distinguishes "no surface" from "digon surface"
sees every curve as invariant and produces nothing. Instead this script
CHAINS N_WEDGES individual revolve() calls, each by 360/N_WEDGES degrees,
each one revolving the *previous* wedge's outer cap face (its mesh is a
sweep of the source face's mesh, so it already carries the BoundaryDouble-
Corner structure — no need to redefine the field per wedge). The last
wedge's cap face merges back (via gmsh's automatic duplicate-entity
coherence) onto the very first source face, closing the loop with no seam.

Usage:
    python standard3d.py           # headless, mesh + checkMesh only
    python standard3d.py --gui     # open the mesh in the gmsh GUI
    python standard3d.py --of      # also run foamRun after checkMesh
"""

import sys, os, math, re, subprocess, shutil
import numpy as np

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")

import gmsh

# ── Geometry ──────────────────────────────────────────────────────────────────
A = 100.0   # ellipse semi-major axis (axial, x)  [m]
B =  27.0   # ellipse semi-minor axis (radial, y) [m]
L = 2 * A   # body length (nose-to-tail chord), the reference length below

# Flow goes from -x (inlet, tail side) to +x (outlet, nose side):
#   upstream extent (tail → inlet)   = 10*L
#   downstream extent (nose → outlet) = 20*L
#   domain height (axis → top)        = 25*L
X_UP   = A + 10 * L   # inlet x  = -X_UP
X_DOWN = A + 20 * L   # outlet x = +X_DOWN
Y_FAR  = 25 * L       # domain height along y

# ── Scatter sampling ──────────────────────────────────────────────────────────
N_SCATTER = 200  # interior sample points on the spline

# ── Mesh sizes ────────────────────────────────────────────────────────────────
_BIGDC_LC_FAR_RATIO = 20.0 / 150.0
LC_FAR  = Y_FAR * _BIGDC_LC_FAR_RATIO
LC_BODY =  0.1
LC_NOSE =  0.01

# ── BoundaryDoubleCorner parameters ──────────────────────────────────────────
H1        = 5.8e-4  # first BL layer height  (y+ ≈ 50, wall-function regime)
RATIO     = 1.20    # BL growth ratio
N_LAY     = 29      # number of BL layers — outer cell ≈ LC_BODY for smooth structured/Delaunay interface
N_COLS_LEN    = 20  # geometric-length columns (outer part of corner zone)
N_COLS_HEIGHT = 100   # height-blend columns (inner part, at the axis), constant ColWidth
COL_WIDTH = 0.01    # arc-length of innermost column at each corner
W0_MAX    = 0.25    # MaxColumnWidth — must be > LC_BODY to avoid non-ortho at structured/Delaunay interface

# ── Revolution / wedge ────────────────────────────────────────────────────────
N_WEDGES    = 24
WEDGE_ANGLE = 2 * math.pi / N_WEDGES   # one slice of a full 24-wedge revolution

# ── Output ────────────────────────────────────────────────────────────────────
CASE_DIR = os.path.join(_root, "testing", "standard3d_of")
MAX_TIME = 200

# ─────────────────────────────────────────────────────────────────────────────
gui    = "--gui" in sys.argv
run_of = "--of"  in sys.argv
gmsh.initialize(["gmsh", "-nopopup"])
gmsh.model.add("standard3d")

# ── Points ────────────────────────────────────────────────────────────────────
p_nose   = gmsh.model.geo.addPoint(  A,       0,     0, LC_BODY)
p_tail   = gmsh.model.geo.addPoint( -A,       0,     0, LC_BODY)
p_ax_r   = gmsh.model.geo.addPoint(  X_DOWN,  0,     0, LC_FAR)
p_ax_l   = gmsh.model.geo.addPoint( -X_UP,    0,     0, LC_FAR)
p_far_tr = gmsh.model.geo.addPoint(  X_DOWN,  Y_FAR, 0, LC_FAR)
p_far_tl = gmsh.model.geo.addPoint( -X_UP,    Y_FAR, 0, LC_FAR)

# ── Single spline: nose → top → tail  (1-cos distribution in x, planar y(x)) ─
js = np.arange(1, N_SCATTER + 1)
s  = (1.0 - np.cos(np.pi * js / (N_SCATTER + 1))) / 2.0
xs = A - 2.0 * A * s   # x: A (nose) → -A (tail), clustered near both ends
interior = []
for x in xs:
    y = B * math.sqrt(max(0.0, 1.0 - (x / A) ** 2))
    interior.append(gmsh.model.geo.addPoint(x, y, 0, LC_BODY))

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

# ── 3D revolution extrusion (geometry, before meshing) ────────────────────────
# Chain N_WEDGES individual revolve() calls of WEDGE_ANGLE = 360/N_WEDGES each,
# every one revolving the PREVIOUS wedge's outer cap face — see module
# docstring for why a single full-2*pi revolve() call doesn't work here, and
# why chaining doesn't need the BoundaryDoubleCorner field redefined per wedge.
#   numElements=[1] : one circumferential cell layer per wedge
#   recombine=True  : quads → hexes, tris → prisms
# revolve output: [0]=(2,cap), [1]=(3,vol), [2..N]=lateral surfaces of that wedge.
#
# Each wedge's own copies of arc_profile/l_right/l_top/l_left get fresh curve
# tags (only the axis points/lines p_nose, p_tail, p_ax_r, p_ax_l, l_ax_r,
# l_ax_l are invariant under rotation and shared by every wedge), so lateral
# surfaces are classified by ROLE via which axis point(s) each curve touches
# (NOT by curve type string — gmsh reports our addSpline curve's type as
# "Nurb", not "Spline"/"BSpline", so a type-string check silently misses it):
#   profile : endpoints are exactly {p_nose, p_tail} (bows off-axis between)
#   right   : touches p_ax_r
#   left    : touches p_ax_l
#   top     : touches neither (both endpoints off-axis, wedge-specific)

def classify_curve(ctag):
    pts = {abs(p) for _, p in gmsh.model.getBoundary([(1, ctag)], oriented=False)}
    if pts == {p_nose, p_tail}:
        return "profile"
    if p_ax_r in pts:
        return "right"
    if p_ax_l in pts:
        return "left"
    return "top"

lat_by_role = {"profile": [], "right": [], "top": [], "left": []}
vols = []
sf_cur = sf

for k in range(N_WEDGES):
    out = gmsh.model.geo.revolve(
        [(2, sf_cur)], 0, 0, 0, 1, 0, 0, WEDGE_ANGLE,
        numElements=[1], recombine=True
    )
    cap_k = out[0][1]
    vol_k = out[1][1]
    vols.append(vol_k)

    gmsh.model.geo.synchronize()  # needed before getBoundary()/getType() below

    # Map this wedge's own boundary curves (of sf_cur) → lateral surfaces.
    # The revolve output order does NOT match getBoundary order, so identify
    # each lateral surface by checking which curve of sf_cur appears in it.
    cur_curves = {abs(c) for _, c in gmsh.model.getBoundary([(2, sf_cur)], oriented=True)}
    for _, s in out[2:]:
        for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
            ctag = abs(c)
            if ctag in cur_curves:
                lat_by_role[classify_curve(ctag)].append(s)
                break

    sf_cur = cap_k

gmsh.model.geo.synchronize()
vol = vols[0]

# ── BoundaryDoubleCorner field ────────────────────────────────────────────────
bdc = gmsh.model.mesh.field.add("BoundaryDoubleCorner")
gmsh.model.mesh.field.setNumbers(bdc, "CurvesList",      [arc_profile])
gmsh.model.mesh.field.setNumbers(bdc, "NosePoint",       [A,  0.0])
gmsh.model.mesh.field.setNumbers(bdc, "TailPoint",       [-A, 0.0])
gmsh.model.mesh.field.setNumber (bdc, "Size",            H1)
gmsh.model.mesh.field.setNumber (bdc, "Ratio",           RATIO)
gmsh.model.mesh.field.setNumber (bdc, "NbLayers",        N_LAY)
gmsh.model.mesh.field.setNumber (bdc, "NbLengthControl", N_COLS_LEN)
gmsh.model.mesh.field.setNumber (bdc, "NbHeightControl", N_COLS_HEIGHT)
gmsh.model.mesh.field.setNumber (bdc, "MaxColumnWidth",  W0_MAX)
gmsh.model.mesh.field.setNumber (bdc, "ColWidth",        COL_WIDTH)

gmsh.model.mesh.field.setAsBackgroundMesh(bdc)
gmsh.option.setNumber("Mesh.BoundaryCornerField",     bdc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

# ── Generate mesh ─────────────────────────────────────────────────────────────
print(f"\nBuilding spline with {N_SCATTER} interior scatter points.")
print(f"Ellipse: {2*A}m × {2*B}m  (L={L}m)  |  BL: h1={H1:.0e}, {N_LAY} layers, r={RATIO}")
print(f"ColWidth={COL_WIDTH:.0e}, NbLengthControl={N_COLS_LEN}, "
      f"NbHeightControl={N_COLS_HEIGHT}, W0_MAX={W0_MAX}")
print(f"Domain: upstream={X_UP-A:.0f}m (={ (X_UP-A)/L:.1f}L), "
      f"downstream={X_DOWN-A:.0f}m (={ (X_DOWN-A)/L:.1f}L), "
      f"height={Y_FAR:.0f}m (={Y_FAR/L:.1f}L), LC_FAR={LC_FAR:.1f}m")
print(f"Revolution: {N_WEDGES} wedges of {math.degrees(WEDGE_ANGLE):.2f} deg each (full 360 deg, closed)\n")
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

print(f"\n  NbLengthControl={N_COLS_LEN}, NbHeightControl={N_COLS_HEIGHT} (each end),  "
      f"NbLayers={N_LAY},  W0_MAX={W0_MAX}")
print(f"  total quads : {n_quads}")

# ── Physical groups (OpenFOAM patch names) ────────────────────────────────────
def pg(dim, tags, name):
    t = gmsh.model.addPhysicalGroup(dim, tags)
    gmsh.model.setPhysicalName(dim, t, name)

# Closed 360-degree revolution: no wedge_front/wedge_back/axis patches — the
# mesh loops back onto itself around x, so body/outlet/top/inlet (each the
# union of all N_WEDGES wedges' lateral surfaces of that role) are the only
# boundary patches.
pg(3, vols,                    "fluid")
pg(2, lat_by_role["profile"],  "body")
pg(2, lat_by_role["right"],    "outlet")
pg(2, lat_by_role["top"],      "top")
pg(2, lat_by_role["left"],     "inlet")

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
        "gradSchemes { default cellLimited Gauss linear 1; }\n\n"
        "divSchemes\n{\n"
        "    default                              none;\n"
        "    div(phi,U)                           Gauss upwind;\n"
        "    div(phi,k)                           Gauss upwind;\n"
        "    div(phi,omega)                       Gauss upwind;\n"
        "    div((nuEff*dev(T(grad(U)))))         Gauss linear;\n"
        "}\n\n"
        "laplacianSchemes { default Gauss linear limited 0.5; }\n\n"
        "interpolationSchemes { default linear; }\n\n"
        "snGradSchemes { default limitedCorrected 0.5; }\n\n"
        "wallDist { method meshWave; }\n"
    )

with open(os.path.join(sys_dir, "fvSolution"), "w") as f:
    f.write(
        foam_header("dictionary", "system", "fvSolution") +
        "solvers\n{\n"
        "    Phi\n    {\n"
        "        solver          GAMG;\n"
        "        smoother        GaussSeidel;\n"
        "        tolerance       1e-7;\n"
        "        relTol          0.01;\n"
        "    }\n"
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
        "    nNonOrthogonalCorrectors 3;\n"
        "    residualControl\n    {\n"
        "        p       1e-4;\n"
        "        U       1e-4;\n"
        "        k       1e-4;\n"
        "        omega   1e-4;\n"
        "    }\n"
        "}\n\n"
        "relaxationFactors\n{\n"
        "    fields      { p 0.2; }\n"
        "    equations   { U 0.5; k 0.1; omega 0.1; }\n"
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

_PATCHES = ("inlet", "outlet", "top", "body")

def _bc(patch, field):
    zeroG  = "{ type zeroGradient; }"
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
# gmshToFoam writes every patch as type=patch. The mesh is a fully closed
# revolution (no axis/wedge patches — every face on y=0 is interior, shared
# between two cells on either side of the axis), so no defaultFaces renaming
# is needed either: every boundary face already belongs to a named patch.
_PATCH_TYPES = {
    "body":   "wall",
    "inlet":  "patch",
    "outlet": "patch",
    "top":    "patch",
}

def fix_of_boundary(boundary_path):
    with open(boundary_path) as f:
        txt = f.read()

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

# No transformPoints step: unlike a single OpenFOAM `wedge` case (which must
# be centered symmetrically about a coordinate plane), this mesh is a fully
# closed revolution with no wedge patches, so there is no angular centering
# requirement.

# ── checkMesh ─────────────────────────────────────────────────────────────────
print("\n── checkMesh " + "─" * 61)
of_run(f"checkMesh -case {CASE_DIR}")

# ── potentialFoam initialization ──────────────────────────────────────────────
if run_of:
    print("\n── potentialFoam " + "─" * 57)
    r = of_run(f"potentialFoam -writePhi -case {CASE_DIR}")
    if r.returncode != 0:
        print("\npotentialFoam failed.")
        sys.exit(r.returncode)

# ── foamRun ───────────────────────────────────────────────────────────────────
if run_of:
    print("\n── foamRun " + "─" * 63)
    r = of_run(f"foamRun -case {CASE_DIR}")
    if r.returncode != 0:
        print("\nfoamRun failed — check OpenFOAM sourcing and log above.")
        sys.exit(r.returncode)
    open(os.path.join(CASE_DIR, "case.foam"), "w").close()

print(f"\nDone — full 3D closed revolution mesh ({N_WEDGES} wedges) in", CASE_DIR)
if not run_of:
    print("(foamRun skipped — pass --of to run the solver)")
