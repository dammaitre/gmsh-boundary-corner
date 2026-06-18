"""
bcexplo.py — BoundaryDoubleCorner parameter sensitivity for slender ellipses.

Sweeps N_COLS × W0_MAX × COL_WIDTH at theta = 10°, 12°, 14° (B=27 fixed).
Goal: find parameter combinations that drive max non-ortho below 85°.

Usage:
    python3 bcexplo.py
"""

import sys, os, math, re, subprocess, shutil, itertools
import numpy as np

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")
import gmsh

# ── Fixed geometry / BL ────────────────────────────────────────────────────────
B         = 27.0
N_SCATTER = 200
LC_BODY   = 0.1
LC_FAR    = 20.0
H1        = 5.8e-4
RATIO     = 1.20
N_LAY     = 29
DZ        = 1.0

CASE_DIR  = os.path.join(_root, "testing", "bcexplo_of")
OF_BASHRC = "/opt/openfoam13/etc/bashrc"

# ── Parameter grid ─────────────────────────────────────────────────────────────
THETAS    = [10.0, 12.0, 14.0]          # deg
N_COLS_G  = [15, 20, 25, 30]
W0_MAX_G  = [0.15, 0.25, 0.40]
COL_WIDTH_G = [0.005, 0.01, 0.02]

# ── OF helpers ─────────────────────────────────────────────────────────────────
def of_run(cmd, cwd=None):
    return subprocess.run(
        ["bash", "-c", f"source {OF_BASHRC} && {cmd}"],
        cwd=cwd, capture_output=True, text=True
    )

_PATCH_TYPES = {
    "front": "empty", "back": "empty",
    "body": "wall", "symmetry": "symmetry", "symAxis": "symmetry",
    "inlet": "patch", "outlet": "patch", "top": "patch",
}

def fix_of_boundary(path):
    with open(path) as f: txt = f.read()
    txt = re.sub(r'\bdefaultFaces\b', 'symAxis', txt)
    def _fix(m):
        name  = m.group(1)
        inner = m.group(2)
        pt    = _PATCH_TYPES.get(name, "patch")
        inner = re.sub(r'(type\s+)\w+(;)',         rf'\g<1>{pt}\2',  inner)
        inner = re.sub(r'(physicalType\s+)\w+(;)', rf'\g<1>{pt}\2',  inner)
        return f'{name}\n    {{\n{inner}    }}'
    txt = re.sub(r'(\w+)\s*\n\s*\{\s*\n(.*?)\n\s*\}', _fix, txt, flags=re.DOTALL)
    with open(path, "w") as f: f.write(txt)

def write_minimal_of(case_dir):
    sys_dir = os.path.join(case_dir, "system")
    os.makedirs(sys_dir, exist_ok=True)
    with open(os.path.join(sys_dir, "controlDict"), "w") as f:
        f.write(
            "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            "    class dictionary;\n    location \"system\";\n    object controlDict;\n}\n\n"
            "application     foamRun;\n"
            "startFrom       startTime;\nstartTime       0;\n"
            "stopAt          endTime;\nendTime         1;\ndeltaT          1;\n"
            "writeControl    timeStep;\nwriteInterval   1;\n"
            "writeFormat     ascii;\n"
            "writePrecision  15;\n"  # full float64 precision — default 6 truncates coords → false non-ortho
        )

def parse_checkmesh(txt):
    ok        = "Mesh OK" in txt
    non_ortho = skewness = None
    m = re.search(r'non-orthogonality Max:\s*([\d.e+\-]+)', txt)
    if m: non_ortho = float(m.group(1))
    m = re.search(r'Max skewness\s*=\s*([\d.e+\-]+)', txt)
    if m: skewness = float(m.group(1))
    return ok, non_ortho, skewness

# ── Mesh one case ──────────────────────────────────────────────────────────────
def mesh_one(A, n_cols, w0_max, col_width):
    X_FAR = 3.0 * A
    Y_FAR = 1.5 * A

    gmsh.clear()
    gmsh.model.add("ellipse")
    gmsh.option.setNumber("General.Verbosity", 1)

    p_nose   = gmsh.model.geo.addPoint( A,      0,     0, LC_BODY)
    p_tail   = gmsh.model.geo.addPoint(-A,      0,     0, LC_BODY)
    p_ax_r   = gmsh.model.geo.addPoint( X_FAR,  0,     0, LC_FAR)
    p_ax_l   = gmsh.model.geo.addPoint(-X_FAR,  0,     0, LC_FAR)
    p_far_tr = gmsh.model.geo.addPoint( X_FAR,  Y_FAR, 0, LC_FAR)
    p_far_tl = gmsh.model.geo.addPoint(-X_FAR,  Y_FAR, 0, LC_FAR)

    js = np.arange(1, N_SCATTER + 1)
    s  = (1.0 - np.cos(np.pi * js / (N_SCATTER + 1))) / 2.0
    ts = np.pi * s
    interior = [gmsh.model.geo.addPoint(A * math.cos(t), B * math.sin(t), 0, LC_BODY)
                for t in ts]

    arc     = gmsh.model.geo.addSpline([p_nose] + interior + [p_tail])
    l_ax_r  = gmsh.model.geo.addLine(p_ax_r,   p_nose)
    l_right = gmsh.model.geo.addLine(p_ax_r,   p_far_tr)
    l_top   = gmsh.model.geo.addLine(p_far_tr, p_far_tl)
    l_left  = gmsh.model.geo.addLine(p_far_tl, p_ax_l)
    l_ax_l  = gmsh.model.geo.addLine(p_ax_l,   p_tail)

    cl = gmsh.model.geo.addCurveLoop([-l_ax_r, l_right, l_top, l_left, l_ax_l, -arc])
    sf = gmsh.model.geo.addPlaneSurface([cl])
    out   = gmsh.model.geo.extrude([(2, sf)], 0, 0, DZ, numElements=[1], recombine=True)
    s_top = out[0][1]
    vol   = out[1][1]
    gmsh.model.geo.synchronize()

    orig = {abs(c) for _, c in gmsh.model.getBoundary([(2, sf)], oriented=True)}
    c2l  = {}
    for _, s in out[2:]:
        for _, c in gmsh.model.getBoundary([(2, s)], oriented=False):
            if abs(c) in orig:
                c2l[abs(c)] = s; break

    bdc = gmsh.model.mesh.field.add("BoundaryDoubleCorner")
    gmsh.model.mesh.field.setNumbers(bdc, "CurvesList",      [arc])
    gmsh.model.mesh.field.setNumbers(bdc, "NosePoint",       [ A,  0.0])
    gmsh.model.mesh.field.setNumbers(bdc, "TailPoint",       [-A,  0.0])
    gmsh.model.mesh.field.setNumber (bdc, "Size",            H1)
    gmsh.model.mesh.field.setNumber (bdc, "Ratio",           RATIO)
    gmsh.model.mesh.field.setNumber (bdc, "NbLayers",        N_LAY)
    gmsh.model.mesh.field.setNumber (bdc, "NbCornerColumns", n_cols)
    gmsh.model.mesh.field.setNumber (bdc, "MaxColumnWidth",  w0_max)
    gmsh.model.mesh.field.setNumber (bdc, "ColWidth",        col_width)
    gmsh.model.mesh.field.setAsBackgroundMesh(bdc)
    gmsh.option.setNumber("Mesh.BoundaryCornerField",     bdc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1 * 0.5)

    gmsh.model.mesh.generate(3)

    def pg(dim, tags, name):
        t = gmsh.model.addPhysicalGroup(dim, tags)
        gmsh.model.setPhysicalName(dim, t, name)

    pg(3, [vol],                           "fluid")
    pg(2, [sf],                            "front")
    pg(2, [s_top],                         "back")
    pg(2, [c2l[arc]],                      "body")
    pg(2, [c2l[l_right]],                  "outlet")
    pg(2, [c2l[l_top]],                    "top")
    pg(2, [c2l[l_left]],                   "inlet")
    pg(2, [c2l[l_ax_r], c2l[l_ax_l]],     "symmetry")

    if os.path.isdir(CASE_DIR): shutil.rmtree(CASE_DIR)
    os.makedirs(CASE_DIR)
    msh = os.path.join(CASE_DIR, "ellipse.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(msh)
    return msh

def run_checkmesh(msh):
    write_minimal_of(CASE_DIR)
    r = of_run(f"gmshToFoam {os.path.basename(msh)}", cwd=CASE_DIR)
    if r.returncode != 0:
        return False, None, None, "gmshToFoam failed"
    bnd = os.path.join(CASE_DIR, "constant", "polyMesh", "boundary")
    if os.path.exists(bnd): fix_of_boundary(bnd)
    r2    = of_run(f"checkMesh -case {CASE_DIR}")
    ok, no, sk = parse_checkmesh(r2.stdout + r2.stderr)
    return ok, no, sk, None

# ── Main ───────────────────────────────────────────────────────────────────────
os.system("clear")
print("bcexplo — BoundaryDoubleCorner parameter sensitivity, slender ellipses")
print(f"  B = {B} m   theta ∈ {THETAS}°")
print(f"  N_COLS   : {N_COLS_G}")
print(f"  W0_MAX   : {W0_MAX_G}")
print(f"  COL_WIDTH: {COL_WIDTH_G}")
print()

# Phase 1: sweep N_COLS × W0_MAX at fixed COL_WIDTH=0.01 for all theta
# Phase 2: sweep COL_WIDTH at best (N_COLS, W0_MAX) from phase 1
# This keeps the run count manageable.

results = []   # (theta, n_cols, w0_max, col_width, ok, no, sk, err)

gmsh.initialize(["gmsh", "-nopopup"])

# ── Phase 1 : N_COLS × W0_MAX, COL_WIDTH=0.01 ────────────────────────────────
COL_WIDTH_FIX = 0.01
combos_p1 = list(itertools.product(THETAS, N_COLS_G, W0_MAX_G))
total_p1   = len(combos_p1)

print(f"Phase 1 — N_COLS × W0_MAX ({total_p1} cases, COL_WIDTH={COL_WIDTH_FIX})")
print("-" * 70)

for i, (theta_deg, nc, w0) in enumerate(combos_p1):
    A = B / math.tan(math.radians(theta_deg))
    tag = f"[{i+1:2d}/{total_p1}] θ={theta_deg:4.1f}° nc={nc:2d} w0={w0:.2f} cw={COL_WIDTH_FIX:.3f}"
    print(tag, end="  ", flush=True)
    try:
        msh        = mesh_one(A, nc, w0, COL_WIDTH_FIX)
        ok, no, sk, err = run_checkmesh(msh)
        no_s = f"{no:5.1f}" if no is not None else "  n/a"
        sk_s = f"{sk:.4f}" if sk is not None else "    n/a"
        flag = " <--" if (no is not None and no > 85) else ""
        print(f"{'PASS' if ok else 'FAIL'}  no={no_s}  sk={sk_s}{flag}")
        results.append((theta_deg, nc, w0, COL_WIDTH_FIX, ok, no, sk, err))
    except Exception as e:
        print(f"ERROR: {e}")
        results.append((theta_deg, nc, w0, COL_WIDTH_FIX, False, None, None, str(e)))

# ── Phase 2 : COL_WIDTH sweep at fixed best N_COLS, W0_MAX ────────────────────
# Pick best (N_COLS, W0_MAX) per theta from phase 1 (min non-ortho), then sweep COL_WIDTH
print()
print(f"Phase 2 — COL_WIDTH sweep at best (N_COLS, W0_MAX) per theta")
print("-" * 70)

for theta_deg in THETAS:
    phase1_theta = [(r[1], r[2], r[5]) for r in results
                    if r[0] == theta_deg and r[5] is not None]
    if not phase1_theta:
        continue
    best_nc, best_w0, _ = min(phase1_theta, key=lambda x: x[2])
    A = B / math.tan(math.radians(theta_deg))

    for j, cw in enumerate(COL_WIDTH_G):
        if cw == COL_WIDTH_FIX:
            continue  # already have this from phase 1
        tag = f"  θ={theta_deg:4.1f}° nc={best_nc:2d} w0={best_w0:.2f} cw={cw:.3f}"
        print(tag, end="  ", flush=True)
        try:
            msh        = mesh_one(A, best_nc, best_w0, cw)
            ok, no, sk, err = run_checkmesh(msh)
            no_s = f"{no:5.1f}" if no is not None else "  n/a"
            sk_s = f"{sk:.4f}" if sk is not None else "    n/a"
            flag = " <--" if (no is not None and no > 85) else ""
            print(f"{'PASS' if ok else 'FAIL'}  no={no_s}  sk={sk_s}{flag}")
            results.append((theta_deg, best_nc, best_w0, cw, ok, no, sk, err))
        except Exception as e:
            print(f"ERROR: {e}")
            results.append((theta_deg, best_nc, best_w0, cw, False, None, None, str(e)))

gmsh.finalize()

# ── Per-theta summary: best parameters ────────────────────────────────────────
W = 72
print()
print("=" * W)
print("SUMMARY — best parameters per theta (lowest non-ortho)")
print("=" * W)
print(f"{'theta':>7}  {'N_COLS':>6}  {'W0_MAX':>6}  {'ColWid':>6}  {'no_max':>6}  {'skew':>7}  status")
print("-" * W)
for theta_deg in THETAS:
    theta_res = [(r[1], r[2], r[3], r[4], r[5], r[6])
                 for r in results if r[0] == theta_deg and r[5] is not None]
    if not theta_res:
        print(f"{theta_deg:7.1f}°  no valid results")
        continue
    nc, w0, cw, ok, no, sk = min(theta_res, key=lambda x: x[4])
    no_s = f"{no:6.1f}" if no is not None else "   n/a"
    sk_s = f"{sk:7.4f}" if sk is not None else "    n/a"
    print(f"{theta_deg:7.1f}°  {nc:6d}  {w0:6.2f}  {cw:6.3f}  {no_s}  {sk_s}  {'PASS' if ok else 'FAIL'}")

# ── Full detail table ──────────────────────────────────────────────────────────
print()
print("=" * W)
print("FULL RESULTS")
print("=" * W)
print(f"{'theta':>7}  {'N_COLS':>6}  {'W0_MAX':>6}  {'ColWid':>6}  {'no_max':>6}  {'skew':>7}  status")
print("-" * W)
for theta_deg in THETAS:
    theta_res = [r for r in results if r[0] == theta_deg]
    theta_res.sort(key=lambda x: (x[5] if x[5] is not None else 999))
    for r in theta_res:
        _, nc, w0, cw, ok, no, sk, err = r
        no_s = f"{no:6.1f}" if no is not None else "   n/a"
        sk_s = f"{sk:7.4f}" if sk is not None else "    n/a"
        note = f"  ({err})" if err else ""
        flag = " <" if (no is not None and no > 85) else ""
        print(f"{theta_deg:7.1f}°  {nc:6d}  {w0:6.2f}  {cw:6.3f}  {no_s}  {sk_s}  {'PASS' if ok else 'FAIL'}{flag}{note}")
    print()
