"""
ellipseexplo.py — BoundaryDoubleCorner robustness sweep.

Sweeps 20 ellipses with theta = arctan(B/A) from 10° to 30°.
B = 27 m is constant; A = B / tan(theta) varies.
For each geometry: mesh, gmshToFoam, checkMesh, collect pass/fail + metrics.

Usage:
    python ellipseexplo.py
"""

import sys, os, math, re, subprocess, shutil
import numpy as np

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")
import gmsh

# ── Fixed parameters ───────────────────────────────────────────────────────────
B         = 27.0
N_THETA   = 20
THETA_MIN = 10.0
THETA_MAX = 30.0

N_SCATTER = 200
LC_BODY   = 0.1
LC_FAR    = 20.0
H1        = 5.8e-4
RATIO     = 1.20
N_LAY     = 29
N_COLS    = 20
COL_WIDTH = 0.01
W0_MAX    = 0.25
DZ        = 1.0

CASE_DIR  = os.path.join(_root, "testing", "ellipseexplo_of")
OF_BASHRC = "/opt/openfoam13/etc/bashrc"

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
    with open(path) as f:
        txt = f.read()
    txt = re.sub(r'\bdefaultFaces\b', 'symAxis', txt)
    def _fix(m):
        name  = m.group(1)
        inner = m.group(2)
        pt    = _PATCH_TYPES.get(name, "patch")
        inner = re.sub(r'(type\s+)\w+(;)',         rf'\g<1>{pt}\2',  inner)
        inner = re.sub(r'(physicalType\s+)\w+(;)', rf'\g<1>{pt}\2',  inner)
        return f'{name}\n    {{\n{inner}    }}'
    txt = re.sub(r'(\w+)\s*\n\s*\{\s*\n(.*?)\n\s*\}', _fix, txt, flags=re.DOTALL)
    with open(path, "w") as f:
        f.write(txt)

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

# ── Mesh one ellipse ───────────────────────────────────────────────────────────
def mesh_one(A):
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
                c2l[abs(c)] = s
                break

    bdc = gmsh.model.mesh.field.add("BoundaryDoubleCorner")
    gmsh.model.mesh.field.setNumbers(bdc, "CurvesList",      [arc])
    gmsh.model.mesh.field.setNumbers(bdc, "NosePoint",       [ A,  0.0])
    gmsh.model.mesh.field.setNumbers(bdc, "TailPoint",       [-A,  0.0])
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

    if os.path.isdir(CASE_DIR):
        shutil.rmtree(CASE_DIR)
    os.makedirs(CASE_DIR)
    msh = os.path.join(CASE_DIR, "ellipse.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(msh)
    return msh

# ── Main ───────────────────────────────────────────────────────────────────────
os.system("clear")
print("ellipseexplo — BoundaryDoubleCorner robustness sweep")
print(f"  B = {B} m (constant)   theta = arctan(B/A) in [{THETA_MIN}°, {THETA_MAX}°]   {N_THETA} steps")
print()

thetas  = np.linspace(THETA_MIN, THETA_MAX, N_THETA)
results = []

gmsh.initialize(["gmsh", "-nopopup"])

for i, theta_deg in enumerate(thetas):
    theta_rad = math.radians(theta_deg)
    A         = B / math.tan(theta_rad)
    print(f"[{i+1:2d}/{N_THETA}] theta={theta_deg:5.2f}°  A={A:7.2f} m  B/A={B/A:.4f}", end="  ", flush=True)

    try:
        msh = mesh_one(A)
        write_minimal_of(CASE_DIR)

        r_conv = of_run(f"gmshToFoam {os.path.basename(msh)}", cwd=CASE_DIR)
        if r_conv.returncode != 0:
            print("GMSH2FOAM FAILED")
            results.append((theta_deg, A, False, None, None, "gmshToFoam failed"))
            continue

        bnd = os.path.join(CASE_DIR, "constant", "polyMesh", "boundary")
        if os.path.exists(bnd):
            fix_of_boundary(bnd)

        r_cm   = of_run(f"checkMesh -case {CASE_DIR}")
        cm_txt = r_cm.stdout + r_cm.stderr
        ok, non_ortho, skewness = parse_checkmesh(cm_txt)

        no_s = f"{non_ortho:5.1f}" if non_ortho is not None else "  n/a"
        sk_s = f"{skewness:.4f}"   if skewness  is not None else "    n/a"
        print(f"{'PASS' if ok else 'FAIL'}   non-ortho={no_s}   skew={sk_s}")
        results.append((theta_deg, A, ok, non_ortho, skewness, None))

    except Exception as e:
        print(f"ERROR: {e}")
        results.append((theta_deg, A, False, None, None, str(e)))

gmsh.finalize()

# ── Summary table ──────────────────────────────────────────────────────────────
W = 70
print()
print("=" * W)
print("SUMMARY")
print("=" * W)
print(f"{'theta':>7}  {'A [m]':>8}  {'B/A':>6}  {'status':6}  {'non-ortho':>9}  {'skewness':>8}")
print("-" * W)
n_pass = n_fail = 0
for theta, A, ok, no, sk, err in results:
    status = "PASS" if ok else "FAIL"
    no_s   = f"{no:9.1f}" if no is not None else "       n/a"
    sk_s   = f"{sk:8.4f}" if sk is not None else "      n/a"
    note   = f"  ({err})"  if err             else ""
    print(f"{theta:7.2f}°  {A:8.2f}  {B/A:6.4f}  {status:6}  {no_s}  {sk_s}{note}")
    if ok: n_pass += 1
    else:  n_fail += 1
print("-" * W)
print(f"Passed: {n_pass}/{N_THETA}    Failed: {n_fail}/{N_THETA}")
print()
