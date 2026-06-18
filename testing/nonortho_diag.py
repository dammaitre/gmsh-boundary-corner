"""
nonortho_diag.py — characterize severely non-orthogonal faces across the theta transition.

Sweeps theta from 10° to 17°.  For each theta:
  - meshes the case
  - runs checkMesh (ground truth for max non-ortho)
  - reads polyMesh directly, computes non-ortho with correct OF centroid algorithm
  - reports count, positions, and face-type analysis

Key finding: the bad faces are TANGENTIAL faces (between adjacent corner columns) where
the cell-to-cell vector is tangential but the face normal is radial — a topology mismatch
at the corner-to-regular-BL transition.

Usage:
    python3 testing/nonortho_diag.py
"""

import sys, os, math, re, subprocess, shutil
import numpy as np
from collections import Counter, defaultdict

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_root, "api"))
os.environ["GMSH_LIB"] = os.path.join(_root, "build", "libgmsh.so")
import gmsh

# ── Fixed params (same as ellipseexplo) ────────────────────────────────────────
B         = 27.0
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
NO_THRESH = 70.0   # match checkMesh "severely non-orthogonal" threshold

CASE_DIR  = os.path.join(_root, "testing", "nonortho_diag_of")
OF_BASHRC = "/opt/openfoam13/etc/bashrc"

# ── Theta sweep: 10-14 coarse, fine near transition, 15.1-17 coarse ───────────
thetas_deg = sorted(set([
    10.0, 11.0, 12.0, 13.0, 14.0,
    14.5, 14.8, 14.9, 15.0,
    15.1, 15.2, 15.5,
    16.0, 17.0,
]))

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
        name = m.group(1); inner = m.group(2)
        pt   = _PATCH_TYPES.get(name, "patch")
        inner = re.sub(r'(type\s+)\w+(;)',         rf'\g<1>{pt}\2', inner)
        inner = re.sub(r'(physicalType\s+)\w+(;)', rf'\g<1>{pt}\2', inner)
        return f'{name}\n    {{\n{inner}    }}'
    txt = re.sub(r'(\w+)\s*\n\s*\{\s*\n(.*?)\n\s*\}', _fix, txt, flags=re.DOTALL)
    with open(path, "w") as f: f.write(txt)

def write_minimal_of(d):
    s = os.path.join(d, "system"); os.makedirs(s, exist_ok=True)
    with open(os.path.join(s, "controlDict"), "w") as f:
        f.write(
            "FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            "    class dictionary;\n    location \"system\";\n    object controlDict;\n}\n\n"
            "application foamRun;\nstartFrom startTime;\nstartTime 0;\n"
            "stopAt endTime;\nendTime 1;\ndeltaT 1;\n"
            "writeControl timeStep;\nwriteInterval 1;\nwriteFormat ascii;\n"
            "writePrecision 15;\n"  # full float64 precision — default 6 truncates coords → false non-ortho
        )

# ── polyMesh readers ───────────────────────────────────────────────────────────
def _skip_header(txt):
    depth = 0; i = 0; in_h = False
    while i < len(txt):
        if txt[i:i+8] == 'FoamFile': in_h = True
        if in_h:
            if txt[i] == '{': depth += 1
            elif txt[i] == '}':
                depth -= 1
                if depth == 0: return txt[i+1:]
        i += 1
    return txt

def read_points(path):
    with open(path) as f: raw = f.read()
    raw = re.sub(r'//[^\n]*', '', raw); raw = _skip_header(raw)
    m = re.search(r'\b(\d+)\b', raw); n = int(m.group(1))
    nums = re.findall(r'[-+]?\d+\.?\d*(?:[eE][+-]?\d+)?', raw[m.end():])
    return np.array(nums[:n*3], dtype=float).reshape(n, 3)

def read_faces(path):
    with open(path) as f: raw = f.read()
    raw = re.sub(r'//[^\n]*', '', raw); raw = _skip_header(raw)
    return [list(map(int, e[1].split())) for e in re.findall(r'(\d+)\(([\d\s]+)\)', raw)]

def read_int_list(path):
    with open(path) as f: raw = f.read()
    raw = re.sub(r'//[^\n]*', '', raw); raw = _skip_header(raw)
    m = re.search(r'\b(\d+)\b', raw); n = int(m.group(1))
    return list(map(int, re.findall(r'\d+', raw[m.end():])[:n]))

def face_sf(verts, pts):
    """Newell area vector. For internal faces: owner→neighbour direction."""
    p = pts[verts]; a = np.zeros(3); n = len(p)
    for i in range(n):
        j = (i+1) % n
        a[0] += (p[i,1]-p[j,1])*(p[i,2]+p[j,2])
        a[1] += (p[i,2]-p[j,2])*(p[i,0]+p[j,0])
        a[2] += (p[i,0]-p[j,0])*(p[i,1]+p[j,1])
    return a * 0.5

def face_centroid(verts, pts):
    """Area-weighted centroid via triangle fan. Uses SCALAR triangle area."""
    p = pts[verts]; n = len(p)
    if n == 3: return p.mean(axis=0)
    p0 = p[0]; tot_a = 0.0; cent = np.zeros(3)
    for i in range(1, n-1):
        e1 = p[i]-p0; e2 = p[i+1]-p0
        a = np.linalg.norm(np.cross(e1, e2)) / 2.0
        cent += a * (p0+p[i]+p[i+1]) / 3.0; tot_a += a
    return cent/tot_a if tot_a > 1e-15 else p.mean(axis=0)

def compute_cell_centroids(pts, faces, owner, neigh):
    """OpenFOAM two-pass pyramid algorithm."""
    n_int = len(neigh); n_tot = len(faces); n_cells = max(owner)+1
    fc_all = np.array([face_centroid(f, pts) for f in faces])
    # Pass 1: simple face-centroid average
    cc0_s = np.zeros((n_cells,3)); cc0_n = np.zeros(n_cells,int)
    for fi in range(n_tot):
        ow = owner[fi]; cc0_s[ow] += fc_all[fi]; cc0_n[ow] += 1
        if fi < n_int: nb = neigh[fi]; cc0_s[nb] += fc_all[fi]; cc0_n[nb] += 1
    cc0 = cc0_s / cc0_n[:,None]
    # Pass 2: pyramid-volume weighted
    cv = np.zeros((n_cells,3)); cw = np.zeros(n_cells)
    for fi in range(n_tot):
        Sf = face_sf(faces[fi], pts); fc = fc_all[fi]; ow = owner[fi]
        d = fc-cc0[ow]; vol = abs(np.dot(d,Sf))/3.0
        cv[ow] += vol*((3*fc+cc0[ow])/4.0); cw[ow] += vol
        if fi < n_int:
            nb = neigh[fi]; d = fc-cc0[nb]; vol = abs(np.dot(d,Sf))/3.0
            cv[nb] += vol*((3*fc+cc0[nb])/4.0); cw[nb] += vol
    return np.where(cw[:,None] > 0, cv/cw[:,None], cc0), fc_all

def no_deg(Sf, d):
    ms = np.linalg.norm(Sf); md = np.linalg.norm(d)
    if ms < 1e-15 or md < 1e-15: return 0.0
    return math.degrees(math.acos(min(1.0, abs(np.dot(Sf,d)/(ms*md)))))

def analyze_mesh(mesh_dir, A_val):
    """Return (n_bad, bad_faces_list, max_no) where bad_faces_list is
    [(no, face_centroid_xy, owner_cc_xy, neigh_cc_xy), ...]."""
    pts   = read_points(os.path.join(mesh_dir, 'points'))
    faces = read_faces (os.path.join(mesh_dir, 'faces'))
    owner = read_int_list(os.path.join(mesh_dir, 'owner'))
    neigh = read_int_list(os.path.join(mesh_dir, 'neighbour'))
    n_int = len(neigh)

    cc, fc_all = compute_cell_centroids(pts, faces, owner, neigh)

    bad = []
    for fi in range(n_int):
        Sf = face_sf(faces[fi], pts)
        d  = cc[neigh[fi]] - cc[owner[fi]]
        no = no_deg(Sf, d)
        if no > NO_THRESH:
            fc_xy  = fc_all[fi, :2]
            ow_xy  = cc[owner[fi], :2]
            nb_xy  = cc[neigh[fi], :2]
            # body normal at face centroid x
            t_ang  = math.acos(max(-1.0, min(1.0, fc_xy[0]/A_val)))
            n_body = np.array([B*math.cos(t_ang), A_val*math.sin(t_ang)])
            n_body /= np.linalg.norm(n_body)
            t_body = np.array([-n_body[1], n_body[0]])
            # body surface y at this x
            y_body = B * math.sin(t_ang)
            # decompose d into normal/tangential
            d2 = d[:2]
            d_norm = np.dot(d2, n_body)
            d_tang = np.dot(d2, t_body)
            bad.append((no, fc_xy, ow_xy, nb_xy, d_norm, d_tang, y_body))

    bad.sort(key=lambda x: -x[0])
    return len(bad), bad

# ── Mesh builder (same as ellipseexplo) ───────────────────────────────────────
def mesh_one(A):
    X_FAR = 3.0*A; Y_FAR = 1.5*A
    gmsh.clear(); gmsh.model.add("e")
    gmsh.option.setNumber("General.Verbosity", 1)
    p_no=gmsh.model.geo.addPoint( A,     0,     0, LC_BODY)
    p_ta=gmsh.model.geo.addPoint(-A,     0,     0, LC_BODY)
    p_ar=gmsh.model.geo.addPoint( X_FAR, 0,     0, LC_FAR)
    p_al=gmsh.model.geo.addPoint(-X_FAR, 0,     0, LC_FAR)
    p_tr=gmsh.model.geo.addPoint( X_FAR, Y_FAR, 0, LC_FAR)
    p_tl=gmsh.model.geo.addPoint(-X_FAR, Y_FAR, 0, LC_FAR)
    js  = np.arange(1, N_SCATTER+1); s=(1-np.cos(np.pi*js/(N_SCATTER+1)))/2; ts=np.pi*s
    interior=[gmsh.model.geo.addPoint(A*math.cos(t), B*math.sin(t), 0, LC_BODY) for t in ts]
    arc  = gmsh.model.geo.addSpline([p_no]+interior+[p_ta])
    l_ar = gmsh.model.geo.addLine(p_ar, p_no); l_ri=gmsh.model.geo.addLine(p_ar, p_tr)
    l_to = gmsh.model.geo.addLine(p_tr, p_tl); l_le=gmsh.model.geo.addLine(p_tl, p_al)
    l_al = gmsh.model.geo.addLine(p_al, p_ta)
    cl   = gmsh.model.geo.addCurveLoop([-l_ar, l_ri, l_to, l_le, l_al, -arc])
    sf   = gmsh.model.geo.addPlaneSurface([cl])
    out  = gmsh.model.geo.extrude([(2,sf)],0,0,DZ,numElements=[1],recombine=True)
    s_top=out[0][1]; vol=out[1][1]
    gmsh.model.geo.synchronize()
    orig={abs(c) for _,c in gmsh.model.getBoundary([(2,sf)],oriented=True)}
    c2l={}
    for _,s2 in out[2:]:
        for _,c in gmsh.model.getBoundary([(2,s2)],oriented=False):
            if abs(c) in orig: c2l[abs(c)]=s2; break
    bdc=gmsh.model.mesh.field.add("BoundaryDoubleCorner")
    gmsh.model.mesh.field.setNumbers(bdc,"CurvesList",[arc])
    gmsh.model.mesh.field.setNumbers(bdc,"NosePoint",[A,0.0])
    gmsh.model.mesh.field.setNumbers(bdc,"TailPoint",[-A,0.0])
    gmsh.model.mesh.field.setNumber(bdc,"Size",H1); gmsh.model.mesh.field.setNumber(bdc,"Ratio",RATIO)
    gmsh.model.mesh.field.setNumber(bdc,"NbLayers",N_LAY)
    gmsh.model.mesh.field.setNumber(bdc,"NbCornerColumns",N_COLS)
    gmsh.model.mesh.field.setNumber(bdc,"MaxColumnWidth",W0_MAX)
    gmsh.model.mesh.field.setNumber(bdc,"ColWidth",COL_WIDTH)
    gmsh.model.mesh.field.setAsBackgroundMesh(bdc)
    gmsh.option.setNumber("Mesh.BoundaryCornerField",     bdc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", LC_FAR)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", H1*0.5)
    gmsh.model.mesh.generate(3)
    def pg(d,t,n): tt=gmsh.model.addPhysicalGroup(d,t); gmsh.model.setPhysicalName(d,tt,n)
    pg(3,[vol],"fluid"); pg(2,[sf],"front"); pg(2,[s_top],"back")
    pg(2,[c2l[arc]],"body"); pg(2,[c2l[l_ri]],"outlet"); pg(2,[c2l[l_to]],"top")
    pg(2,[c2l[l_le]],"inlet"); pg(2,[c2l[l_ar],c2l[l_al]],"symmetry")
    if os.path.isdir(CASE_DIR): shutil.rmtree(CASE_DIR)
    os.makedirs(CASE_DIR)
    msh = os.path.join(CASE_DIR,"ellipse.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2); gmsh.write(msh)
    return msh

# ── Main sweep ─────────────────────────────────────────────────────────────────
os.system("clear")
print("nonortho_diag — characterization of severe non-ortho across theta transition")
print(f"  B={B}m  threshold={NO_THRESH}°  N_COLS={N_COLS}  W0_MAX={W0_MAX}")
print()

gmsh.initialize(["gmsh","-nopopup"])

all_results = []  # (theta, A, n_bad, cm_max, bad_list)

for theta_deg in thetas_deg:
    A = B / math.tan(math.radians(theta_deg))
    print(f"theta={theta_deg:5.1f}°  A={A:7.2f}m ... ", end="", flush=True)
    try:
        msh = mesh_one(A)
        write_minimal_of(CASE_DIR)
        r = of_run(f"gmshToFoam {os.path.basename(msh)}", cwd=CASE_DIR)
        if r.returncode != 0:
            print("gmshToFoam FAILED"); all_results.append((theta_deg,A,-1,None,[])); continue
        bnd = os.path.join(CASE_DIR,"constant","polyMesh","boundary")
        if os.path.exists(bnd): fix_of_boundary(bnd)
        cm = of_run("checkMesh", cwd=CASE_DIR)
        m = re.search(r'non-orthogonality Max:\s*([\d.e+\-]+)', cm.stdout+cm.stderr)
        cm_max = float(m.group(1)) if m else None
        n_bad, bad = analyze_mesh(os.path.join(CASE_DIR,"constant","polyMesh"), A)
        my_max = bad[0][0] if bad else 0.0
        print(f"checkMesh={cm_max:.2f}°  mine: {n_bad:4d} bad  max={my_max:.2f}°")
        all_results.append((theta_deg, A, n_bad, cm_max, bad))
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback; traceback.print_exc()
        all_results.append((theta_deg, A, -1, None, []))

gmsh.finalize()

# ── Summary ────────────────────────────────────────────────────────────────────
print()
print("=" * 90)
print(f"  {'theta':>7}  {'A [m]':>8}  {'cm_max':>8}  {'n_bad':>6}  {'my_max':>8}"
      f"  {'nose':>5}  {'tail':>5}  {'body':>5}  y_range")
print("-" * 90)
for theta, A, n_bad, cm_max, bad in all_results:
    if n_bad <= 0:
        cms = f"{cm_max:.2f}°" if cm_max else "  n/a"
        print(f"  {theta:7.1f}°  {A:8.2f}  {cms:>8}  {0:6d}  {'0.00°':>8}")
        continue
    my_max = bad[0][0]
    xs  = [b[1][0] for b in bad]; ys = [b[1][1] for b in bad]
    nose_n = sum(1 for x in xs if x >  0.9*A)
    tail_n = sum(1 for x in xs if x < -0.9*A)
    body_n = n_bad - nose_n - tail_n
    cms = f"{cm_max:.2f}°" if cm_max else "  n/a"
    print(f"  {theta:7.1f}°  {A:8.2f}  {cms:>8}  {n_bad:6d}  {my_max:8.2f}°"
          f"  {nose_n:5d}  {tail_n:5d}  {body_n:5d}  y=[{min(ys):.3f},{max(ys):.3f}]")

# ── Detailed face analysis for key theta values ────────────────────────────────
print()
print("=" * 90)
print("DETAILED FACE ANALYSIS — worst 5 faces per theta (internal faces only)")
print("  d_norm = cell-to-cell component along body normal")
print("  d_tang = cell-to-cell component along body tangent")
print("  d_norm >> d_tang = radially adjacent cells (expected for radial faces)")
print("  d_tang >> d_norm = tangentially adjacent cells (TOPOLOGY MISMATCH)")
print()

for theta, A, n_bad, cm_max, bad in all_results:
    if n_bad <= 0 or theta not in (10.0, 13.0, 15.0, 15.1, 16.0):
        continue
    R_nose = B**2 / A
    print(f"── theta={theta:.1f}°  A={A:.1f}m  R_nose={R_nose:.2f}m  n_bad={n_bad} ──")
    print(f"  {'no':>7}  {'x':>9}  {'y':>8}  {'y_body':>8}  {'dist':>6}  {'d_norm':>8}  {'d_tang':>8}  loc")
    for no, fc, ow_cc, nb_cc, d_norm, d_tang, y_body in bad[:5]:
        dist = fc[1] - y_body  # approximate y-distance from body surface
        loc  = 'NOSE' if fc[0] > 0.9*A else ('TAIL' if fc[0] < -0.9*A else 'body')
        print(f"  {no:7.2f}°  {fc[0]:+9.2f}  {fc[1]:8.4f}  {y_body:8.4f}  "
              f"{dist:6.4f}  {d_norm:8.5f}  {d_tang:8.5f}  {loc}")
    print()
