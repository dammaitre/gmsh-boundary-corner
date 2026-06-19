#!/usr/bin/env python3
"""
Standalone reproducer for the BoundaryDoubleCorner degenerate-hex bug.
=======================================================================
Self-contained: does NOT import anything from core/ or geom/. Recreates the
exact same 2D mesh + extrusion problem AEROD's pipeline builds (geom/ellipse.py
profile, theta=10deg, default V0, k=2 refinement level with r=sqrt(2)), using
only the gmsh Python API against the BoundaryDoubleCorner fork build.

Background (see AEROD mesh-quality investigation, 2026-06-19):
  At AEROD's k=2 refinement level (mesh_scale=1/sqrt(2)^2=0.5), ~972 of the
  generated 3D hexahedra have a node tag repeated within their own 8-node
  connectivity list (i.e. a structurally degenerate cell — one vertex used
  twice instead of two distinct corners). gmshToFoam then fails to convert
  the mesh (segfault or FATAL ERROR depending on geometry).

  Root-cause tracing ruled out:
    - the body curve's own 1D mesh (clean, no near-duplicate points)
    - BoundaryDoubleCornerField::operator()'s background-size discontinuity
      at S_corner (patched it, zero effect on the bug)
    - BoundaryDoubleCornerField::buildForFace's own grid construction
      (column-to-column spacing at the exact defect location is normal,
      ~0.27 m, nothing near-zero)

  The degenerate hexes were traced to a vertex CONFLATION: e.g. the quad
  that should span (grid[18][0], grid[19][0], grid[19][1], grid[18][1])
  instead gets grid[19][1] substituted in place of grid[19][0], duplicating
  it. Both vertices are real, distinct, well-separated MVertex objects at
  construction time in buildForFace — so the conflation happens downstream,
  most likely in GMSH's generic extrude+recombine 3D mesher (triggered by
  the BL quads' extreme aspect ratio: ~0.27 m tangential vs ~2e-4 m radial
  at k=2, i.e. ~1300:1), not in the BoundaryDoubleCorner field code itself.

Usage:
    python3 toolbox/repro_bc_degenerate.py            # k=2 (fails)
    python3 toolbox/repro_bc_degenerate.py --k 1       # k=1 (mostly OK, 6 cells)
    python3 toolbox/repro_bc_degenerate.py --k 0       # k=0 (clean)
"""

import argparse
import math
import os
import sys

import numpy as np

GMSH_FORK_LIB = "/home/dam/DEV/gmsh-boundary-corner/build/libgmsh.so"


def ensure_gmsh_lib():
    if os.path.exists(GMSH_FORK_LIB):
        os.environ["GMSH_LIB"] = GMSH_FORK_LIB


# --- geom/ellipse.py, inlined (no AEROD import) ---------------------------

V0_DEFAULT = 200_000.0  # m^3
U_INF_DEFAULT = 44.44   # m/s
NU = 1.5e-5             # m^2/s


_TRESHOLD_5BLNS = (np.sqrt(5) - 1) / 4


def ellipse_compute_R0(theta, V0=None):
    Vc0 = (V0 if V0 is not None else V0_DEFAULT) ** (1 / 3)
    alpha = 1 - 2 * np.tan(theta) ** 2
    if np.tan(theta) <= _TRESHOLD_5BLNS:
        beta = 8 * np.tan(theta) ** 4 - 8 * np.tan(theta) ** 2 + 1
        return Vc0 * (3 / (4 * np.pi * (1 + 2 * alpha ** 3 + 2 * beta ** 3))) ** (1 / 3)
    return Vc0 * (3 / (4 * np.pi * (1 + 2 * alpha ** 3))) ** (1 / 3)


def ellipse_lim(theta, V0=None):
    return 2 * ellipse_compute_R0(theta, V0) / np.tan(theta)


def ellipse_profile(x, theta, V0=None):
    b = ellipse_compute_R0(theta, V0)
    a = b / np.tan(theta)
    return b * np.sqrt(np.maximum(1 - (x - a) ** 2 / a ** 2, 0.0))


def sample_profile(theta_rad, V0=None, n_pts=200, n_fine=20_000):
    """Cosine arc-length-clustered profile sampling — mirrors
    core/phase_geometry.py's geometry_step exactly."""
    L = ellipse_lim(theta_rad, V0)
    x_fine = np.linspace(0, L, n_fine)
    y_fine = np.maximum(ellipse_profile(x_fine, theta_rad, V0), 0.0)

    ds = np.hypot(np.diff(x_fine), np.diff(y_fine))
    s_fine = np.concatenate([[0.0], np.cumsum(ds)])

    t = np.linspace(0, np.pi, n_pts)
    s_target = 0.5 * s_fine[-1] * (1 - np.cos(t))
    x_c = np.interp(s_target, s_fine, x_fine)
    y_c = ellipse_profile(x_c, theta_rad, V0)
    pts = np.column_stack([x_c, y_c])
    pts[:, 1] = np.maximum(pts[:, 1], 0.0)
    pts[0, 1] = 0.0
    pts[-1, 1] = 0.0
    return pts, L


# --- core/constants.py compute_bl_params, inlined --------------------------

def compute_bl_params(u_inf, nu, L, y_plus=30, r_bl=1.2):
    Re = u_inf * L / nu
    Cf = 0.074 / Re ** 0.2
    u_tau = u_inf * np.sqrt(Cf / 2.0)
    h1 = y_plus * nu / u_tau
    delta_bl = 0.37 * L / Re ** 0.2
    n_bl = int(np.ceil(np.log(delta_bl / h1 * (r_bl - 1.0) + 1.0) / np.log(r_bl)))
    h_total = h1 * (r_bl ** n_bl - 1.0) / (r_bl - 1.0)
    return h1, n_bl, r_bl, h_total


# --- degenerate-hex scan ----------------------------------------------------

def scan_degenerate(gmsh):
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = {t: np.array(node_coords[3 * i:3 * i + 3]) for i, t in enumerate(node_tags)}
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim=3)
    degenerate = []
    for et, etags, ents in zip(elem_types, elem_tags, elem_node_tags):
        _, _, _, nnodes, *_ = gmsh.model.mesh.getElementProperties(et)
        if nnodes != 8:
            continue
        ents = np.array(ents).reshape(-1, nnodes)
        for i, tag in enumerate(etags):
            nds = ents[i]
            if len(set(nds.tolist())) < nnodes:
                p = np.array([coords[n] for n in nds])
                degenerate.append((int(tag), nds.tolist(), p[:, 0].mean(), p[:, 1].mean()))
    return degenerate


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=2, choices=[0, 1, 2],
                     help="Refinement level (mesh_scale = 1/sqrt(2)^k). Default 2 (reproduces the bug).")
    ap.add_argument("--theta-deg", type=float, default=10.0)
    ap.add_argument("--gui", action="store_true", help="Open the result in the gmsh GUI after meshing.")
    args = ap.parse_args()

    ensure_gmsh_lib()
    import gmsh

    theta_rad = math.radians(args.theta_deg)
    pts, L = sample_profile(theta_rad, V0=V0_DEFAULT)
    pts = np.maximum(pts, 0.0)

    x_upstream_rel, x_downstream_rel, domain_y_rel = 10.0, 20.0, 25.0
    x_min = -x_upstream_rel * L
    x_max = x_downstream_rel * L
    y_max = domain_y_rel * L

    r = math.sqrt(2)
    mesh_scale = 1.0 / (r ** args.k)
    lc_far = y_max * (20.0 / 150.0) * mesh_scale

    h1_ref, n_bl_ref, r_bl, h_total_ref = compute_bl_params(U_INF_DEFAULT, NU, L)
    if mesh_scale != 1.0:
        h1 = h1_ref * mesh_scale
        n_bl = int(np.ceil(np.log(h_total_ref / h1 * (r_bl - 1.0) + 1.0) / np.log(r_bl)))
    else:
        h1, n_bl = h1_ref, n_bl_ref
    lc_wall = h1 * r_bl ** (n_bl - 1)

    print(f"theta={args.theta_deg} deg  L={L:.4f} m  k={args.k}  mesh_scale={mesh_scale:.4f}")
    print(f"h1={h1:.4e} m  n_bl={n_bl}  r_bl={r_bl}  lc_wall={lc_wall:.4f} m  lc_far={lc_far:.2f} m")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.ExpertMode", 1)
    gmsh.model.add("repro")

    x_nose, x_tail = float(pts[0, 0]), float(pts[-1, 0])
    profile_ptags = [gmsh.model.geo.addPoint(x, y, 0.0, lc_wall) for x, y in pts]
    ptag_nose, ptag_tail = profile_ptags[0], profile_ptags[-1]

    P_inlet_axis = gmsh.model.geo.addPoint(x_min, 0.0, 0.0, lc_far)
    P_inlet_top = gmsh.model.geo.addPoint(x_min, y_max, 0.0, lc_far)
    P_farfield_top = gmsh.model.geo.addPoint(x_max, y_max, 0.0, lc_far)
    P_outlet_axis = gmsh.model.geo.addPoint(x_max, 0.0, 0.0, lc_far)

    CURVE_BODY = gmsh.model.geo.addSpline(profile_ptags)
    CURVE_AXIS_FWD = gmsh.model.geo.addLine(P_inlet_axis, ptag_nose)
    CURVE_AXIS_AFT = gmsh.model.geo.addLine(P_outlet_axis, ptag_tail)
    CURVE_INLET = gmsh.model.geo.addLine(P_inlet_axis, P_inlet_top)
    CURVE_FARFIELD = gmsh.model.geo.addLine(P_inlet_top, P_farfield_top)
    CURVE_OUTLET = gmsh.model.geo.addLine(P_farfield_top, P_outlet_axis)

    loop = gmsh.model.geo.addCurveLoop([
        -CURVE_AXIS_FWD, CURVE_INLET, CURVE_FARFIELD, CURVE_OUTLET, CURVE_AXIS_AFT, -CURVE_BODY,
    ])
    sf = gmsh.model.geo.addPlaneSurface([loop])

    out = gmsh.model.geo.extrude([(2, sf)], 0.0, 0.0, 1.0, numElements=[1], recombine=True)
    gmsh.model.geo.synchronize()

    bdc = gmsh.model.mesh.field.add("BoundaryDoubleCorner")
    gmsh.model.mesh.field.setNumbers(bdc, "CurvesList", [CURVE_BODY])
    gmsh.model.mesh.field.setNumbers(bdc, "NosePoint", [x_nose, 0.0])
    gmsh.model.mesh.field.setNumbers(bdc, "TailPoint", [x_tail, 0.0])
    gmsh.model.mesh.field.setNumber(bdc, "Size", h1)
    gmsh.model.mesh.field.setNumber(bdc, "Ratio", r_bl)
    gmsh.model.mesh.field.setNumber(bdc, "NbLayers", n_bl)
    gmsh.model.mesh.field.setNumber(bdc, "NbCornerColumns", 20)
    gmsh.model.mesh.field.setNumber(bdc, "MaxColumnWidth", lc_wall)
    gmsh.model.mesh.field.setNumber(bdc, "ColWidth", 0.1 * lc_wall)
    gmsh.model.mesh.field.setAsBackgroundMesh(bdc)
    gmsh.option.setNumber("Mesh.BoundaryCornerField", bdc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc_far)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h1 * 0.5)

    gmsh.model.mesh.generate(3)

    degenerate = scan_degenerate(gmsh)
    print(f"\n=== RESULT: {len(degenerate)} degenerate hexahedra (repeated node tag) ===")
    for tag, nds, xc, yc in degenerate[:5]:
        print(f"  elem {tag}: nodes={nds}  centroid=({xc:.4f},{yc:.4f})")
    if len(degenerate) > 5:
        print(f"  ... and {len(degenerate) - 5} more")

    if args.gui:
        gmsh.fltk.run()
    gmsh.finalize()

    sys.exit(1 if degenerate else 0)


if __name__ == "__main__":
    main()
