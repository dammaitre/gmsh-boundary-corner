#!/usr/bin/env python3
"""
BoundaryCorner visualization — external axisymmetric CFD flow.

Geometry:
  - Half-ellipse (right half, cut by y-axis at x=0) = the body surface / wall
  - Solid body is INSIDE the ellipse
  - Fluid is OUTSIDE (above and to the right of the profile)
  - Symmetry axis: y = 0 (horizontal)
  - The profile meets y=0 at 90° at the nose = AxisPoint

BoundaryCorner zone:
  - From StartPoint (on profile, y>0) to AxisPoint (nose, y=0)
  - Structured quad fan growing OUTWARD from the body surface (into the fluid)
  - Each column = outward-normal direction; layers = BL thickness distribution
  - The fan sweeps 90° as the outward normal rotates from ~vertical to horizontal
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial import Delaunay

# ── Ellipse geometry ─────────────────────────────────────────────────────────
a, b = 2.0, 1.2          # semi-axes: x (major), y (minor)
theta_start = 55 * np.pi / 180   # angle where BC zone begins on the profile

# ── BoundaryCorner parameters ────────────────────────────────────────────────
N_COL = 6        # nbCornerColumns
N_LAY = 5        # nbLayers
RATIO = 1.3      # geometric ratio (fine near wall, coarse outward)
H_TOTAL = 0.32   # total BL thickness

# ── Derived ──────────────────────────────────────────────────────────────────
x_sp = a * np.cos(theta_start)
y_sp = b * np.sin(theta_start)

def ellipse_y(x):
    return b * np.sqrt(np.maximum(0.0, 1.0 - (x / a) ** 2))

def outward_normal(theta_val):
    """Unit outward normal to ellipse at parametric angle theta."""
    nx = b * np.cos(theta_val)
    ny = a * np.sin(theta_val)
    return nx / np.hypot(nx, ny), ny / np.hypot(nx, ny)

# Layer fractions: 0 = surface (wall), 1 = BL outer edge
k_arr = np.arange(N_LAY + 1, dtype=float)
hk = (RATIO ** k_arr - 1.0) / (RATIO ** N_LAY - 1.0)

# ── Column positions: geometric arc-length progression along the profile ─────
# Arc-length integrand for ellipse: ds = sqrt((a sinθ)² + (b cosθ)²) dθ
theta_fine = np.linspace(theta_start, 0.0, 8000)
integ = np.sqrt((a * np.sin(theta_fine))**2 + (b * np.cos(theta_fine))**2)
dtheta = np.abs(np.diff(theta_fine))
arc_cum = np.zeros(len(theta_fine))
arc_cum[1:] = np.cumsum(dtheta * (integ[:-1] + integ[1:]) / 2)
total_arc = arc_cum[-1]

# Geometric progression of column arc positions:
# widths w_0, w_0*eps, w_0*eps², ..., w_0*eps^(N_COL-1)
# eps < 1 → columns compress toward AxisPoint (nose), matching high curvature
EPS_TAN = 0.68
w_0 = total_arc * (1 - EPS_TAN) / (1 - EPS_TAN ** N_COL)
arc_pos = [0.0]
w = w_0
for _ in range(N_COL):
    arc_pos.append(arc_pos[-1] + w)
    w *= EPS_TAN
arc_pos = np.array(arc_pos)
arc_pos[-1] = total_arc   # pin exactly at AxisPoint

# Map arc positions back to theta via interpolation
theta_cols = np.interp(arc_pos, arc_cum, theta_fine)
xc = a * np.cos(theta_cols)
yc = b * np.sin(theta_cols)

# ── BC fan grid ───────────────────────────────────────────────────────────────
# gx[i, k], gy[i, k]:  column i (along profile), layer k (outward from surface)
gx = np.zeros((N_COL + 1, N_LAY + 1))
gy = np.zeros((N_COL + 1, N_LAY + 1))
for i, th in enumerate(theta_cols):
    xi, yi = xc[i], yc[i]
    nx, ny = outward_normal(th)
    gx[i, :] = xi + hk * H_TOTAL * nx
    gy[i, :] = yi + hk * H_TOTAL * ny

# ── Regular BL to the left of StartPoint (illustrative, 4 columns) ───────────
N_REG = 4
theta_reg = np.linspace(np.pi / 2, theta_start, N_REG + 1)  # uniform spacing ok here
xr = a * np.cos(theta_reg)
yr = b * np.sin(theta_reg)
rx = np.zeros((N_REG + 1, N_LAY + 1))
ry = np.zeros((N_REG + 1, N_LAY + 1))
for i, th in enumerate(theta_reg):
    nx, ny = outward_normal(th)
    rx[i, :] = xr[i] + hk * H_TOTAL * nx
    ry[i, :] = yr[i] + hk * H_TOTAL * ny

# ── BL top curve (for triangulation boundary) ─────────────────────────────────
theta_bl = np.linspace(np.pi / 2, 0.0, 200)
bl_top = np.array([
    [a * np.cos(th) + H_TOTAL * outward_normal(th)[0],
     b * np.sin(th) + H_TOTAL * outward_normal(th)[1]]
    for th in theta_bl
])
bl_top_x = bl_top[:, 0]
bl_top_y = bl_top[:, 1]

def above_bl(cx, cy):
    """True if (cx, cy) is above the BL outer edge (in the far-field fluid)."""
    if cx < bl_top_x[0]:
        return cy > b + H_TOTAL + 0.02
    if cx > bl_top_x[-1]:
        return cy > 0.02
    return cy > np.interp(cx, bl_top_x, bl_top_y) + 0.02

# ── Figure ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 7))

# ── Solid body (inside ellipse) ───────────────────────────────────────────────
theta_fill = np.linspace(0.0, np.pi / 2, 400)
ax.fill_between(
    a * np.cos(theta_fill), 0.0, b * np.sin(theta_fill),
    color='#636e72', alpha=0.40, zorder=1, label='Solid body')
ax.text(0.65, 0.32, 'solid\nbody', fontsize=10, color='white',
        ha='center', style='italic', zorder=3, fontweight='bold')

# ── Outer fluid: Delaunay triangles above BL ──────────────────────────────────
y_far = 1.78
x_max = 2.55
rng = np.random.default_rng(7)
pts = []
for _ in range(220):
    cx = rng.uniform(0.0, x_max)
    lower = np.interp(cx, bl_top_x, bl_top_y) if cx <= bl_top_x[-1] else 0.0
    cy = rng.uniform(lower + 0.06, y_far)
    if cy < y_far:
        pts.append([cx, cy])

# Boundary constraints
for xi in np.linspace(0.0, x_max, 18):
    pts.append([xi, y_far])
for yi in np.linspace(b + H_TOTAL * 0.5, y_far, 8):
    pts.append([0.0, yi])
for i in range(0, len(bl_top), 4):   # BL top edge points
    pts.append(bl_top[i].tolist())
# Axis region to the right of nose
for xi in np.linspace(a + H_TOTAL * 0.1, x_max, 10):
    pts.append([xi, 0.02])

pts = np.array(pts)
tri = Delaunay(pts)
for simplex in tri.simplices:
    verts = pts[simplex]
    cx, cy = verts.mean(axis=0)
    if (0 <= cx <= x_max + 0.1 and 0 < cy <= y_far
            and above_bl(cx, cy)
            and cy > ellipse_y(cx) + 0.01):
        ax.add_patch(mpatches.Polygon(
            verts, closed=True,
            fc='#dfe6e9', ec='#b2bec3', lw=0.45, alpha=0.80, zorder=2))

ax.text(0.45, 1.42, 'fluid (triangles)', fontsize=9,
        color='#2980b9', ha='center', style='italic', zorder=4, fontweight='bold')

# ── Regular BL cells (left of BC zone, lighter) ───────────────────────────────
for i in range(N_REG):
    for k in range(N_LAY):
        xs_q = [rx[i, k], rx[i+1, k], rx[i+1, k+1], rx[i, k+1]]
        ys_q = [ry[i, k], ry[i+1, k], ry[i+1, k+1], ry[i, k+1]]
        frac = k / max(N_LAY - 1, 1)
        grey = 0.65 + 0.25 * frac
        ax.add_patch(mpatches.Polygon(
            list(zip(xs_q, ys_q)), closed=True,
            fc=(min(grey, 1.0), min(grey + 0.05, 1.0), min(grey + 0.15, 1.0)),
            ec='#74b9ff', lw=0.7, alpha=0.75, zorder=3))

# ── BC fan cells (orange-to-yellow, by BL depth) ──────────────────────────────
for i in range(N_COL):
    for k in range(N_LAY):
        xs_q = [gx[i, k], gx[i+1, k], gx[i+1, k+1], gx[i, k+1]]
        ys_q = [gy[i, k], gy[i+1, k], gy[i+1, k+1], gy[i, k+1]]
        frac = k / max(N_LAY - 1, 1)    # 0 = near wall, 1 = near BL top
        r_c = 0.96 - 0.28 * frac
        g_c = 0.52 + 0.38 * frac
        b_c = 0.10 + 0.60 * frac
        ax.add_patch(mpatches.Polygon(
            list(zip(xs_q, ys_q)), closed=True,
            fc=(r_c, g_c, b_c), ec='#e17055', lw=0.9, alpha=0.90, zorder=4))

# ── Scatter: BC quad nodes ────────────────────────────────────────────────────
ax.scatter(gx.ravel(), gy.ravel(), c='#d63031', s=28, zorder=8,
           label='BC quad nodes', marker='o')
ax.scatter(rx.ravel(), ry.ravel(), c='#2980b9', s=15, zorder=7,
           label='Regular BL nodes', marker='o', alpha=0.7)

# ── Profile curves ────────────────────────────────────────────────────────────
theta_all = np.linspace(np.pi / 2, 0.0, 600)
xp = a * np.cos(theta_all)
yp = b * np.sin(theta_all)

# Non-BC portion
ax.plot(xp[theta_all >= theta_start - 0.001],
        yp[theta_all >= theta_start - 0.001],
        '-', color='#2d3436', lw=3, zorder=6, label='Wall profile (ellipse)')
# BC zone
ax.plot(xp[theta_all <= theta_start + 0.001],
        yp[theta_all <= theta_start + 0.001],
        '-', color='#e17055', lw=3.5, zorder=6, label='Wall (BC zone)')

# ── Axis y=0 and y-axis x=0 ───────────────────────────────────────────────────
ax.axhline(0.0, color='#2d3436', lw=2.2, ls='-.', zorder=5,
           label='Symmetry axis  y = 0')
ax.axvline(0.0, color='#636e72', lw=1.2, ls=':', zorder=4, alpha=0.55)

# ── Key points ────────────────────────────────────────────────────────────────
ax.plot(a, 0.0, 'r^', ms=15, zorder=9, label=f'AxisPoint  ({a:.1f}, 0)')
ax.plot(x_sp, y_sp, 'bs', ms=13, zorder=9,
        label=f'StartPoint  ({x_sp:.2f}, {y_sp:.2f})')

# ── Annotations ───────────────────────────────────────────────────────────────
ax.annotate('AxisPoint\n(nose — wall ⊥ axis)',
            xy=(a, 0.0), xytext=(2.42, 0.52),
            arrowprops=dict(arrowstyle='->', color='#d63031', lw=2.2),
            color='#d63031', fontsize=10, ha='center', fontweight='bold', zorder=10)

ax.annotate('StartPoint',
            xy=(x_sp, y_sp), xytext=(x_sp - 0.72, y_sp + 0.30),
            arrowprops=dict(arrowstyle='->', color='#2980b9', lw=1.6),
            color='#2980b9', fontsize=10, zorder=10)

# 90° arc at nose (between wall tangent and axis)
arc_theta = np.linspace(np.pi / 2, np.pi, 40)
r_arc = 0.14
ax.plot(a + r_arc * np.cos(arc_theta), r_arc * np.sin(arc_theta),
        '-', color='#d63031', lw=1.7, zorder=8)
ax.text(a - 0.24, 0.18, '90°', fontsize=8, color='#d63031', zorder=9)

# Column count arrow
ax.annotate('', xy=(gx[-1, 0], -0.17), xytext=(gx[0, 0], -0.17),
            arrowprops=dict(arrowstyle='<->', color='#6c5ce7', lw=1.6), zorder=5)
ax.text((gx[0, 0] + gx[-1, 0]) / 2, -0.25,
        f'nbCornerColumns = {N_COL}', ha='center', color='#6c5ce7', fontsize=9)

# Layer count arrow on first column (StartPoint side)
top_pt = (gx[0, -1], gy[0, -1])
bot_pt = (gx[0, 0], gy[0, 0])
mid_x = top_pt[0] + 0.10
ax.annotate('', xy=(mid_x, bot_pt[1]), xytext=(mid_x, top_pt[1]),
            arrowprops=dict(arrowstyle='<->', color='navy', lw=1.5), zorder=5)
ax.text(mid_x + 0.06, (bot_pt[1] + top_pt[1]) / 2,
        f'nbLayers={N_LAY}\nratio={RATIO}\nhTotal={H_TOTAL}',
        color='navy', fontsize=8, va='center', zorder=5)

# Outward normal arrow at StartPoint
nx_sp, ny_sp = outward_normal(theta_start)
ax.annotate('', xy=(x_sp + nx_sp * 0.34, y_sp + ny_sp * 0.34),
            xytext=(x_sp, y_sp),
            arrowprops=dict(arrowstyle='->', color='#00b894', lw=2.0), zorder=9)
ax.text(x_sp + nx_sp * 0.44, y_sp + ny_sp * 0.44 + 0.06,
        'outward\nnormal', fontsize=8, color='#00b894', ha='center', zorder=9)

# Fan zone label
ax.text(1.82, 0.62, 'BC fan\n(quads)', fontsize=9,
        color='#c0392b', ha='center', style='italic', zorder=5, fontweight='bold')

# Regular BL label
ax.text(0.28, b + H_TOTAL * 0.6, 'regular\nBL', fontsize=8,
        color='#2980b9', ha='center', zorder=5)

# y-axis label
ax.text(-0.06, b / 2, 'y-axis\n(cut)', ha='right', fontsize=8,
        color='#636e72', va='center')

# axis label to the right
ax.text(2.65, 0.05, 'axis\n(y=0)', ha='center', fontsize=8,
        color='#2d3436', va='bottom')

# ── Formatting ────────────────────────────────────────────────────────────────
ax.set_xlabel('x   →   axial direction', fontsize=11)
ax.set_ylabel('y   →   radial direction', fontsize=11)
ax.set_title(
    'BoundaryCorner — half-ellipse profile cut by the y-axis\n'
    'Structured quad fan OUTSIDE the body (external axisymmetric CFD)',
    fontsize=12, fontweight='bold')
ax.legend(loc='upper right', fontsize=9, framealpha=0.95)
ax.set_aspect('equal')
ax.set_xlim(-0.3, 2.85)
ax.set_ylim(-0.38, 1.88)
ax.grid(alpha=0.3)

plt.tight_layout()
out_path = 'local/boundary_corner_sketch.png'
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"Saved: {out_path}")
plt.show()
