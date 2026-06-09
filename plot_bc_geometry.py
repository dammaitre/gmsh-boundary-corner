"""
Visualisation de la géométrie BoundaryCorner.
Scatter + grille de quads sur les deux nappes (inférieure et supérieure).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# ─── Paramètres ──────────────────────────────────────────────────────────────
N_col = 5        # NbCornerColumns
N_lay = 5        # NbLayers
h1    = 0.012
ratio = 1.30
omega = 1.0
x_S   = 0.68     # StartPoint.x
x_A   = 1.00     # AxisPoint.x  (coin sur l'axe y=0)

# ─── Profil courbe (descend de y>0 en x_S vers y=0 en x_A) ─────────────────
def y_prof(x):
    t = np.clip((x - x_S) / (x_A - x_S), 0.0, 1.0)
    return 0.15 * (1.0 - t) ** 0.7

def dy_dx(x, dx=1e-6):
    return (y_prof(x + dx) - y_prof(x - dx)) / (2.0 * dx)

def unit_normal(x):
    """Normale unitaire au profil, pointant vers y+ (intérieur domaine)."""
    tx, ty = 1.0, float(dy_dx(x))
    L = np.hypot(tx, ty)
    nx, ny = -ty / L, tx / L
    return (nx, ny) if ny >= 0.0 else (-nx, -ny)

# ─── Hauteur totale BL ───────────────────────────────────────────────────────
hTotal = h1 * (ratio ** N_lay - 1.0) / (ratio - 1.0)

# ─── Colonnes : trouver eps + l0 par point fixe ──────────────────────────────
# Système : l0 * eps^(N-1) = delta1 = h1
#           sum_{i=0}^{N-1} l0*eps^i = x_A - x_S
# => l0 = (x_A-x_S)*(eps-1)/(eps^N - 1)  et  l0 = h1/eps^(N-1)
# Résolution numérique simple par dichotomie sur eps
L_span = x_A - x_S
delta1 = h1

def residual(eps):
    if abs(eps - 1.0) < 1e-12:
        l0_a = L_span / N_col
        l0_b = delta1
    else:
        l0_a = L_span * (eps - 1.0) / (eps**N_col - 1.0)
        l0_b = delta1 / eps**(N_col - 1)
    return l0_a - l0_b

# dichotomie : eps ∈ (0.3, 0.999) pour convergence vers le coin
lo, hi = 0.30, 0.999
for _ in range(60):
    mid = 0.5 * (lo + hi)
    (lo if residual(mid) < 0 else hi)
    if residual(mid) < 0:
        lo = mid
    else:
        hi = mid
eps = 0.5 * (lo + hi)
l0  = L_span * (eps - 1.0) / (eps**N_col - 1.0)

# positions x des N_col+1 colonnes
x_cols = [x_S]
for i in range(N_col):
    x_cols.append(x_cols[-1] + l0 * eps**i)
x_cols[-1] = x_A   # forcer exactement sur le coin

# ─── 3 rangées de contrôle ───────────────────────────────────────────────────
ptsProfil = np.array([[x, float(y_prof(x))] for x in x_cols])
ptsAxe    = np.array([[x, 0.0]              for x in x_cols])
ptsBLTop  = np.array([[x + unit_normal(x)[0]*hTotal,
                        float(y_prof(x)) + unit_normal(x)[1]*hTotal]
                       for x in x_cols])

# ─── Grilles de nœuds (nappe inf et sup) ────────────────────────────────────
def hk(k):
    return 0.0 if k == 0 else h1 * omega * (ratio**k - 1.0) / (ratio - 1.0)

g_inf = np.zeros((N_col + 1, N_lay + 1, 2))  # profil (k=0) → axe (k=K)
g_sup = np.zeros((N_col + 1, N_lay + 1, 2))  # profil (k=0) → BLTop (k=K)

for i in range(N_col + 1):
    xp, yp = ptsProfil[i]
    xa, ya = ptsAxe[i]
    xb, yb = ptsBLTop[i]
    for k in range(N_lay + 1):
        f = hk(k) / hTotal
        g_inf[i, k] = [xp + (xa - xp)*f, yp + (ya - yp)*f]
        g_sup[i, k] = [xp + (xb - xp)*f, yp + (yb - yp)*f]

# ─── Figure ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 6))
ax.set_aspect("equal")
ax.set_facecolor("#f5f7fa")

C_INF = "#2b7fb8"   # bleu  – nappe inférieure
C_SUP = "#d95f02"   # orange – nappe supérieure

# ── Lignes de grille ─────────────────────────────────────────────────────────
def grid_lines(g, color, lw=0.9, alpha=0.5):
    segs = []
    nc, nl, _ = g.shape
    for i in range(nc):
        segs.append(g[i, :, :])
    for k in range(nl):
        segs.append(g[:, k, :])
    ax.add_collection(LineCollection(segs, colors=color, linewidths=lw,
                                     alpha=alpha, zorder=2))

grid_lines(g_inf, C_INF)
grid_lines(g_sup, C_SUP)

# ── Nœuds internes ───────────────────────────────────────────────────────────
def scatter_nodes(g, color, exclude_k0=True):
    pts = g[:, (1 if exclude_k0 else 0):, :].reshape(-1, 2)
    ax.scatter(pts[:, 0], pts[:, 1], s=14, c=color, zorder=4,
               alpha=0.75, linewidths=0)

scatter_nodes(g_inf, C_INF)
scatter_nodes(g_sup, C_SUP)

# ── 3 rangées de contrôle ─────────────────────────────────────────────────────
ax.scatter(*ptsAxe.T,    s=55, c="#1a5276",  zorder=6, marker="o",
           label="ptsAxe  (y = 0, axe de symétrie)")
ax.scatter(*ptsProfil.T, s=70, c="#c0392b",  zorder=6, marker="D",
           label="ptsProfil  (sur le profil)")
ax.scatter(*ptsBLTop.T,  s=55, c="#1e8449",  zorder=6, marker="s",
           label="ptsBLTop  (offset normal × hTotal)")

# ── Profil complet (contexte) ─────────────────────────────────────────────────
x_ctx = np.linspace(0.3, x_A, 400)
y_ctx = np.where(x_ctx >= x_S, y_prof(x_ctx),
                 float(y_prof(x_S)) * (x_ctx - 0.3) / (x_S - 0.3))
ax.plot(x_ctx, y_ctx, c="#888", lw=1.4, ls="--", alpha=0.55,
        zorder=1, label="profil (contexte)")

# ── Axe de symétrie ───────────────────────────────────────────────────────────
ax.axhline(0, color="#444", lw=1.0, ls="-", alpha=0.35, zorder=1)
ax.text(0.32, -0.008, "axe de symétrie  (y = 0)", fontsize=7.5,
        color="#444", alpha=0.7)

# ── Annotations ───────────────────────────────────────────────────────────────
ax.annotate("AxisPoint\n(coin)", xy=(x_A, 0),
            xytext=(x_A + 0.04, -0.025), fontsize=8.5, color="#1a5276",
            arrowprops=dict(arrowstyle="-|>", color="#1a5276", lw=0.9))

ax.annotate("StartPoint", xy=(x_S, float(y_prof(x_S))),
            xytext=(x_S - 0.07, float(y_prof(x_S)) + 0.03),
            fontsize=8.5, color="#c0392b",
            arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=0.9))

# flèche double hTotal
x_arr = x_A + 0.03
ax.annotate("", xy=(x_arr, 0), xytext=(x_arr, hTotal),
            arrowprops=dict(arrowstyle="<->", color="#1e8449", lw=1.2))
ax.text(x_arr + 0.012, hTotal * 0.5, f"hTotal\n= {hTotal:.3f}",
        fontsize=7.5, color="#1e8449", va="center")

# labels nappes
x_mid = 0.5 * (x_S + x_A)
ax.text(x_mid, -0.022, "nappe inférieure\n(profil → axe)",
        ha="center", fontsize=7.5, color=C_INF, style="italic")
ax.text(x_mid, hTotal + 0.018, "nappe supérieure\n(profil → offset BL)",
        ha="center", fontsize=7.5, color=C_SUP, style="italic")

# progression eps entre colonnes
ax.annotate("", xy=(x_cols[1], -0.038), xytext=(x_cols[0], -0.038),
            arrowprops=dict(arrowstyle="<->", color="#555", lw=0.9))
ax.annotate("", xy=(x_cols[2], -0.038), xytext=(x_cols[1], -0.038),
            arrowprops=dict(arrowstyle="<->", color="#555", lw=0.9))
ax.text(0.5*(x_cols[0]+x_cols[1]), -0.048, "l₀", ha="center",
        fontsize=8, color="#555")
ax.text(0.5*(x_cols[1]+x_cols[2]), -0.048, "l₀·ε", ha="center",
        fontsize=8, color="#555")
ax.text(0.5*(x_cols[-2]+x_cols[-1]) - 0.01, -0.048,
        f"… δ₁=h₁\n(ε≈{eps:.2f})", ha="center", fontsize=7, color="#555")

# progression h1, ratio entre rangées (sur la première colonne)
for k in range(1, N_lay + 1):
    xa_k  = g_inf[0, k - 1, 0] - 0.022
    ya_k  = g_inf[0, k - 1, 1]
    xa_k1 = g_inf[0, k,     0] - 0.022
    ya_k1 = g_inf[0, k,     1]
    ax.annotate("", xy=(xa_k1, ya_k1), xytext=(xa_k, ya_k),
                arrowprops=dict(arrowstyle="<->", color=C_INF, lw=0.7))

ax.text(g_inf[0, 0, 0] - 0.04, 0.5*(g_inf[0,0,1]+g_inf[0,1,1]),
        "h₁", va="center", fontsize=7.5, color=C_INF)

# ── Légende & finalisation ────────────────────────────────────────────────────
ax.legend(loc="upper left", fontsize=8.5, framealpha=0.92,
          edgecolor="#ccc", borderpad=0.8)
ax.set_title(
    f"BoundaryCorner — géométrie structurée  "
    f"(N_col={N_col}, N_lay={N_lay}, h₁={h1}, ratio={ratio})",
    fontsize=10.5, pad=10)
ax.set_xlabel("x", fontsize=9)
ax.set_ylabel("y", fontsize=9)
ax.set_xlim(0.28, 1.17)
ax.set_ylim(-0.065, 0.24)
ax.tick_params(labelsize=8)
ax.grid(True, ls=":", alpha=0.3, zorder=0)

plt.tight_layout()
out = "boundary_corner_geometry.png"
plt.savefig(out, dpi=160, bbox_inches="tight")
print(f"Sauvegardé : {out}")
