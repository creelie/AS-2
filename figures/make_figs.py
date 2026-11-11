#!/usr/bin/env python3
"""Figures for the paper. Requires numpy and matplotlib."""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

plt.rcParams.update({"font.family": "serif", "font.size": 10, "mathtext.fontset": "cm"})
NA = 200
a = [0] * (NA + 1); a[0] = 1
for _ in range(6):
    for n in range(1, NA + 1):
        for m in range(n, NA + 1):
            a[m] += a[m - n]

def cabs(D):
    if D < -1 or D % 4 in (1, 2): return 0
    s, j = 0, 0
    while j * j <= D + 1:
        if (D + 1 - j * j) % 4 == 0:
            s += a[(D + 1 - j * j) // 4] * (1 if j == 0 else 2)
        j += 1
    return s
def c(D): return (1 if D % 2 else -1) * cabs(D)
def W(E):
    if E < 1: return 0
    return sum(s * s * cabs(E // (s * s)) for s in range(2, math.isqrt(E) + 1) if E % (s * s) == 0)

# Figure 1: the coefficients c(4k - l^2) as a 3D stem plot, coloured by sign
fig = plt.figure(figsize=(6.4, 4.6))
ax = fig.add_subplot(111, projection="3d")
for k in range(0, 9):
    for l in range(-6, 7):
        D = 4 * k - l * l
        if D < -1: continue
        v = c(D)
        h = math.log10(abs(v)) if v else 0
        col = "#b2182b" if v > 0 else "#2166ac"
        ax.plot([k, k], [l, l], [0, h], color=col, lw=2.2)
        ax.scatter([k], [l], [h], color=col, s=14, depthshade=False)
ax.set_xlabel("$k$"); ax.set_ylabel("$l$"); ax.set_zlabel(r"$\log_{10}|\hat c(4k-l^2)|$", labelpad=6)
ax.view_init(elev=24, azim=-58)
ax.text2D(0.02, 0.93, "red: $\\hat c>0$ ($D$ odd)   blue: $\\hat c<0$ ($D$ even)", transform=ax.transAxes, fontsize=9)
fig.tight_layout()
fig.savefig("fig_coefficients.png", dpi=220); plt.close(fig)

# Figure 2: the margin 2W*(D)/G(D) and (W* + (D+4)(D+8)/16)/G(D) for 3 <= D <= 359
Ds, ra, rb = [], [], []
Wv = [0] + [W(E) for E in range(1, 400)]
for D in range(3, 360):
    if D % 4 in (1, 2): continue
    G = cabs(D) - cabs(D - 4)
    Ws = max(Wv[1:D + 1])
    Ds.append(D); ra.append(2 * Ws / G); rb.append((Ws + (D + 4) * (D + 8) / 16) / G)
fig, ax = plt.subplots(figsize=(6.2, 3.4))
ax.semilogy(Ds, rb, "o", ms=2.6, color="#5e3c99", label=r"$(W^*(D)+(D+4)(D+8)/16)/G(D)$")
ax.semilogy(Ds, ra, "s", ms=2.4, color="#e66101", label=r"$2W^*(D)/G(D)$")
ax.axhline(1, color="k", lw=0.8)
ax.annotate("$D=12$: $64/152$", xy=(12, ra[Ds.index(12)]), xytext=(60, 0.5),
            arrowprops=dict(arrowstyle="->", lw=0.7), fontsize=9)
ax.set_xlabel("$D$"); ax.set_ylabel("ratio"); ax.set_ylim(1e-12, 3)
ax.legend(loc="lower left", fontsize=8.5, frameon=False)
fig.tight_layout(); fig.savefig("fig_margins.png", dpi=220); plt.close(fig)

# Figure 3: the exponents of the analytic range
def Lk(k): return 2 * math.sqrt(6 * k) - 6 - math.log(k) + math.log(5 / (k + 5))
def S2(k): return math.log(4) + math.log(1.0001 * 4 * (math.sqrt(k + 1) + 1) / (k + 1) ** 1.5) + math.pi * math.sqrt(k + 1)
def S3(k): return math.log(4) + 1.5 * math.log(4 * k) + math.log(math.sqrt(4 * k / 9 + 1) + 1) + 2 * math.pi / 3 * math.sqrt(k + 2.25)
ks = np.arange(50, 241)
fig, ax = plt.subplots(figsize=(6.2, 3.3))
ax.plot(ks, [Lk(k) for k in ks], color="#1b7837", lw=1.6, label=r"$\Lambda(k)$ (lower bound for $\log G$)")
ax.plot(ks, [S3(k) for k in ks], color="#762a83", lw=1.4, ls="--", label=r"$\log 4+\log S_3(k)$")
ax.plot(ks, [S2(k) for k in ks], color="#e08214", lw=1.4, ls=":", label=r"$\log 4+\log S_2(k)$")
ax.axvline(90, color="0.4", lw=0.8); ax.text(92, 22, "$k=90$", fontsize=9)
ax.set_xlabel("$k$"); ax.legend(loc="upper left", fontsize=8.5, frameon=False)
fig.tight_layout(); fig.savefig("fig_exponents.png", dpi=220); plt.close(fig)

# Figure 4: the five-dimensional index at fixed J, log10 d(J), g = 1, as 3D bars
def F(M, g, l):
    E = 4 * M - l * l
    if E < -1: return 0
    return sum(s * c(E // (s * s)) for s in range(1, g + 1) if g % s == 0 and l % s == 0 and E % (s * s) == 0)
pref = {2: 1, 1: -4, 0: 6, -1: -4, -2: 1}
def dt(M, g, J): return (1 if J % 2 else -1) * sum(v * F(M, g, J - r) for r, v in pref.items())
fig = plt.figure(figsize=(6.4, 4.6))
ax = fig.add_subplot(111, projection="3d")
xs, ys, zs, cs = [], [], [], []
for M in range(1, 21):
    for J in range(0, 2 * math.isqrt(M) + 2):
        if 4 * M - J * J <= 0: continue
        val = dt(M, 1, J) - dt(M, 1, J + 2)
        xs.append(M); ys.append(J); zs.append(math.log10(val)); cs.append(J)
cmap = plt.get_cmap("viridis")
ax.bar3d(np.array(xs) - 0.35, np.array(ys) - 0.35, 0, 0.7, 0.7, zs,
         color=[cmap(j / max(cs)) for j in cs], shade=True, edgecolor="none")
ax.set_xlabel("$M=Q_1Q_5n$"); ax.set_ylabel("$J$"); ax.set_zlabel(r"$\log_{10} d(J)$", labelpad=4)
ax.view_init(elev=30, azim=48); ax.set_xticks([1, 5, 10, 15, 20]); ax.set_yticks(range(0, 9, 2))
fig.tight_layout(); fig.savefig("fig_fixedJ.png", dpi=220); plt.close(fig)
print("figures written")

# Figure 5: a(n) against its two bounds, drawn as three ribbons in depth
from fractions import Fraction
def bcoef(n):
    if n == 0: return 1
    return sum(Fraction(6 ** j, math.factorial(j)) * math.comb(n - 1, j - 1) for j in range(1, n + 1))
ns = np.arange(1, 161)
la = np.array([math.log(a[n]) for n in ns])
lu = np.array([math.log(1.0001) + 2 * math.pi * math.sqrt(n) - math.log(8 * n ** 1.5) for n in ns])
lb = np.array([math.log(float(bcoef(n))) for n in ns])
le = np.array([2 * math.sqrt(6 * n) - 6 - math.log(n) for n in ns])
fig = plt.figure(figsize=(6.6, 4.4))
ax = fig.add_subplot(111, projection="3d")
for depth, arr, col, lab in [(0, lu, "#d6604d", r"upper bound for $a(n)$"), (1, la, "#222222", r"$\log a(n)$"),
                             (2, lb, "#4393c3", r"$\log b(n)$"), (3, le, "#1b7837", r"$2\sqrt{6n}-6-\log n$")]:
    ax.plot(ns, [depth] * len(ns), arr, color=col, lw=1.6, label=lab)
    ax.fill_between(ns, np.full(len(ns), depth), arr, ns, np.full(len(ns), depth), np.zeros(len(ns)), color=col, alpha=0.10)
ax.set_xlabel("$n$"); ax.set_yticks([]); ax.set_zlabel("natural logarithm", labelpad=4)
ax.view_init(elev=18, azim=-70)
ax.legend(loc="upper left", fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig("fig_abounds.png", dpi=220); plt.close(fig)

# Figure 6: the four-dimensional index over (Q^2/2, P^2/2) at Q.P = 0 and l1 = 1, as 3D bars
fig = plt.figure(figsize=(6.4, 4.6))
ax = fig.add_subplot(111, projection="3d")
for n in range(1, 8):
    for m in range(1, 8):
        l2 = math.gcd(n, m); D = 4 * n * m
        v = -sum(s * c(D // (s * s)) for s in range(1, l2 + 1) if l2 % s == 0)
        h = math.log10(v)
        col = plt.cm.viridis(min(1, (l2 - 1) / 6))
        ax.bar3d(n - 0.35, m - 0.35, 0, 0.7, 0.7, h, color=col, shade=True, edgecolor="k", linewidth=0.2)
ax.set_xlabel(r"$Q^2/2$"); ax.set_ylabel(r"$P^2/2$"); ax.set_zlabel(r"$\log_{10} d(Q,P)$", labelpad=6)
ax.view_init(elev=26, azim=-130)
ax.text2D(0.02, 0.93, r"$Q\cdot P=0$, $\ell_1=1$; colour: $\ell_2=\gcd(Q^2/2,P^2/2)$", transform=ax.transAxes, fontsize=9)
fig.tight_layout(); fig.savefig("fig_index4d.png", dpi=220); plt.close(fig)

# Figure 7: the pairing of Theorem C on the parabola D(l) = 4M - l^2, with the polar point of M = 242
M, J = 242, 29
ls_ = np.linspace(20, 36, 400)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(ls_, 4 * M - ls_ ** 2, color="0.35", lw=1.2)
ax.axhline(0, color="k", lw=0.6); ax.axhline(-1, color="0.6", lw=0.5, ls=":")
pts = {J - 2: "#1b7837", J - 1: "#1b7837", J: "#1b7837", J + 2: "#b2182b", J + 3: "#b2182b", J + 4: "#b2182b"}
for l, col in pts.items():
    ax.plot([l], [4 * M - l * l], "o", color=col, ms=6)
for r, (lo, hi) in enumerate([(J, J + 2), (J - 1, J + 3), (J - 2, J + 4)]):
    y = 260 + 70 * r
    ax.annotate("", xy=(hi, y), xytext=(lo, y), arrowprops=dict(arrowstyle="<->", lw=0.9, color="#5e3c99"))
    ax.text((lo + hi) / 2, y + 12, f"$c_{r}={[5, 4, 1][r]}$", ha="center", fontsize=9, color="#5e3c99")
ax.plot([33], [4 * M - 33 * 33], "*", color="#e66101", ms=13)
ax.annotate(r"polar point $l=33$: $D=-121=-11^2$, $\varepsilon(33)=11$", xy=(33, -121), xytext=(22.5, -330), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.7))
ax.set_xlabel("$l$"); ax.set_ylabel(r"$D(l)=4M-l^2$"); ax.set_xlim(20, 36); ax.set_ylim(-420, 520)
ax.set_title(r"$M=242$, $g=11$, $J=29$: brackets $\varepsilon(|J-r|)-\varepsilon(J+2+r)$", fontsize=10)
fig.tight_layout(); fig.savefig("fig_pairing.png", dpi=220); plt.close(fig)
