#!/usr/bin/env python3
"""Worked examples for the paper: writes ex4d.tex and ex5d.tex (table bodies)."""
import math
NA = 300
a = [0] * (NA + 1); a[0] = 1
for _ in range(6):
    for n in range(1, NA + 1):
        for m in range(n, NA + 1):
            a[m] += a[m - n]
def A(n): return a[n] if n >= 0 else 0
def cabs(D):
    if D < -1 or D % 4 in (1, 2): return 0
    return sum(A((D + 1 - j * j) // 4) * (1 if j == 0 else 2) for j in range(0, math.isqrt(D + 1) + 1) if (D + 1 - j * j) % 4 == 0)
def c(D): return (1 if D % 2 else -1) * cabs(D)

rows = []
for (l1, n, m, p) in [(1,1,1,0),(1,1,1,1),(1,2,1,1),(1,2,2,0),(1,2,2,1),(1,2,2,2),(1,3,3,3),(1,2,4,2),(1,4,4,4),(1,4,4,0),(1,6,6,6),(2,1,1,1),(3,1,1,1),(2,1,3,1)]:
    Q2h, P2h, QP = l1 * l1 * n, m, l1 * p
    l2 = math.gcd(math.gcd(Q2h, P2h), QP)
    assert math.gcd(l1, l2) == 1
    D = 4 * Q2h * P2h - QP * QP
    N = l1 * l2
    terms = [(s, (-1) ** (QP + 1) * s * c(D // (s * s))) for s in range(1, N + 1) if N % s == 0]
    tot = sum(t for _, t in terms)
    ts = " ".join(("$%+d$" % t) if i else ("$%d$" % t) for i, (_, t) in enumerate(terms))
    rows.append(f"{Q2h} & {P2h} & {QP} & {l1} & {l2} & {D} & {ts} & {tot}\\\\")
open("ex4d.tex", "w").write("\n".join(rows) + "\n")

def F(M, g, l):
    E = 4 * M - l * l
    return sum(s * c(E // (s * s)) for s in range(1, g + 1) if g % s == 0 and l % s == 0 and E % (s * s) == 0)
def eps(M, g, l): return (1 if abs(l) % 2 else -1) * F(M, g, abs(l))
def dt(M, g, J): return sum(b * eps(M, g, J - r) for r, b in zip(range(-2, 3), (1, 4, 6, 4, 1)))
rows = []
for (Q1, Q5, n) in [(1,1,1),(2,1,1),(1,1,3),(2,1,2),(3,1,3),(2,1,4),(4,1,4)]:
    M = Q1 * Q5 * n; g = math.gcd(Q1 * Q5, n)
    for J in range(0, math.isqrt(4 * M - 1) + 1):
        if 4 * M - J * J <= 0: continue
        x, y = dt(M, g, J), dt(M, g, J + 2)
        rows.append(f"$({Q1},{Q5},{n})$ & {M} & {g} & {J} & {4*M-J*J} & {eps(M,g,J)} & {x} & {x-y}\\\\")
for J in (29, 30, 31):
    M, g = 242, 11
    x, y = dt(M, g, J), dt(M, g, J + 2)
    rows.append(f"$(11,1,22)$ & {M} & {g} & {J} & {4*M-J*J} & {eps(M,g,J)} & {x} & {x-y}\\\\")
open("ex5d.tex", "w").write("\n".join(rows) + "\n")
# polar example
M, g = 242, 11
print("polar example eps:", {l: eps(M, g, l) for l in range(29, 38)})
print("dt, d at J=31:", dt(M, g, 31), dt(M, g, 31) - dt(M, g, 33), " J=29:", dt(M,g,29), dt(M,g,29)-dt(M,g,31))
# even s0 negative example
M, g = 8, 2
print("M=8 eps:", {l: eps(M, g, l) for l in range(0, 8)}, "dt(4)=", dt(M, g, 4), "d(4)=", dt(M,g,4)-dt(M,g,6))
