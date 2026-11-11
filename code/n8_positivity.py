#!/usr/bin/env python3
"""Exact checks for the positivity of the N=8 black hole indices.

Standard library only. Checks:
 (1) the theta formula (Lemma 3.1) against a direct expansion of -theta_1^2/eta^6,
     and against the product formula (7) of the paper;
 (2) the bounds of Section 4 on the coefficients a(n) of prod(1-q^n)^-6;
 (3) Proposition 5.1 for 3 <= D <= 4*K0 - 1 (K0 = 90) in exact arithmetic;
 (4) Theorems A-C by brute force on all charges in a box, and on every 5D charge
     with M <= 600 in which a polar term (Lemma 7.1(3)) enters, with Lemma 7.2;
 (5) the same for the U-duality invariant ansatz of Maldacena-Moore-Strominger.
Prints a summary and exits with status 1 on any failure.
"""
import sys, math
from fractions import Fraction

K0 = 90
DMAX = 4 * K0 + 40
MPOL = 600          # range of M = Q1 Q5 n for the polar-term test
NA = max(DMAX, 4 * MPOL + 8) // 4 + 3
fails = []

def check(cond, msg):
    if not cond:
        fails.append(msg)

# a(n): coefficients of prod_{n>=1} (1-q^n)^{-6}
a = [0] * (NA + 1); a[0] = 1
for _ in range(6):
    for n in range(1, NA + 1):
        for m in range(n, NA + 1):
            a[m] += a[m - n]

def A(n):
    return a[n] if n >= 0 else 0

def cabs(D):
    """|c(D)| from the theta formula (Lemma 3.1)."""
    if D < -1 or D % 4 in (1, 2):
        return 0
    s, j = 0, 0
    while j * j <= D + 1:
        if (D + 1 - j * j) % 4 == 0:
            s += A((D + 1 - j * j) // 4) * (1 if j == 0 else 2)
        j += 1
    return s

def c(D):
    return (1 if D % 2 else -1) * cabs(D)

# (1) direct expansion of -theta_1(tau,z)^2 / eta(tau)^6 = sum c(4k-l^2) q^k y^l
KD = 30
# theta_1^2 = -sum_{r,s in Z+1/2} (-1)^{r+s-1} q^{(r^2+s^2)/2} y^{r+s}; work with 2r, 2s odd
th = {}
R = int(math.isqrt(8 * KD)) + 3
for r2 in range(-R, R + 1, 2):
    pass
for r2 in range(-2 * R - 1, 2 * R + 2):
    if r2 % 2 == 0: continue
    for s2 in range(-2 * R - 1, 2 * R + 2):
        if s2 % 2 == 0: continue
        e = (r2 * r2 + s2 * s2) // 8  # (r^2+s^2)/2 - 1/4 ; q-power shifted by 1/4
        if e > KD: continue
        sign = (-1) ** ((r2 + s2) // 2 - 1)
        l = (r2 + s2) // 2
        th[(e, l)] = th.get((e, l), 0) - sign   # theta_1^2 * q^{-1/4}
# multiply by q^{1/4}/eta^6 * q^{-1/4} = prod (1-q^n)^-6, and by -1
ser = {}
for (e, l), v in th.items():
    for n in range(0, KD - e + 1):
        ser[(e + n, l)] = ser.get((e + n, l), 0) + v * a[n]
for (k, l), v in ser.items():
    if k <= KD:
        check(-v == c(4 * k - l * l), f"theta formula at k={k}, l={l}")
# (1') the product formula (zeta - 2 + 1/zeta) prod (1-q^n zeta)^2 (1-q^n/zeta)^2 (1-q^n)^-4, up to q^12
KP = 12
def pmul(A_, B_):
    C_ = {}
    for (i, j), u in A_.items():
        for (k, l), v in B_.items():
            if i + k <= KP: C_[(i + k, j + l)] = C_.get((i + k, j + l), 0) + u * v
    return C_
prod = {(0, 1): 1, (0, 0): -2, (0, -1): 1}
for n in range(1, KP + 1):
    for _ in range(2):
        prod = pmul(prod, {(0, 0): 1, (n, 1): -1}); prod = pmul(prod, {(0, 0): 1, (n, -1): -1})
    inv = {(0, 0): 1}
    for _ in range(4):
        inv = pmul(inv, {(n * e, 0): 1 for e in range(0, KP // n + 1)})
    prod = pmul(prod, inv)
for (k, l), v in prod.items():
    check(v == c(4 * k - l * l), f"product formula at k={k}, l={l}")
check(c(-1) == 1 and c(0) == -2 and c(3) == 8 and c(4) == -12, "first values")

# (2) bounds of Section 4
for k in range(1, NA + 1):
    b = sum(Fraction(6 ** j, math.factorial(j)) * math.comb(k - 1, j - 1) for j in range(1, k + 1))
    check(a[k] >= b, f"a(k) >= b(k) at {k}")
    check(5 * a[k] <= (k + 5) * (a[k] - a[k - 1]), f"difference bound at {k}")
    check(a[k] <= math.exp(2 * math.pi * math.sqrt(k)), f"crude upper bound at {k}")
    check(a[k] <= 1.0001 * math.exp(2 * math.pi * math.sqrt(k)) / (8 * k ** 1.5), f"upper bound at {k}")
    if k >= 50:
        check(b >= Fraction(math.exp(2 * math.sqrt(6 * k) - 6)) / k, f"b(k) lower bound at {k}")

# (3) key proposition, exact, 3 <= D <= 4 K0 - 1
def W(E):
    if E < 1: return 0
    return sum(s * s * cabs(E // (s * s)) for s in range(2, math.isqrt(E) + 1) if E % (s * s) == 0)
Wv = [0] + [W(E) for E in range(1, DMAX + 1)]
worst_a, worst_b = (0, None), (0, None)
for D in range(3, 4 * K0):
    if D % 4 in (1, 2):
        continue
    G = cabs(D) - cabs(D - 4)
    k = (D + 1) // 4
    check(G >= (2 if D % 4 == 0 else 1) * (a[k] - a[k - 1]), f"monotonicity at D={D}")
    Ws = max(Wv[1:D + 1])
    check(G > 2 * Ws, f"key proposition (a) at D={D}")
    check(16 * (G - Ws) > (D + 4) * (D + 8), f"key proposition (b) at D={D}")
    check(cabs(D) <= (math.isqrt(D + 1) + 1) * a[k], f"upper bound for |c| at D={D}")
    ra, rb = Fraction(2 * Ws, G), Fraction(16 * Ws + (D + 4) * (D + 8), 16 * G)
    if ra > worst_a[0]: worst_a = (ra, D)
    if rb > worst_b[0]: worst_b = (rb, D)
# analytic range: h_i(K0) > 0 (the derivatives are positive for k >= 9, see the paper)
def Lk(k): return 2 * math.sqrt(6 * k) - 6 - math.log(k) + math.log(5 / (k + 5))
def S2(k): return math.log(1.0001 * 4 * (math.sqrt(k + 1) + 1) / (k + 1) ** 1.5) + math.pi * math.sqrt(k + 1)
def S3(k): return 1.5 * math.log(4 * k) + math.log(math.sqrt(4 * k / 9 + 1) + 1) + 2 * math.pi / 3 * math.sqrt(k + 2.25)
def S4(k): return math.log(2 * (k + 1) * (k + 2))
for f in (S2, S3, S4):
    check(Lk(K0) - math.log(4) - f(K0) > 0.5, f"analytic inequality {f.__name__} at K0")

# (4) brute force on charges
def d4(n, m, p, l1):
    """Sen's formula: Q = l1 Q0, P = P0, Q0^2=2n, P0^2=2m, Q0.P0=p, Q0 ^ P0 primitive."""
    Q2, P2, QP = 2 * n * l1 * l1, 2 * m, p * l1
    l2 = math.gcd(math.gcd(Q2 // 2, P2 // 2), QP)
    if math.gcd(l1, l2) != 1: return None
    Dl = Q2 * P2 - QP * QP
    if Dl <= 0: return None
    N = l1 * l2
    tot = sum(s * c(Dl // (s * s)) for s in range(1, N + 1) if N % s == 0)
    return (-1) ** (QP + 1) * tot
cnt4 = 0
for l1 in range(1, 7):
    for n in range(0, 9):
        for m in range(0, 9):
            for p in range(-8, 9):
                if 4 * n * m * l1 * l1 - p * p * l1 * l1 > DMAX: continue
                v = d4(n, m, p, l1)
                if v is None: continue
                cnt4 += 1
                check(v > 0, f"4D index at n={n}, m={m}, p={p}, l1={l1}: {v}")

def F(M, g, l):
    E = 4 * M - l * l
    # every s | gcd(g, l) is kept, including the polar terms E/s^2 = -1 with E <= -4
    return sum(s * c(E // (s * s)) for s in range(1, g + 1) if g % s == 0 and l % s == 0 and E % (s * s) == 0)
pref = {2: 1, 1: -4, 0: 6, -1: -4, -2: 1}
def dt(Q1, Q5, N, J):
    M = Q1 * Q5 * N; g = math.gcd(Q1 * Q5, N)
    return (1 if J % 2 else -1) * sum(v * F(M, g, J - r) for r, v in pref.items())
cnt5 = 0
for Q1 in range(1, 9):
    for Q5 in range(1, 9):
        if math.gcd(Q1, Q5) != 1: continue
        for N in range(1, 25):
            M = Q1 * Q5 * N
            if 4 * M > DMAX: continue
            for J in range(0, math.isqrt(4 * M) + 1):
                if 4 * M - J * J <= 0: continue
                x, y = dt(Q1, Q5, N, J), dt(Q1, Q5, N, J + 2)
                cnt5 += 1
                check(x > 0, f"5D index at {Q1},{Q5},{N},{J}")
                check(dt(Q1, Q5, N, -J) == x, f"5D symmetry at {Q1},{Q5},{N},{J}")
                check(x - y > 0, f"5D fixed-J index at {Q1},{Q5},{N},{J}")

# (4b) charges whose 5D index contains polar terms 4M - l^2 = -s^2 with s >= 2
cntp, polar = 0, 0
for M in range(1, MPOL + 1):
    for Q1Q5 in [d for d in range(1, M + 1) if M % d == 0]:
        N = M // Q1Q5; g = math.gcd(Q1Q5, N)
        if g == 1: continue          # no polar term with s >= 2 unless s | g
        # gcd(Q1, Q5) = 1 is possible for every Q1Q5 (take Q5 = 1)
        for J in range(0, math.isqrt(4 * M) + 1):
            if 4 * M - J * J <= 0: continue
            hit = False
            for t in range(1, 5):
                l = J + t; E = 4 * M - l * l
                s0 = math.isqrt(-E) if E < 0 else 0
                if E <= -4 and s0 * s0 == -E and g % s0 == 0 and l % s0 == 0:
                    hit = True
                    check(s0 <= (4 * M - J * J) + t * t, f"polar bound at M={M},J={J},t={t}")
            if not hit: continue
            polar += 1
            x, y = dt(Q1Q5, 1, N, J), dt(Q1Q5, 1, N, J + 2)
            check(x > 0 and x - y > 0, f"5D index with polar terms at M={M}, Q1Q5={Q1Q5}, J={J}: {x}, {x-y}")
check(dt(11, 1, 22, 31) - dt(11, 1, 22, 33) > 0 and F(242, 11, 33) == 11, "example M = 242")

# (5) the U-duality invariant ansatz of Maldacena-Moore-Strominger for gcd(Q1,Q5,N)=1
def wmms(Q1, Q5, N, s):
    vals = [N * Q1, N * Q5, Q1 * Q5]
    if any(v % s for v in vals) or (N * Q1 * Q5) % (s * s): return 0
    lst = [N, Q1, Q5, s, N * Q1 // s, N * Q5 // s, Q1 * Q5 // s, N * Q1 * Q5 // (s * s)]
    g = 0
    for v in lst: g = math.gcd(g, v)
    return s * sum(1 for dlt in range(1, g + 1) if g % dlt == 0)
def Fm(Q1, Q5, N, l):
    M = Q1 * Q5 * N; E = 4 * M - l * l
    top = math.isqrt(abs(E)) if E != 0 else math.isqrt(M)
    return sum(wmms(Q1, Q5, N, s) * c(E // (s * s)) for s in range(1, top + 1)
               if (l % s == 0) and E % (s * s) == 0)
def dtm(Q1, Q5, N, J):
    return (1 if J % 2 else -1) * sum(v * Fm(Q1, Q5, N, J - r) for r, v in pref.items())
cntm = 0
for Q1 in range(1, 7):
    for Q5 in range(1, 7):
        for N in range(1, 13):
            if math.gcd(math.gcd(Q1, Q5), N) != 1: continue
            M = Q1 * Q5 * N
            if 4 * M > DMAX: continue
            for J in range(0, math.isqrt(4 * M) + 1):
                if 4 * M - J * J <= 0: continue
                x, y = dtm(Q1, Q5, N, J), dtm(Q1, Q5, N, J + 2)
                cntm += 1
                check(x > 0 and x - y > 0, f"MMS ansatz at {Q1},{Q5},{N},{J}")

print(f"theta formula checked for k <= {KD}; product formula for k <= {KP}")
print(f"bounds checked for 1 <= k <= {NA}")
print(f"key proposition checked for 3 <= D <= {4*K0-1}; worst ratio (a) {float(worst_a[0]):.4f} at D = {worst_a[1]}, (b) {float(worst_b[0]):.4f} at D = {worst_b[1]}")
print(f"5D charges with polar terms checked (M <= {MPOL}): {polar}")
print(f"4D charges checked: {cnt4}; 5D (Q1,Q5,N,J) checked: {cnt5}; MMS ansatz checked: {cntm}")
if fails:
    print("FAILURES:", len(fails)); print("\n".join(fails[:20])); sys.exit(1)
print("all checks passed")
