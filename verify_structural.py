"""
verify_structural.py -- numerically confirm each load-bearing STEP of the
structural proof that, for integers A,B,C >= 1 with b_k = A k^2 + B k + C and the
standard continuant q_0=1, q_1=b_1, q_n = b_n q_{n-1} + q_{n-2},

    log q_n = n log A + 2 log(n!) + (B/A) log n + (K_Gamma + delta) + O(1/n),
    K_Gamma = -log( Gamma(1-r1) Gamma(1-r2) ),   delta = log R_inf,

with r1,r2 the roots of A x^2+B x+C (r1+r2=-B/A, r1 r2=C/A) and
R_inf = lim q_n / prod_{k=1}^n b_k.

This script does NOT prove anything; it checks that every algebraic/asymptotic
step the proof relies on holds numerically, so the written proof rests on
confirmed identities (the "verify before asserting" discipline). Each step prints
PASS/FAIL.
"""
import mpmath as mp

mp.mp.dps = 60

CASES = [(1, 0, 1), (2, 0, 1), (1, 0, 3), (1, 1, 1), (1, 3, 1), (3, 2, 5), (2, 3, 7)]


def roots(A, B, C):
    Af, Bf, Cf = mp.mpf(A), mp.mpf(B), mp.mpf(C)
    s = mp.sqrt(Bf * Bf - 4 * Af * Cf)            # mpc if disc < 0
    return (-Bf + s) / (2 * Af), (-Bf - s) / (2 * Af)


def b(A, B, C, k):
    return (mp.mpf(A) * k + mp.mpf(B)) * k + mp.mpf(C)


def P(A, B, C, n):
    pr = mp.mpf(1)
    for k in range(1, n + 1):
        pr *= b(A, B, C, k)
    return pr


def q_seq(A, B, C, N):
    qs = [mp.mpf(1), b(A, B, C, 1)]                # q_0=1, q_1=b_1
    for n in range(2, N + 1):
        qs.append(b(A, B, C, n) * qs[-1] + qs[-2])
    return qs


def R_inf(A, B, C, N=300000):
    Rm2, Rm1 = mp.mpf(1), mp.mpf(1)               # R_0, R_1
    for n in range(2, N + 1):
        Rn = Rm1 + Rm2 / (b(A, B, C, n - 1) * b(A, B, C, n))
        Rm2, Rm1 = Rm1, Rn
    return Rm1


def K_gamma(A, B, C):
    r1, r2 = roots(A, B, C)
    return -mp.log(mp.re(mp.gamma(1 - r1) * mp.gamma(1 - r2)))


TOL = mp.mpf(10) ** (-40)


def step1_product_gamma_identity():
    """P_n == A^n Gamma(n+1-r1)Gamma(n+1-r2)/(Gamma(1-r1)Gamma(1-r2)) exactly."""
    ok = True
    for (A, B, C) in CASES:
        r1, r2 = roots(A, B, C)
        denom = mp.gamma(1 - r1) * mp.gamma(1 - r2)
        for n in [1, 2, 3, 5, 8]:
            lhs = P(A, B, C, n)
            rhs = mp.mpf(A) ** n * mp.gamma(n + 1 - r1) * mp.gamma(n + 1 - r2) / denom
            if abs(mp.re(lhs - rhs)) > TOL * (1 + abs(mp.re(lhs))):
                ok = False
    return ok


def step2_gamma_ratio_rate():
    """g(n)=logGamma(n+1-r)-logGamma(n+1)+r log n = r(r-1)/(2n)+O(1/n^2):
    verify n*g(n) -> r(r-1)/2 (closed form) via 2-point Richardson (kills 1/n)."""
    ok = True
    for (A, B, C) in CASES:
        for r in roots(A, B, C):
            pred = r * (r - 1) / 2
            def ng(n):
                g = mp.loggamma(n + 1 - r) - mp.loggamma(n + 1) + r * mp.log(n)
                return mp.mpf(n) * g
            rich = 2 * ng(8000) - ng(4000)            # = pred + O(1/n^2)
            if abs(rich - pred) > mp.mpf('1e-6') * (1 + abs(pred)):
                ok = False
    return ok


def step3_product_constant_is_Kgamma():
    """[log P_n - (n log A + 2 log n! + (B/A) log n)] -> K_Gamma (product alone, NO delta)."""
    ok = True
    for (A, B, C) in CASES:
        kg = K_gamma(A, B, C)
        res = []
        for n in [2000, 4000, 8000]:
            skel = n * mp.log(A) + 2 * mp.loggamma(n + 1) + (mp.mpf(B) / A) * mp.log(n)
            res.append(mp.log(P(A, B, C, n)) - skel)
        # residual -> kg
        if abs(res[-1] - kg) > mp.mpf('1e-3'):
            ok = False
    return ok


def step4_R_monotone_bounded_rate():
    """R_n strictly increasing; R_inf - R_n = O(1/n^3) (n^3 * (R_inf - R_n) bounded)."""
    ok = True
    for (A, B, C) in CASES:
        Rinf = R_inf(A, B, C)
        Rm2, Rm1 = mp.mpf(1), mp.mpf(1)
        prev = Rm1
        scaled = []
        for n in range(2, 4001):
            Rn = Rm1 + Rm2 / (b(A, B, C, n - 1) * b(A, B, C, n))
            if Rn <= prev:                            # monotonicity
                ok = False
            prev = Rn
            Rm2, Rm1 = Rm1, Rn
            if n in (1000, 2000, 4000):
                scaled.append(mp.mpf(n) ** 3 * (Rinf - Rn))
        if Rinf <= 1:
            ok = False
        # n^3*(Rinf - Rn) should be bounded / converge to a constant, not grow
        if abs(scaled[-1] - scaled[-2]) > abs(scaled[-1]) * mp.mpf('0.1') + mp.mpf('1e-6'):
            ok = False
    return ok


def step5_full_law_error_O_1_over_n():
    """log q_n - skeleton - (K_Gamma+delta) = kappa/n + O(1/n^2) with the CLOSED FORM
    kappa = (1/2) sum_i r_i(r_i-1) = (1/2)(B^2/A^2 + B/A - 2C/A); the R_n side is
    O(1/n^3) so it does not touch the 1/n term. Verify n*resid -> kappa via Richardson."""
    ok = True
    for (A, B, C) in CASES:
        r1, r2 = roots(A, B, C)
        kappa = mp.re((r1 * (r1 - 1) + r2 * (r2 - 1)) / 2)
        kappa_alt = (mp.mpf(B) ** 2 / A ** 2 + mp.mpf(B) / A - 2 * mp.mpf(C) / A) / 2
        if abs(kappa - kappa_alt) > TOL:              # the two closed forms agree
            ok = False
        kt = K_gamma(A, B, C) + mp.log(R_inf(A, B, C))
        N = 8000
        qs = q_seq(A, B, C, N)
        def nresid(n):
            skel = n * mp.log(A) + 2 * mp.loggamma(n + 1) + (mp.mpf(B) / A) * mp.log(n)
            return mp.mpf(n) * (mp.log(qs[n]) - skel - kt)
        rich = 2 * nresid(8000) - nresid(4000)        # -> kappa + O(1/n^2)
        if abs(rich - kappa) > mp.mpf('1e-5') * (1 + abs(kappa)):
            ok = False
    return ok


if __name__ == "__main__":
    steps = [
        ("S1  exact product/Gamma identity  P_n = A^n Gamma(.)Gamma(.)/Gamma(.)Gamma(.)",
         step1_product_gamma_identity),
        ("S2  Gamma ratio  logGamma(n+1-r)-logGamma(n+1)+r log n = r(r-1)/(2n)+O(1/n^2)",
         step2_gamma_ratio_rate),
        ("S3  product constant = K_Gamma exactly (no delta in the product alone)",
         step3_product_constant_is_Kgamma),
        ("S4  R_n increasing, bounded, R_inf>1, R_inf-R_n = O(1/n^3)",
         step4_R_monotone_bounded_rate),
        ("S5  full law residual = kappa/n+O(1/n^2), kappa=(1/2)sum r_i(r_i-1) closed form",
         step5_full_law_error_O_1_over_n),
    ]
    print(f"Structural-proof step verification  [mpmath dps={mp.mp.dps}]  cases={CASES}\n")
    allok = True
    for label, fn in steps:
        res = fn()
        allok = allok and res
        print(f"  [{'PASS' if res else 'FAIL'}]  {label}")
    print(f"\nALL STEPS: {'PASS' if allok else 'FAIL'}")
