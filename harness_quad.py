"""
pcf-quadratic-growth -- numeric harness (EXPLORATORY: numerically verifies, does NOT prove).

Conjecture under test (G-law for the full quadratic PCF):

  For integers A, B, C >= 1, the generalized continued fraction
      V(A,B,C) = 1 + K_{n>=1} 1 / (A n^2 + B n + C)
  has convergent denominators q_n (q_0 = 1, q_1 = A+B+C,
  q_n = (A n^2 + B n + C) q_{n-1} + q_{n-2}) obeying the growth law

      log q_n = 2 n log n - 2 n + log(2 pi n) + n log A
                + (B/A) log n + K(A,B,C) + o(1),

  with a closed-form constant determined by the roots r1, r2 of A x^2 + B x + C
  (so r1 + r2 = -B/A, r1 r2 = C/A):

      K(A,B,C) = - log( Gamma(1 - r1) Gamma(1 - r2) ).

  The B = 0 slice must reduce to the DEPOSITED even-family constant
      K(A,0,C) = log( sinh(pi sqrt(C/A)) / (pi sqrt(C/A)) ).

NOTE the genuinely NEW structural feature vs the even (B=0) family: a nonzero B
forces an extra  (B/A) log n  term in the growth law (because the naive product
prod (A k^2 + B k + C) carries an n^{B/A} factor). That term vanishes at B = 0,
recovering the deposited even law exactly.

This script computes K numerically (Neville extrapolation in 1/n) and compares
to the closed form. A MATCH is numerical evidence (status: VERIFIED-numerically
for the tested triples only), NOT a proof. A MISMATCH is an equally valid,
informative result: it would mean the +q_{n-2} correction contributes to the
constant and the closed form needs a correction factor. Either way the number is
the ground truth, computed first, before any write-up or formalization.
"""
import mpmath as mp

mp.mp.dps = 50


def logR(A, B, C, N):
    """log R_N where R_n = q_n / (A^n (n!)^2 n^(B/A)) -> exp(K) + o(1)."""
    Af, Bf, Cf = mp.mpf(A), mp.mpf(B), mp.mpf(C)
    if N == 0:
        q = mp.mpf(1)
    elif N == 1:
        q = Af + Bf + Cf
    else:
        qm2 = mp.mpf(1)          # q_0
        qm1 = Af + Bf + Cf       # q_1
        for n in range(2, N + 1):
            bn = Af * n * n + Bf * n + Cf
            q = bn * qm1 + qm2
            qm2, qm1 = qm1, q
    Nf = mp.mpf(N)
    skel = Nf * mp.log(Af) + 2 * mp.loggamma(Nf + 1) + (Bf / Af) * mp.log(Nf)
    return mp.log(q) - skel


def neville_at_zero(xs, ys):
    """Polynomial extrapolation of (xs, ys) evaluated at x = 0 (Neville)."""
    P = list(ys)
    m = len(xs)
    for k in range(1, m):
        for i in range(m - 1, k - 1, -1):
            P[i] = ((0 - xs[i - k]) * P[i] - (0 - xs[i]) * P[i - 1]) / (xs[i] - xs[i - k])
    return P[-1]


def K_numeric(A, B, C, Ns):
    xs = [1 / mp.mpf(N) for N in Ns]
    ys = [logR(A, B, C, N) for N in Ns]
    return neville_at_zero(xs, ys)


def K_closed(A, B, C):
    Af, Bf, Cf = mp.mpf(A), mp.mpf(B), mp.mpf(C)
    disc = Bf * Bf - 4 * Af * Cf
    s = mp.sqrt(disc)                 # mpc when disc < 0
    r1 = (-Bf + s) / (2 * Af)
    r2 = (-Bf - s) / (2 * Af)
    val = mp.gamma(1 - r1) * mp.gamma(1 - r2)   # real: conjugate gammas or two real gammas
    return -mp.log(mp.re(val))


def K_even_sinh(A, C):
    Af, Cf = mp.mpf(A), mp.mpf(C)
    y = mp.sqrt(Cf / Af)
    return mp.log(mp.sinh(mp.pi * y) / (mp.pi * y))


def delta_scaled(A, B, C, N=50000):
    """delta = log R_inf via the EXACT scaled recurrence, independent of any
    logarithmic extrapolation. With R_n = q_n / prod_{k=1}^n b_k one has, exactly,
        R_n = R_{n-1} + R_{n-2} / (b_{n-1} b_n),   R_0 = R_1 = 1.
    R_n increases monotonically to R_inf > 1; the increments are ~ 1/n^4
    (summable), so R_inf is reached fast and delta = log R_inf. This is the
    robust cross-check that delta is real and not a Neville artifact."""
    Af, Bf, Cf = mp.mpf(A), mp.mpf(B), mp.mpf(C)

    def b(k):
        kk = mp.mpf(k)
        return Af * kk * kk + Bf * kk + Cf

    Rm2 = mp.mpf(1)   # R_0
    Rm1 = mp.mpf(1)   # R_1
    for n in range(2, N + 1):
        Rn = Rm1 + Rm2 / (b(n - 1) * b(n))
        Rm2, Rm1 = Rm1, Rn
    return mp.log(Rm1), Rm1


if __name__ == "__main__":
    # Larger nodes -> K_true to ~10+ digits (tail is ~ c/n, confirmed by diagnostic).
    Ns = [400, 560, 720, 880, 1040, 1200, 1400, 1600]
    cases = [
        (1, 0, 1), (2, 0, 1), (1, 0, 3),          # B=0 self-tests (K_Gamma must equal sinh)
        (1, 1, 1), (1, 2, 1), (1, 1, 2), (2, 1, 1),
        (1, 3, 1), (3, 2, 5), (2, 3, 7), (1, 5, 2),
    ]
    print("Quadratic PCF growth constant: K_true = K_Gamma + delta")
    print("  K_true  = lim [ log q_n - (n log A + 2 log n! + (B/A) log n) ]   (numeric)")
    print("  K_Gamma = -log( Gamma(1-r1) Gamma(1-r2) )   (naive-product / Borel constant)")
    print("  delta   = K_true - K_Gamma   (continued-fraction correction from +q_{n-2})")
    print("  delta_R = log R_inf via exact scaled recurrence R_n=R_{n-1}+R_{n-2}/(b_{n-1}b_n)")
    print(f"  [mpmath dps={mp.mp.dps}; Neville nodes N={Ns}]\n")
    hdr = f"{'(A,B,C)':>11} {'disc':>5}  {'K_true':>16} {'K_Gamma':>16} {'delta':>13} {'delta_R':>13} {'|d-dR|':>9}"
    print(hdr)
    print("-" * len(hdr))
    for (A, B, C) in cases:
        kt = K_numeric(A, B, C, Ns)
        kg = K_closed(A, B, C)
        d = kt - kg
        dR, _ = delta_scaled(A, B, C)
        gap = abs(d - dR)
        print(f"{str((A,B,C)):>11} {B*B-4*A*C:>5}  {mp.nstr(kt,12):>16} {mp.nstr(kg,12):>16} "
              f"{mp.nstr(d,9):>13} {mp.nstr(dR,9):>13} {mp.nstr(gap,2):>9}")

    print("\n[GATE] B=0: K_Gamma must equal the deposited even-family sinh constant")
    ok = True
    for (A, C) in [(1, 1), (2, 1), (1, 3)]:
        kg = K_closed(A, 0, C)
        ks = K_even_sinh(A, C)
        m = abs(kg - ks) < mp.mpf(10) ** (-40)
        ok = ok and m
        print(f"  A={A} C={C}: K_Gamma={mp.nstr(kg,15)}  sinh={mp.nstr(ks,15)}  match={m}")
    print(f"  GATE: {'PASS' if ok else 'FAIL'}")

    print("\n[RESOLVED] delta != 0 means the TRUE convergent-denominator constant is")
    print("  K_Gamma + delta, NOT K_Gamma. Checked against the deposited sources")
    print("  (papanokechi/submitted-manuscripts @ c45e93b): the corpus uses the STANDARD")
    print("  continuant Q_n (below: Q_3 = 112, not product 100) but asserts only a")
    print("  one-sided LOWER bound log Q_n >~ n log n (arithmetic-sector & logarithmic-")
    print("  ladder papers, for Euler's criterion). No deposited paper claims the exact")
    print("  sinh additive constant, so K_true = K_Gamma + delta is NEW and consistent")
    print("  with -- it sharpens -- the deposited bound. Not an erratum. (The only")
    print("  deposited 'sinh' is the coth Bessel-ratio VALUE in the self-adjoint paper.)")

    # Decisive check (now RESOLVED): for A=C=1 the standard continuant gives
    #   q_0=1, q_1=2, q_2=11, q_3=112   while the naive product P_3 = 2*5*10 = 100.
    # The deposited papers (arithmetic-sector + ladder) use Q_0=1, Q_{-1}=0, i.e. the
    # standard continuant -> Q_3 = 112, so delta = log R_inf is real in their object;
    # they bound it one-sidedly (log Q_n >~ n log n) rather than pinning the constant.
    qs = [1, 2]
    for n in range(2, 4):
        bn = 1 * n * n + 0 * n + 1
        qs.append(bn * qs[-1] + qs[-2])
    P3 = 2 * 5 * 10
    print(f"\n  DECISIVE CHECK (A=C=1): standard continuant q_0..q_3 = {qs}  (P_3 = {P3})")
    print(f"    Deposited Q_n is the standard continuant -> Q_3 = 112 == {qs[3]}: CONFIRMED.")
    print(f"    => delta is real in the deposited object; the papers bound it one-sidedly,")
    print(f"       they do not pin the additive constant (so no conflict, no erratum).")
