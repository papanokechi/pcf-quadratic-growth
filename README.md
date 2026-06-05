# pcf-quadratic-growth

Numeric + (later) Lean study of the convergent-denominator growth law for the
quadratic polynomial continued fraction

    V(A,B,C) = 1 + K_{n>=1} 1 / (A n^2 + B n + C),   integers A,B,C >= 1,

the B != 0 generalization of the deposited even (B=0) family. EXPLORATORY until a
core is machine-checked: numerics verify, they do not prove.

## Object

Standard continuant denominators
    q_0 = 1,  q_1 = A+B+C,  q_n = (A n^2 + B n + C) q_{n-1} + q_{n-2}.

Conjectured growth law
    log q_n = 2n log n - 2n + log(2 pi n) + n log A + (B/A) log n + K(A,B,C) + o(1).

The `(B/A) log n` term is the genuinely new feature versus the even family; it
vanishes at B=0.

## The constant (numerically established; see harness_quad.py)

    K_true(A,B,C) = K_Gamma(A,B,C) + delta(A,B,C)

  - K_Gamma = -log( Gamma(1-r1) Gamma(1-r2) ),  r1,r2 roots of A x^2+B x+C.
    The elementary naive-product constant prod_{k>=1} (b_k/(A k^2)).
    At B=0 it equals EXACTLY log( sinh(pi sqrt(C/A)) / (pi sqrt(C/A)) )
    (the deposited even constant; matched to 50 digits).
  - delta = log R_inf > 0,  R_inf = lim q_n / prod_{k=1}^n b_k,  a continued-
    fraction correction from the `+q_{n-2}` term. Two independent computations
    (Neville extrapolation; exact scaled recurrence R_n = R_{n-1} +
    R_{n-2}/(b_{n-1} b_n)) agree to ~15 digits.

## Epistemic status (four-grade SIARC convention)

- PROVEN  = clean Lean axiom cone {propext, Classical.choice, Quot.sound}, zero sorry.
- VERIFIED = numerically checked (mpmath/sympy) for the stated finite cases only.
- STRUCTURAL = a hand/symbolic argument, not machine-checked.
- CONJECTURED = believed, not yet established.

Current grades here:
- The full growth law  log q_n = n log A + 2 log n! + (B/A) log n + (K_Gamma+delta)
  + kappa/n + O(1/n^2):  STRUCTURAL (complete elementary hand proof in
  growth_law_note.tex) + VERIFIED numerically (verify_structural.py, 5 steps PASS).
  Upgraded from VERIFIED-numerically; PROVEN still pending a Lean core.
- K_Gamma = -log(Gamma(1-r1)Gamma(1-r2)): STRUCTURAL (exact product/Gamma identity
  P_n = A^n Gamma(n+1-r1)Gamma(n+1-r2)/(Gamma(1-r1)Gamma(1-r2)) + the standard
  Gamma-ratio asymptotic). At B=0 equals the sinh form (reflection formula).
- delta = log R_inf > 0: STRUCTURAL (scaled recurrence R_n monotone increasing,
  bounded by a convergent product, so R_n -> R_inf > 1; tail O(1/n^3)).
- kappa = (1/2)(B^2/A^2 + B/A - 2C/A): STRUCTURAL closed form for the 1/n term
  (= (1/2) sum r_i(r_i-1)); a bonus refinement, VERIFIED via Richardson (step S5).
- delta non-elementary in general: CONJECTURED (PSLQ null at 80 digits, (1,0,1)).

Standard input (confirmed): the only non-self-contained step is the Gamma-ratio
asymptotic in growth_law_note.tex Lemma 2 -- DLMF 5.11(iii) eqs. 5.11.13 and
5.11.17, equation numbers verified against the live DLMF release. No longer a
pre-deposit gate.

## RESOLVED (was STOP-AND-SURFACE): no tension with deposited work

The earlier worry -- that delta contradicts a deposited even-family claim of
"K = sinh-form" -- was checked against the actual deposited sources
(papanokechi/submitted-manuscripts @ commit c45e93b) and is a NON-issue:

- Decisive q_3 check: the corpus uses the STANDARD continuant
  Q_0=1, Q_{-1}=0, Q_n=(A n^2+B n+C)Q_{n-1}+Q_{n-2}. For A=C=1, B=0 this gives
  Q_3 = 112 (not the product-normalized 100), so delta = log R_inf is genuinely
  present in the deposited object.
- BUT no deposited paper asserts the exact additive constant. Both quadratic
  papers prove only a ONE-SIDED LOWER bound (log Q_n >~ n log n), which is all
  Euler's irrationality criterion (Wronskian = +/-1) needs:
    * arithmetic-sector-bifurcation-quadratic-continued-fractions/manuscript:
      "Q_n grows at least exponentially: log Q_n >~ n log n".
    * polynomial-continued-fractions-proved-logarithmic-laddera4casoratian/source.tex,
      Lemma "Denominator growth" (lem:Qgrowth): "log Q_n >= n log n - O(n)".
- The only "sinh" anywhere in the deposited corpus is the classical Stieltjes
  S-fraction VALUE K = coth(1) = I_{-1/2}(1)/I_{1/2}(1) (a Bessel-ratio value of
  one even CF, in the self-adjoint paper) -- NOT a q_n growth constant.

Conclusion: K_true = K_Gamma + delta is a genuinely NEW exact-constant result. It
is consistent with -- and sharpens, from a one-sided bound to an exact constant --
the deposited lower bound; it contradicts nothing deposited. NOT an erratum.
(The sinh form is simply K_Gamma at B=0, a correct naive-product identity; the
B=0 gate below verifies that identity, not a deposited q_n claim.)

## OPEN / out of scope

- No closed form for delta is asserted (a PSLQ search is expected to return null;
  a null result is a complete, honest answer).
- Lean formalization is deferred to where `lake` runs; only the finitary core
  (e.g. K_Gamma's product structure) is a formalization candidate.

## Layout

    harness_quad.py        numeric harness (run: python harness_quad.py)
    verify_structural.py   numeric confirmation of every proof step (S1..S5 PASS)
    growth_law_note.tex    the deposit write-up + STRUCTURAL proof (-> growth_law_note.pdf)
    CITATION.cff           citation metadata (from CITATION_quad.cff; DOIs FLAGGED)
    .zenodo.json           Zenodo deposit metadata (from zenodo_quad.json; DOIs FLAGGED)
    verify/                (later) Lean finitary core + cone_check
    LICENSE                Apache-2.0

Supplements the deposited even-quadratic family and the M9 corpus.
License: Apache-2.0.

## Deposit framing (what this artifact claims)

The deposited paper *"Polynomial Continued Fractions: a Proved Logarithmic Ladder,
a 4/pi Casoratian Identity, and 482 Irrational Constants"* (Papanokechi; source.tex
in papanokechi/submitted-manuscripts @ c45e93b) proves only the ONE-SIDED
denominator-growth bound `log Q_n >= n log n - O(n)` (lemma lem:Qgrowth), which is
all Euler's irrationality criterion needs for its 482-constant catalogue.

growth_law_note.tex (the deposit write-up) SHARPENS that bound, on the integer
sub-family A,B,C >= 1, to the exact two-sided law

    log q_n = n log A + 2 log(n!) + (B/A) log n + (K_Gamma + delta) + kappa/n + O(1/n^2),

recovering the deposited leading order `2 log(n!) = 2n log n - O(n)` as a corollary
while additionally pinning every lower-order term and the additive constant
`K = K_Gamma + delta`. It contradicts nothing deposited (no deposited manuscript
pins the constant). Honest scope: A,B,C >= 1 (the deposited bound covers the wider
range A n^2+B n+C > 0); and no closed form is claimed for delta.

Before deposit (FLAG operator): insert the deposited ladder paper's own DOI in the
bibliography -- do NOT reuse the M9 DOI 10.5281/zenodo.20542161 (different
artifact). This is the ONLY remaining pre-deposit FLAG; the DLMF equation numbers
(5.11.13, 5.11.17) are already confirmed against the live release. Then
create/confirm the github.com/papanokechi/pcf-quadratic-growth repo and set the tag
date. Commit/push/tag stay operator-gated (see setup_next_repos.ps1).
