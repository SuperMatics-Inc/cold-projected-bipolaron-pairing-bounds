# CLAIMS: what the preprint asserts, and on what evidence

Vocabulary, as in Section 1 of the paper: *proved* means a written derivation in
the manuscript; *certified* means a finite arithmetic identity or inequality on
which the derivation rests is recomputed in exact rational or outward-rounded
interval arithmetic by `replay/replay.py`; *numerical* means a finite-matrix
diagonalisation; *open* means neither. Labels below are the paper's; the source
`paper/main.tex` carries each claim under the named label, and the replay check
groups are those printed by the program.

## C1. Theorem 1 (`thm:main`): cold-projected pairing ceiling with `e^{-7U/(8\Omega)}` narrowing

**Proved, with a certified replay.**

- C1a. The rotated cold Hamiltonian expands in a norm-convergent series on the
  cold manifold for `|\theta_u - \theta_v| <= 4 ln 2` (BCH and Lie series, with
  cluster bounds). The constants `17/8 = cosh(2 ln 2)`,
  `||S_{1,e}||_op < 7x`, `||ad_{S_{1,e}}||_op < 14x`,
  `||V_mix|| + ||V_diag|| < 12Ux` and `r = 56x` lie in the certified chain
  (default run, check groups 1, 2 and 3).
- C1b. Charge-diagonal monomials carry no twist, and the `m = 0, 1` terms and the
  `U(1)`-invariant terms cancel identically. Proved.
- C1c. The cold projection between phonon shells `M` and `M+1` costs
  `sqrt 2 e^{-7X/16}` per site (Bargmann-Fock), so two displacement operators
  give `2 e^{-7X/8}`. Proved; certified in outward-rounded interval arithmetic.
- C1d. Correlations decay as `(d+1)^{-xi_cold}` with
  `xi_cold = 4 - 16 A(beta)` for `A(beta) < 1/8` and `xi_cold = 1/(4 A(beta))`
  otherwise, where `A(beta) = 8 beta U K_4(x) e^{-7X/8}` and
  `K_4(x) <= 2^20 x^2` (`eq:A`, `eq:decay`, `eq:K4def`). Proved; the series
  constant is certified at `x = 2^-14` as
  `64 K_4(x)/x^2 ~ 5.562e7 < 2^26 ~ 6.711e7`.
- C1e. If `tau > 2^26 (t^2/U) e^{-7U/(8\Omega)}` then the polaron-pair
  correlations are uniformly 2D summable in the torus side `L >= 4` and in `mu`,
  so `k_B T^cold_pair <= 2^26 (t^2/U) e^{-7U/(8\Omega)}` in the sense of
  `def:tpair` (`eq:ceiling`, `eq:summ`). Proved. The threshold is the algebraic
  decay `A_0 (d(x,y)+1)^{-xi}` with `xi > 2` at fixed `X = U/\Omega >= 1` and
  `A_0`, `xi` uniform in `L` and `mu`; the exponential-decay and
  `\Omega`-uniform statement belongs to the full-Gibbs ceiling and is not
  claimed for the cold state.
- Scope: `x = t/U <= 2^-14`, `X = U/\Omega >= 1`, and the state is the Gibbs
  state of the projected effective Hamiltonian `H_cold`, not of the Holstein
  Hamiltonian. The prefactor `2^26` is the dominating power of two of the
  certified chain, not an optimum.

## C2. Corollary 1.1 (`cor:orig`): the original electronic pair field in the cold state

Proved under the hypotheses of C1. The cold-projected original pair field obeys
`4 (d+1)^{-xi_cold}`, so the identical ceiling gives 2D summability of the
original electronic pair correlations inside the cold manifold (`eq:corbound`).
The prefactor `4` comes from two rotated clouds of ratio `r = 56x`, that is
`(1-r)^{-2} < 2` per cloud, certified in check group 4.

## C3. Proposition 2 (`prop:floor`): neutral-loop floor for a direct microscopic rotation

Proved as a two-site statement, with a numerical two-site cross-check
(`eq:fc`, `eq:gain`, `eq:hubbard`).

- C3a. Franck-Condon resolvent floor: for all `\Omega > 0`, `U > 0`,
  `sum_n e^{-X} X^n / (n! (U + n\Omega)) >= 1/(2U)` (Jensen). Proved; certified.
- C3b. Free-energy gain of the Hermitian twist: for all `\Omega > 0`,
  `beta U >= 4`, `0 < t/U <= 1/16` and `eta = cosh(\delta/2)` in `[1, 2]`,
  `-(F_sym(\delta) - F_sym(0)) >= (1/400)(t^2/U) sinh^2(\delta/2)`. Proved;
  certified (the chain `e^{-1/400-12/100} = e^{-49/400} > 0.88` and
  `\Phi_3' >= 0.0176 beta (t^2/U) u` are checked in exact rational arithmetic).
- C3c. Zero-displacement attractive-Hubbard reference splitting
  `E_0(t) - E_0(t cosh(\delta/2)) >= 4 t^2 sinh^2(\delta/2) / sqrt(U^2 + 16 t^2 cosh^2(\delta/2))`;
  certified as an algebraic identity.
- C3d. Numerical cross-check: five independent two-site parameter points across
  `X < 1` and `X > 1`, spanning the allowed `t/U` and `beta U`, satisfy C3a and
  C3b with wide margin; the positivity statistic holds in every row. These rows
  come from an independent two-site implementation, not from `replay.py --numerical`,
  and the paper says so in the text and caption.
- Scope. The proved content is the two-site floor at `N = 2`, `S_z = 0`,
  `mu = 0`, for every `\Omega` and every `beta U >= 4`. The reading that a
  rotation applied directly to the microscopic Hamiltonian cannot reach a
  ceiling below a constant multiple of `t^2/U` on a lattice is the *mechanism*
  identified by the floor: the paper marks it as heuristic and states in Section
  4 and in the Limitations that a lattice-wide lower bound over all rotation
  profiles and all fillings is not proved. No second theorem is claimed.

## C4. Proposition 3 (`prop:thermal`): thermal pair transfer and the Gibbs weight of the cold manifold

Proved, with a certified peak identity and a numerical crossover check
(`eq:peak`, `eq:weight`).

- C4a. Peak identity: for all `t, U > 0`,
  `max_{\tau>0} beta (t^2/U) e^{-U/\tau} = t^2/(e U^2)`, attained at
  `\tau = U`. Proved; certified in exact rational arithmetic. In the effective
  thermal pair hopping the exponent crosses over from `2U/\Omega`
  (`\tau << \Omega`) to `U/\tau` (`\tau >> \Omega`).
- C4b. Cold-manifold Gibbs weight: for the atomic oscillator Hamiltonian, at
  every `mu`, `beta > 0` and `L >= 4`, the ratio is at most
  `(1 - e^{-beta\Omega(M+1)})^{|\Lambda|} <= exp(-|\Lambda| e^{-beta\Omega(M+1)})`,
  so the projection cannot be undone by restricting the Gibbs trace, and it
  gives nothing at fixed temperature once occupations above `M` are admitted.
  Proved.
- C4c. Numerical: the Boltzmann-averaged parity splitting at `t = 0.5`,
  `U = 20`, `\Omega = 1` rises from `1.40e-10` (`\tau = 0.5\Omega`) to
  `2.02e-2` (`\tau = 8\Omega`), reproduced by `replay.py --numerical`.

## C5. Open: control of the full Gibbs state for `n^ph_z > M`

The paper does not lift Theorem 1 to `Tr(e^{-beta H_\Lambda})`. Admitting
occupations above `M` in the Gibbs state reintroduces the neutral excursions
whose free-energy gain Proposition 2 bounds below, and the two mechanisms have
not been reconciled in one estimate. Whether the full Gibbs pair correlations
obey an exponentially narrowed ceiling in `U/\Omega`, or a sub-`t^2/U` ceiling
such as `C\Omega`, or neither, is not decided here. This is also the substance
of the title's obstructions, and it is what the site record and the deposit
descriptions call open.

## Evidence trail

- `paper/main.tex` is the reviewed source: `sha256
  54577845049ab85c2c3cb1e5a0acc9acbaae88526510a695ce155a07445ac84a`.
- `paper/main.pdf` is the reviewed build: `sha256
  b92ee523c4947cce05e9f1b7906cf6abe47685568982759c4ec185c98df7633c`.
- `replay/replay.py` (`sha256
  4de1a9451dab5715a74001bad155b33da4a8b11f88b474dc9ed43640036d54fb`) prints
  `CERTIFIED REPLAY PASSED (5 exact-rational & outward-interval check groups).`
  in its default mode and, with `--all-mutants`, corrupts sixteen guarded
  claims one at a time and prints `ALL 16 GUARDED MUTANTS REJECTED (each exited 1).`
  It certifies finite arithmetic only: it does not formalise the operator
  inequalities or the analytic steps in infinite dimensions.
- `MANIFEST.sha256` pins every released file; the hosted CI reruns both replay
  modes on the released tree and reads the released PDF.
