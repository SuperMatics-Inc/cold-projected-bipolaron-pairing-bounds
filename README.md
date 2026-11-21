# Complex-rotation pairing bounds for a cold-projected bipolaron-phonon effective model and obstructions to microscopic extensions in the two-dimensional Holstein model

Davis Rens, SuperMatics Inc. ([ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)).
Version 1.0.0.

Paper DOI: <https://doi.org/10.5281/zenodo.24415833>
Replay supplement DOI: <https://doi.org/10.5281/zenodo.24415835>
Source: <https://github.com/SuperMatics-Inc/cold-projected-bipolaron-pairing-bounds>

## What this deposit is

This repository carries the preprint `paper/main.pdf` (20 pages) and the replay supplement for it. The paper
applies a charge rotation (complex rotation) to the cold-projected Lang-Firsov plus first-order
Schrieffer-Wolff effective model of the spinful two-dimensional Holstein model, proves a pairing bound for that
effective model, and then proves two obstructions that locate the boundary of the method: a two-site floor for
a rotation applied directly to the microscopic Hamiltonian, and the vanishing Gibbs weight of the cold manifold.

The reviewed manuscript is released byte for byte. `\SELFDOI` is left empty in `paper/main.tex` and the text
reads "the record reserved for it"; the DOI of this record appears only in this file, in `CITATION.cff` and in
the deposit metadata, so the deposited PDF is exactly the document that was reviewed.

## What is proved

Write `x = t/U`, `X = U/Omega`, `M = floor(X/64)`, and let `H_cold` be the projection of the effective
Hamiltonian onto the sector `N_odd = 0`, `n^ph_z <= M`. All statements hold for `x <= 2^-14` and `X >= 1`.

1. **Theorem 1 (cold-projected pairing bound).** The polaron-pair correlation of the Gibbs state of `H_cold`
   obeys

   `|<Delta^dagger_{x_0} Delta_{y_0}>_cold| <= (d+1)^{-xi_cold(beta)}`,

   where `d` is the torus distance and

   `xi_cold(beta) = 4 - 16 A(beta)` for `A(beta) < 1/8`, and `xi_cold(beta) = 1/(4 A(beta))` otherwise,

   with `A(beta) = 8 beta U K_4(x) e^{-7X/8}` and `K_4(x) <= 2^20 x^2`. Hence the correlations are summable in
   two dimensions whenever `tau > 2^26 (t^2/U) e^{-7U/(8 Omega)}`, uniformly in the torus side `L >= 4` and in
   the chemical potential, and the clustering temperature of the effective state obeys
   `k_B T^cold_pair <= 2^26 (t^2/U) e^{-7U/(8 Omega)}`.
2. **Corollary 1.1 (original pair field).** In the same effective state the original equal-time electronic pair
   field obeys the same bound with prefactor `4`, because two rotated clouds each cost less than a factor `2`.
3. **Proposition 2 (neutral-loop floor for a direct microscopic rotation).** Already on two sites (`N = 2`,
   `S_z = 0`) virtual neutral excursions lower the Hermitian-twist free energy by at least
   `(1/400)(t^2/U) sinh^2(delta/2)` for every `Omega > 0`, `beta U >= 4` and `0 < t/U <= 1/16`, and the
   zero-temperature Franck-Condon resolvent obeys `sum_n e^{-X} X^n / (n! (U + n Omega)) >= 1/(2U)` for every
   `Omega > 0`.
4. **Proposition 3 (thermal transfer and the Gibbs weight of the cold manifold).** The Boltzmann factor
   `beta (t^2/U) e^{-U/tau}` peaks at `tau = U` with the value `t^2/(e U^2)`, and at every `mu`, `beta > 0` and
   `L >= 4` the Gibbs weight of the cold manifold is at most `(1 - e^{-beta Omega (M+1)})^{|Lambda|}`, which
   behaves like `exp(-|Lambda| e^{-beta Omega (M+1)})`.
5. **Open.** Control of the full Holstein Gibbs state for phonon occupations above `M` is not proved here, nor
   is the fixed-scale problem `k_B T_pair <= C Omega`.

## What this adds, and what it does not claim

The narrowed ceiling `2^26 (t^2/U) e^{-7U/(8 Omega)}` sits below the released full-Gibbs ceiling `2^27 t^2/U`
of the companion preprint (<https://doi.org/10.5281/zenodo.24414616>) for every `X > (8/7) ln(2^26/6000)`,
about `10.65`. It is not an improvement of that theorem: the object bounded here is an effective model, not
the Holstein Gibbs state, and Proposition 3 quantifies why the projection cannot be removed for free.

The limitations are stated in the paper and hold here as well.

- The bounds are on a clustering temperature of one two-point function inside an effective state. They are not
  bounds on a critical temperature, no phase transition is proved or excluded, and nothing here asserts that a
  material is a superconductor.
- Theorem 1 and Corollary 1.1 are statements about `H_cold`. The reading that a rotation applied directly to
  the microscopic Hamiltonian cannot reach a ceiling below a constant multiple of `t^2/U` on a lattice is the
  mechanism identified by the two-site floor of Proposition 2; a lattice-wide lower bound over all rotation
  profiles and all fillings is not proved, and no second theorem is claimed.
- The prefactor `2^26` is the dominating power of two of the certified chain, not an optimum, and no optimality
  is claimed for `2^26`, for `K_4` or for the exponent `xi_cold`.
- The certificate program proves finite arithmetic, not analysis. It does not replace any operator inequality,
  each of which has a written derivation.

## Contents

```
MANIFEST.sha256          SHA-256 of every tracked file except this one
CITATION.cff             citation metadata for the preprint
LICENSE                  MIT License, for the code
paper/main.pdf           the reviewed manuscript, released unmodified
paper/main.tex           LaTeX source of that manuscript
paper/refs.bib           bibliography
replay/replay.py         the standalone certificate program
replay/README.md         what the program verifies and how to run it
replay/SHA256SUMS        SHA-256 of replay.py and replay/README.md, pinned since review
```

`paper/main.pdf` is the reviewed artifact and is released byte for byte; its SHA-256 digest is
`b92ee523c4947cce05e9f1b7906cf6abe47685568982759c4ec185c98df7633c` and `paper/main.tex` is
`54577845049ab85c2c3cb1e5a0acc9acbaae88526510a695ce155a07445ac84a`. Continuous integration verifies those
bytes, the text they carry and the replay program; it does not rebuild the PDF from `main.tex`, so the deposit
is exactly the reviewed document.

## Replay

`replay/replay.py` is one standalone program, no data file and no network access. It uses `mpmath` for
outward-rounded interval arithmetic, `numpy` for the finite-level operators, and `scipy` for the matrix
exponentials of the two-site checks; the proof steps it certifies are in exact rational arithmetic or in
outward-rounded intervals, never in floating point.

```
pip install mpmath numpy scipy
python replay/replay.py
python replay/replay.py --all-mutants
python replay/replay.py --numerical
```

The default run checks five groups of statements and ends with
`CERTIFIED REPLAY PASSED (5 exact-rational & outward-interval check groups).`: the sixteen-state operator
identities on a bond and the rotation identity; the Bargmann-Fock and Poisson bounds with the edge operator
norms; the torus shell telescoping and the two series that give `K_4(x) <= 2^20 x^2` together with
`64 K_4(x)/x^2 < 2^26`; the rotated-cloud convergence of Corollary 1.1; and the two-site bounds of
Proposition 2. The `--all-mutants` run corrupts the guarded claims one at a time, sixteen variants in all, each
in a separate child process, and requires every variant to exit `1`; it ends with
`ALL 16 GUARDED MUTANTS REJECTED (each exited 1).` The `--numerical` run adds the two two-site
diagonalisations described in `replay/README.md`. The program proves arithmetic, not analysis.

`replay/README.md` is the guide pinned since review: it records the commands and what each one certifies.

## Licences

The code in this repository, that is `replay/replay.py` and the workflow, is released under the MIT License
(`LICENSE`). The manuscript text, `paper/main.pdf`, `paper/main.tex` and `paper/refs.bib`, is released under
the Creative Commons Attribution 4.0 International license (CC-BY-4.0), as on the paper record.

## Status and disclosures

AI agents drafted the proofs, the code and the drafts under D. Rens' direction, and separate AI systems carried
out the independent reviews of the manuscript and of the certificate program. This preprint has not been peer
reviewed. It is a companion note to the released full-Gibbs ceilings and records an obstruction, not a
resolution of it.
