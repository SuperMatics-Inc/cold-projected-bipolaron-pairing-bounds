# Replay supplement

Certificate program for the companion note *Complex-rotation pairing bounds for a
cold-projected bipolaron-phonon effective model and obstructions to microscopic
extensions in the 2D Holstein model*.

## Contents

| File | Description |
| --- | --- |
| `replay.py` | The certificate program. Standalone, single file, no data files. |
| `README.md` | This file. |
| `SHA256SUMS` | SHA-256 digests of `replay.py` and `README.md`. |

## Requirements

Python 3.10 or newer, with `sympy`, `mpmath`, `numpy` and `scipy`. No other
package and no network access is needed. The program writes nothing to disk.

## Commands and expected results

Run from this directory. With `uv`, which installs the dependencies into an
isolated cache:

```
UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with sympy --with mpmath --with numpy --with scipy python -B replay.py
UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with sympy --with mpmath --with numpy --with scipy python -B replay.py --all-mutants
UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with sympy --with mpmath --with numpy --with scipy python -B replay.py --numerical
```

or, in an environment where the four packages are installed:

```
python -B replay.py
python -B replay.py --all-mutants
python -B replay.py --numerical
```

All three commands must exit `0`.

### 1. Default mode

Checks five groups of statements in exact rational arithmetic or in
outward-rounded interval arithmetic: the sixteen-state operator identities on a
bond and the rotation identity, the Bargmann-Fock and Poisson bounds together
with the edge operator norms, the torus shell telescoping and the two series that
give `K_4(x) <= 2^20 x^2` together with `64 K_4(x)/x^2 < 2^26`, the rotated-cloud
convergence of the corollary, and the two-site bounds of the obstruction
proposition. Expected final line:

```
CERTIFIED REPLAY PASSED (5 exact-rational & outward-interval check groups).
```

### 2. Mutant mode

Corrupts the claimed inequalities one at a time, sixteen variants in all, and
requires every variant to be rejected. Each variant is run as a separate
subprocess with the same interpreter, so this mode takes about sixteen times as
long as the default mode. Expected final line:

```
ALL 16 GUARDED MUTANTS REJECTED (each exited 1).
```

### 3. Numerical mode

Runs the default mode, then two diagonalisations of the two-site problem.

Check A is the thermal splitting of the two-site `N=2`, `S_z=0` spectrum at
`t=0.5`, `U=20`, phonon frequency `Omega=1` and `2`, and `tau/Omega` from `0.5`
to `8`. At `Omega=1` the Boltzmann-averaged parity splitting rises from
`1.402e-10` at `tau/Omega=0.5` to `2.018e-2` at `tau/Omega=8`.

Check B is the Hermitian-twist free-energy shift of the same two-site problem at
`t=0.25`, `U=20`, twist `delta=0.5`, and `(Omega,tau)` in
`{(2.0,0.5),(2.0,1.0),(2.0,2.0),(0.5,0.5),(0.5,1.0),(0.5,2.0)}`. The printed
`ratio` is the free-energy gain divided by `(t^2/U) sinh^2(delta/2)`; it lies
between `1.0267` and `1.1662`, above the proved lower bound `1/400` of the note
(the written chain records the stronger intermediate `1/200` first and then
weakens it).

## What is certified and what is not

The program certifies finite arithmetic identities and inequalities: operator
norms and traces on finite-dimensional spaces, exact rational prefactors, the
summation of two convergent series plus a geometric tail, and finite matrix
spectra. It does not formalise the operator inequalities of the paper, which rest
on the written argument, and it does not check the analytic steps in the
infinite-dimensional Hilbert space.

## Digest

```
4de1a9451dab5715a74001bad155b33da4a8b11f88b474dc9ed43640036d54fb  replay.py
```

The full list is in `SHA256SUMS`; verify it with `sha256sum -c SHA256SUMS` or with
Python's `hashlib`.
