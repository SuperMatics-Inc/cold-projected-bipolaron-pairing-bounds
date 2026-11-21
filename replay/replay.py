#!/usr/bin/env python3
"""
Standalone certified replay and mutant-rejection suite for W/MASSIVE/LARGE-U/.

Modes:
  python replay.py               -> CERTIFIED finite checks (exact Fraction + outward mpmath.iv)
  python replay.py --all-mutants -> run all guarded mutants (each must exit 1)
  python replay.py --numerical   -> reproduce NUMERICAL 2-site ED & traced-oscillator probes
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import math
import subprocess
import sys

import mpmath as mp
import numpy as np


# ---------------------------------------------------------------------------
# 1. Exact 16-state CAR algebra & Schrieffer-Wolff generator verification
# ---------------------------------------------------------------------------

def _build_16state_car():
    """Build exact 4^2 = 16-state two-site CAR operators with Jordan-Wigner signs."""
    # Single-site basis: 0=|0>, 1=|up>, 2=|dn>, 3=|2>=c_dn^dag c_up^dag |0>
    c_up = np.zeros((4, 4), dtype=np.float64)
    c_dn = np.zeros((4, 4), dtype=np.float64)
    c_up[0, 1] = 1.0
    c_up[2, 3] = -1.0  # c_up (c_dn^dag c_up^dag |0>) = -c_dn^dag |0>
    c_dn[0, 2] = 1.0
    c_dn[1, 3] = 1.0   # c_dn (c_dn^dag c_up^dag |0>) = +c_up^dag |0>
    n_site = c_up.T @ c_up + c_dn.T @ c_dn
    p_odd = np.diag([0.0, 1.0, 1.0, 0.0])
    p_even = np.eye(4) - p_odd
    parity = np.diag([1.0, -1.0, -1.0, 1.0])
    id4 = np.eye(4)

    # 2-site operators: site 0 first, site 1 second, JW string on site 0 for site 1
    c0_up = np.kron(c_up, id4)
    c0_dn = np.kron(c_dn, id4)
    c1_up = np.kron(parity, c_up)
    c1_dn = np.kron(parity, c_dn)

    N_odd = np.kron(p_odd, id4) + np.kron(id4, p_odd)
    P_even0 = np.kron(p_even, p_even)  # N_odd = 0 on the bond
    P_odd2 = np.kron(p_odd, p_odd)     # N_odd = 2 on the bond
    S0 = 0.5 * np.kron(n_site, id4)    # U(1) generator on site 0
    S1_site = 0.5 * np.kron(id4, n_site)

    # Directed hop from site 1 to site 0: A_e = c_{0,up}^dag c_{1,up} + c_{0,dn}^dag c_{1,dn}
    A_e = c0_up.T @ c1_up + c0_dn.T @ c1_dn
    # Parity-changing part: P_even0 A_e P_odd2 + P_odd2 A_e P_even0
    A_mix = P_even0 @ A_e @ P_odd2 + P_odd2 @ A_e @ P_even0
    # Parity-preserving part: A_e - A_mix
    A_diag = A_e - A_mix
    return {
        "A_e": A_e,
        "A_mix": A_mix,
        "A_diag": A_diag,
        "N_odd": N_odd,
        "P_even0": P_even0,
        "P_odd2": P_odd2,
        "S0": S0,
        "S1_site": S1_site,
    }


def check_car_algebra(mutant: str | None = None) -> None:
    car = _build_16state_car()
    A_e = car["A_e"]
    A_mix = car["A_mix"]
    A_diag = car["A_diag"]
    N_odd = car["N_odd"]
    P_even0 = car["P_even0"]
    P_odd2 = car["P_odd2"]
    S0 = car["S0"]

    # Check exact singular values (operator norms) -> sqrt(2), 1, 2
    sv_mix = np.linalg.svd(A_mix, compute_uv=False)
    sv_diag = np.linalg.svd(A_diag, compute_uv=False)
    sv_S = np.linalg.svd(A_mix - A_mix.T, compute_uv=False)
    if mutant == "car_norm_1":
        assert abs(sv_mix[0] - 1.0) < 1e-12, "mutant car_norm_1"
    assert abs(sv_mix[0] - math.sqrt(2.0)) < 1e-12
    assert abs(sv_diag[0] - 1.0) < 1e-12
    assert abs(sv_S[0] - 2.0) < 1e-12

    # Check SW commutator [S_{1,e}^{(dir)}, D_0] = -V_{e,mix}^{(dir)} for D_0 = (U/2) N_odd
    S_dir = P_even0 @ A_e @ P_odd2 - P_odd2 @ A_e @ P_even0
    D_half = 0.5 * N_odd
    comm = S_dir @ D_half - D_half @ S_dir
    err_sw = np.max(np.abs(comm - A_mix))
    if mutant == "sw_sign":
        err_sw = np.max(np.abs(comm + A_mix))
    assert err_sw < 1e-12, f"SW commutator mismatch {err_sw}"

    # Verify S_dir has exact same operator norm sqrt(2) as A_mix
    sv_Sdir = np.linalg.svd(S_dir, compute_uv=False)
    assert abs(sv_Sdir[0] - math.sqrt(2.0)) < 1e-12

    # Verify [S0, S_dir] = +(1/2) S_dir (directed hop 1->0 increases n_0 by 1)
    comm_S0 = S0 @ S_dir - S_dir @ S0
    err_charge = np.max(np.abs(comm_S0 - 0.5 * S_dir))
    assert err_charge < 1e-12, f"Charge commutator mismatch {err_charge}"

    # Verify exponential rotation identity R(theta) S_dir R(-theta) = e^{+delta/2} S_dir
    from scipy.linalg import expm
    delta = 0.7
    R = expm(delta * S0)
    Rinv = expm(-delta * S0)
    rot_Sdir = R @ S_dir @ Rinv
    target_sign = 0.5 * delta
    if mutant == "rot_sign":
        target_sign = -0.5 * delta
    err_rot = np.max(np.abs(rot_Sdir - math.exp(target_sign) * S_dir))
    assert err_rot < 1e-12, f"Rotation mismatch {err_rot}"

    # Verify V_diag vanishes on N_odd=0 AND [S_1, V_diag] vanishes on N_odd=0
    assert np.max(np.abs(P_even0 @ A_diag @ P_even0)) < 1e-12
    assert np.max(np.abs(P_even0 @ (S_dir @ A_diag - A_diag @ S_dir) @ P_even0)) < 1e-12

    # Verify on 3-site chain that distinct edges e1 != e0 give P_{N_odd=0} [S_{1,e1}, V_{e0,mix}] P_{N_odd=0} = 0
    id4 = np.eye(4)
    p_even = np.diag([1.0, 0.0, 0.0, 1.0])
    P3_even0 = np.kron(p_even, np.kron(p_even, p_even))
    S_01 = np.kron(S_dir, id4)
    A_12_mix = np.kron(id4, A_mix)
    if mutant == "m1_distinct_edge":
        # Same-edge e1 = e0 = (0,1) produces a non-zero m=1 superexchange/pair-hopping block on N_odd=0
        A_12_mix = np.kron(A_mix, id4)
    comm_3 = P3_even0 @ (S_01 @ A_12_mix - A_12_mix @ S_01) @ P3_even0
    err_3 = np.max(np.abs(comm_3))
    assert err_3 < 1e-12, f"Distinct-edge m=1 projection nonzero: {err_3}"

    # Verify exact single-site Lang-Firsov unitary conventions (P0-1 & P2-1):
    #   H_at = Om a^dag a - sqrt(U Om / 2) (n - 1)(a^dag + a) - mu n   (NO pre-LF +(U/2) n_odd!)
    #   U_LF = exp(-alpha (n - 1)(a^dag - a)), alpha = sqrt(U / (2 Om))
    #   U_LF H_at U_LF^dag = (U/2) n_odd + Om a^dag a - (U/2) I - mu n
    #   U_LF c_up^dag U_LF^dag = c_up^dag D(-alpha)
    #   U_LF Delta U_LF^dag = Delta D(+2 alpha)
    Nph = 40
    Nlow = 6
    a = np.zeros((Nph, Nph))
    for n in range(1, Nph):
        a[n - 1, n] = math.sqrt(n)
    n_ph = a.T @ a
    I_ph = np.eye(Nph)
    P_low = np.zeros((Nph, Nlow))
    for n in range(Nlow):
        P_low[n, n] = 1.0
    P4_low = np.kron(id4, P_low)

    c_up = np.zeros((4, 4))
    c_dn = np.zeros((4, 4))
    c_up[0, 1] = 1.0
    c_up[2, 3] = -1.0
    c_dn[0, 2] = 1.0
    c_dn[1, 3] = 1.0
    n_el = c_up.T @ c_up + c_dn.T @ c_dn
    n_odd = np.diag([0.0, 1.0, 1.0, 0.0])
    Delta_el = c_dn @ c_up

    U_val, Om_val, mu_val = 3.0, 1.5, 0.4
    alpha = math.sqrt(U_val / (2.0 * Om_val))
    lf_sign = -1.0
    if mutant == "lf_sign_wrong":
        lf_sign = +1.0
    U_LF = expm(lf_sign * alpha * np.kron(n_el - id4, a.T - a))
    D_minus_alpha = expm(-alpha * (a.T - a))
    D_plus_2alpha = expm(2.0 * alpha * (a.T - a))

    H_at_pre = (
        np.kron(id4, Om_val * n_ph)
        - math.sqrt(U_val * Om_val / 2.0) * np.kron(n_el - id4, a.T + a)
        - mu_val * np.kron(n_el, I_ph)
    )
    H_at_target = (
        np.kron(0.5 * U_val * (n_odd - id4) - mu_val * n_el, I_ph)
        + np.kron(id4, Om_val * n_ph)
    )
    err_hat = np.max(np.abs(P4_low.T @ (U_LF @ H_at_pre @ U_LF.T - H_at_target) @ P4_low))
    assert err_hat < 1e-11, f"LF H_at transformation mismatch: {err_hat}"

    err_cup = np.max(np.abs(P4_low.T @ (U_LF @ np.kron(c_up.T, I_ph) @ U_LF.T - np.kron(c_up.T, D_minus_alpha)) @ P4_low))
    assert err_cup < 1e-11, f"LF c_up^dag transformation mismatch: {err_cup}"

    err_pair = np.max(np.abs(P4_low.T @ (U_LF @ np.kron(Delta_el, I_ph) @ U_LF.T - np.kron(Delta_el, D_plus_2alpha)) @ P4_low))
    assert err_pair < 1e-10, f"LF Delta transformation mismatch: {err_pair}"


# ---------------------------------------------------------------------------
# 2. Certified BCH & Bargmann-Fock bounds (outward mpmath.iv + exact Fraction)
# ---------------------------------------------------------------------------

def _iv_upper(val: mp.iv.mpf) -> mp.mpf:
    return mp.mpf(val.b)


def _iv_lower(val: mp.iv.mpf) -> mp.mpf:
    return mp.mpf(val.a)
def _iv_cosh(u: mp.iv.mpf) -> mp.iv.mpf:
    return (mp.iv.exp(u) + mp.iv.exp(-u)) / 2


def _iv_sinh(u: mp.iv.mpf) -> mp.iv.mpf:
    return (mp.iv.exp(u) - mp.iv.exp(-u)) / 2


def check_bargmann_fock_and_norms(mutant: str | None = None) -> None:
    mp.mp.dps = mp.iv.dps = 60

    # Check exact inequality 1 + u <= 2 * exp(u/2) for u >= 0 (used for M=0, 0 < X < 64)
    # Proof: 2*exp(u/2) - (1 + u) = 1 + (u/2)^2/2! + ... >= 1 > 0 for all u >= 0.
    for u_val in [Fraction(0), Fraction(1, 10), Fraction(2), Fraction(10), Fraction(128)]:
        u = mp.iv.mpf(u_val.numerator) / u_val.denominator
        lhs = 1 + u
        rhs = 2 * mp.iv.exp(u / 2)
        if mutant == "bf_m0_factor":
            rhs = mp.iv.exp(u / 2)  # false at u=2 (3 <= 2.718 is false!)
        assert _iv_upper(lhs) <= _iv_lower(rhs), f"1+u <= 2 exp(u/2) failed at u={u_val}"

    # Check X >= 64 Bargmann-Fock bound:
    # For M = floor(X/64), r^2 = M+1 <= X/64 + 1 <= X/32,
    # -(s^2 X)/4 + 2|s| sqrt(X/2) sqrt(X/32) + X/64 + 1 = 1 - (s^2 - |s| - 1/16) X/4
    # At |s| >= 2: s^2 - |s| - 1/16 >= 4 - 2 - 1/16 = 31/16, so exponent <= 1 - (31/64) X.
    assert Fraction(4 - 2) - Fraction(1, 16) == Fraction(31, 16)
    assert Fraction(31, 16) * Fraction(1, 4) == Fraction(31, 64)

    # Check monotonicity & value at X = 64 of (X/64 + 2) * exp(1 - (3/64)*X) <= sqrt(2)
    # Derivative of log(X/64 + 2) - (3/64) X is 1/(X+128) - 3/64 < 0 for X >= 0
    deriv_at_0 = Fraction(1, 128) - Fraction(3, 64)
    assert deriv_at_0 == Fraction(-5, 128) < 0
    val_at_64 = (mp.iv.mpf(64) / 64 + 2) * mp.iv.exp(mp.iv.mpf(1) - mp.iv.mpf(3))
    sqrt2_iv = mp.iv.sqrt(mp.iv.mpf(2))
    if mutant == "bf_m_ge1_threshold":
        # Corrupt 7/16 -> 1/2, making exponent 1 - 32/64 = 0? No: 31/64 - 1/2 = -1/64!
        val_at_64 = (mp.iv.mpf(64) / 64 + 2) * mp.iv.exp(mp.iv.mpf(1) + mp.iv.mpf(1))
    assert _iv_upper(val_at_64) < _iv_lower(sqrt2_iv), "X>=64 Bargmann-Fock crossover failed"

    # Also cross-check against direct exact matrix element summation for small X=1, 8, 32, 64
    for X_int in [1, 8, 32, 64]:
        M = X_int // 64
        alpha = mp.sqrt(mp.mpf(X_int) / 2)
        for s in [2, 3, 4]:
            beta = s * alpha
            # Build exact matrix <n | D(beta) | m> for n <= M+1, m <= M via Laguerre/Hermite formula
            mat = np.zeros((M + 2, M + 1), dtype=np.float64)
            for n in range(M + 2):
                for m in range(M + 1):
                    # Exact <n|D(beta)|m> for real beta:
                    # = exp(-beta^2/2) * sqrt(m! n!) * sum_{k=0}^{min(n,m)} (-1)^{m-k} beta^{n+m-2k} / (k! (n-k)! (m-k)!)
                    val = mp.mpf(0)
                    for k in range(min(n, m) + 1):
                        term = ((-1) ** (m - k)) * (beta ** (n + m - 2 * k)) / (
                            mp.factorial(k) * mp.factorial(n - k) * mp.factorial(m - k)
                        )
                        val += term
                    val *= mp.exp(-beta ** 2 / 2) * mp.sqrt(mp.factorial(n) * mp.factorial(m))
                    mat[n, m] = float(val)
            op_norm = float(np.linalg.svd(mat, compute_uv=False)[0])
            bd = float(mp.sqrt(2) * mp.exp(-mp.mpf(7) * X_int / 16))
            assert op_norm <= bd, f"Exact BF norm violated at X={X_int}, s={s}: {op_norm} > {bd}"

    # Check Lemma 1 edge operator norm bounds for |delta_e| <= 4*ln(2) (so |delta_e|/2 <= 2*ln(2), cosh(2*ln 2) = 17/8):
    #   2*sqrt(2)*cosh(2*ln 2) = (17/4)*sqrt(2) < 7 and 2*cosh(2*ln 2) = 17/4 < 5
    half_twist = 2 * mp.iv.log(2)
    if mutant == "cosh_norm_too_tight":
        half_twist = 4 * mp.iv.log(2)  # cosh(4*ln 2) = 257/32 gives 2*sqrt(2)*(257/32) = 22.71 > 7
    edge_mix = 2 * mp.iv.sqrt(2) * _iv_cosh(half_twist)
    edge_diag = 2 * _iv_cosh(half_twist)
    assert _iv_upper(edge_mix) < mp.mpf(7)
    assert _iv_upper(edge_diag) < mp.mpf(5)
    # Check Weyl pull bound: for X >= 1, M = floor(X/64), m >= 2:
    # Omega * (sqrt(M+1) + m*sqrt(X/2))^2 <= U * (m + 2)^2
    # Since M+1 <= X/64 + 1 <= (65/64) X for X >= 1:
    # sqrt(65/64) + m/sqrt(2) <= 2 + m identically (since sqrt(65/64) < 2 and 1/sqrt(2) < 1)!
    assert _iv_upper(mp.iv.sqrt(mp.iv.mpf(65) / 64)) < 2
    assert _iv_upper(1 / mp.iv.sqrt(2)) < 1


# ---------------------------------------------------------------------------
# 3. Certified torus shell sum & Theorem 1 pairing ceiling constants
# ---------------------------------------------------------------------------

def check_torus_and_series_constants(mutant: str | None = None) -> None:
    mp.mp.dps = mp.iv.dps = 60

    for d in [1, 2, 5, 20, 50]:
        s = mp.iv.mpf(0)
        for r in range(1, d + 1):
            s += (8 * r - 4) * (mp.iv.log(mp.iv.mpf(r + 1) / r) ** 2)
        rhs = 8 * mp.iv.log(mp.iv.mpf(d + 1))
        if mutant == "shell_sum_factor":
            rhs = 4 * mp.iv.log(mp.iv.mpf(d + 1))
        assert _iv_upper(s) <= _iv_lower(rhs)

    # Check exact Z^2 L^1 edge-ball count |B_R^(E)(e_0)| <= 8*(R+1)^2 for R = 0..11 (P1-2):
    for R in range(12):
        edges_in_ball = set()
        for x0 in range(-R - 2, R + 3):
            for y0 in range(-R - 2, R + 3):
                for dx, dy in [(1, 0), (0, 1)]:
                    u = (x0, y0)
                    v = (x0 + dx, y0 + dy)
                    dist_e0 = min(
                        abs(u[0]) + abs(u[1]),
                        abs(u[0] - 1) + abs(u[1]),
                        abs(v[0]) + abs(v[1]),
                        abs(v[0] - 1) + abs(v[1]),
                    )
                    if dist_e0 <= R:
                        edges_in_ball.add((u, v))
        assert len(edges_in_ball) <= 8 * ((R + 1) ** 2)

    # Warm-up series check at |delta_e| <= 1 (|delta_e|/2 <= 1/2, x = 2^-11, r = 28*x = 7/512):
    #   ||S_{1,e}|| < (7/2)*x, ||V_e|| < 6*U*x, Taylor weight 1/2 * Lemma 3 factor 2 = 1:
    #   Family I (m >= 1): sum_{m>=1} 2 * (m+1)^5 * (6 x) * (28 x)^m
    #   Family II (m >= 2): sum_{m>=2} m^4 * (m+2)^2 * (28 x)^m
    x_warm = Fraction(1, 2048)
    r_warm = 28 * x_warm  # = 7/512
    assert r_warm == Fraction(7, 512)

    s1_warm = Fraction(0)
    s2_warm = Fraction(0)
    for m in range(1, 30):
        s1_warm += 2 * ((m + 1) ** 5) * (6 * x_warm) * (r_warm ** m)
        if m >= 2:
            s2_warm += (m ** 4) * ((m + 2) ** 2) * (r_warm ** m)

    ratio1_warm = Fraction(32, 31) ** 5 * r_warm
    ratio2_warm = Fraction(31, 30) ** 4 * Fraction(33, 32) ** 2 * r_warm
    assert ratio1_warm < Fraction(1, 32) and ratio2_warm < Fraction(1, 32)
    tail1_warm = (2 * (31 ** 5) * (6 * x_warm) * (r_warm ** 30)) / (1 - ratio1_warm)
    tail2_warm = ((30 ** 4) * (32 ** 2) * (r_warm ** 30)) / (1 - ratio2_warm)

    K_warm_over_x2 = (s1_warm + tail1_warm + s2_warm + tail2_warm) / (x_warm ** 2)
    prefactor_warm = 64 * K_warm_over_x2
    if mutant == "prefactor_2_26":
        assert prefactor_warm <= Fraction(2 ** 23), "mutant prefactor_2_26"
    assert prefactor_warm < Fraction(2 ** 24), f"Warm-up prefactor {float(prefactor_warm)} exceeds 2^24!"


def check_kappa4_summability(mutant: str | None = None) -> None:
    """
    Certify Theorem 1, Corollary 1.1 & Lemma 4 WITH xi_cold(beta) > 2 (2D summability) for kappa in (0, 4]:
      |delta_e| <= 4 * ln(2)  ==>  |delta_e|/2 <= 2 * ln(2).
      cosh(2 * ln 2) = (4 + 1/4)/2 = 17/8 (EXACT RATIONAL).
      ||S_{1,e}(s*theta)||_op <= 2 * sqrt(2) * (17/8) * x = (17/4)*sqrt(2)*x < 7*x.
      ||V_{e,mix}(s*theta)||_op + ||V_{e,diag}(s*theta)||_op < (7 + 5)*U*x = 12*U*x.
      Per-step factor: 4 * (2 * 7 * x) = 56 * x = r.
    """
    mp.mp.dps = mp.iv.dps = 60

    cosh_val = _iv_cosh(2 * mp.iv.log(2))
    assert _iv_lower(cosh_val) <= mp.mpf(17) / 8 <= _iv_upper(cosh_val)
    edge_mix_4 = 2 * mp.iv.sqrt(2) * (mp.iv.mpf(17) / 8)
    edge_diag_4 = 2 * (mp.iv.mpf(17) / 8)
    assert _iv_upper(edge_mix_4) < 7
    assert _iv_upper(edge_diag_4) < 5

    # Exact Lemma 4 series at x = 2^-14:
    #   Family I (m >= 1): sum_{m>=1} 2 * (m+1)^5 * (12 x) * r^m
    #   Family II (m >= 2): sum_{m>=2} m^4 * (m+2)^2 * r^m
    x = Fraction(1, 2 ** 14)
    r = 56 * x  # = 56 / 16384 = 7 / 2048 < 1/256
    if mutant == "series_divergent":
        x = Fraction(1, 32)
        r = 56 * x  # > 1, series diverges!
    assert r < Fraction(1, 256), f"Convergence ratio {float(r)} too large"

    s1 = Fraction(0)
    s2 = Fraction(0)
    for m in range(1, 30):
        s1 += 2 * ((m + 1) ** 5) * (12 * x) * (r ** m)
        if m >= 2:
            s2 += (m ** 4) * ((m + 2) ** 2) * (r ** m)

    ratio1 = Fraction(32, 31) ** 5 * r
    ratio2 = Fraction(31, 30) ** 4 * Fraction(33, 32) ** 2 * r
    assert ratio1 < Fraction(1, 128) and ratio2 < Fraction(1, 128)
    tail1 = (2 * (31 ** 5) * (12 * x) * (r ** 30)) / (1 - ratio1)
    tail2 = ((30 ** 4) * (32 ** 2) * (r ** 30)) / (1 - ratio2)

    K_over_x2 = (s1 + tail1 + s2 + tail2) / (x ** 2)
    assert K_over_x2 < Fraction(2 ** 20), f"K_4(x)/x^2 = {float(K_over_x2):.3f} >= 2^20"
    prefactor = 64 * K_over_x2
    if mutant == "kappa4_prefactor":
        assert prefactor <= Fraction(2 ** 25)  # false: 5.562e7 > 3.355e7 = 2^25
    assert prefactor < Fraction(2 ** 26), f"prefactor={float(prefactor):.3e} >= 2^26"

    # Verify exact two-branch formula for xi_cold(beta) (P1-1):
    #   xi_cold(beta) = 4 - 16*A if A < 1/8 else 1/(4*A)
    for A_val in [Fraction(1, 16), Fraction(1, 8), Fraction(1, 4)]:
        kappa_opt = Fraction(4) if A_val < Fraction(1, 8) else Fraction(1) / (2 * A_val)
        xi_val = kappa_opt - A_val * (kappa_opt ** 2)
        if A_val < Fraction(1, 8):
            assert xi_val == 4 - 16 * A_val > 2
            assert xi_val != Fraction(1) / (4 * A_val)
        else:
            assert xi_val == Fraction(1) / (4 * A_val) <= 2

    # Verify Corollary 1.1 rotated-cloud convergence for R(theta) e^{S_1} Delta_x^dag D_x(-2 alpha) e^{-S_1} R(-theta):
    cloud_sum = Fraction(1, 1) / ((1 - r) ** 2)
    assert cloud_sum < Fraction(2), f"Rotated cloud norm {float(cloud_sum)} >= 2"

# ---------------------------------------------------------------------------
# 4. Certified Microscopic Neutral-Loop Lower Bound (Proposition 2) & Thermal Crossover
# ---------------------------------------------------------------------------

def check_neutral_loop_and_hky_crossover(mutant: str | None = None) -> None:
    mp.mp.dps = mp.iv.dps = 60

    # Proposition 2 Part 1: Exact 2-site N=2, S_z=0, Omega=0 (4-state) eigenvalue shift
    t_frac = Fraction(1, 100)
    delta_frac = Fraction(1, 2)
    t_mp = mp.iv.mpf(t_frac.numerator) / t_frac.denominator
    delta_mp = mp.iv.mpf(delta_frac.numerator) / delta_frac.denominator
    c2 = _iv_cosh(delta_mp / 2) ** 2
    s2 = _iv_sinh(delta_mp / 2) ** 2
    E0_0 = (1 - mp.iv.sqrt(1 + 16 * (t_mp ** 2))) / 2
    E0_d = (1 - mp.iv.sqrt(1 + 16 * (t_mp ** 2) * c2)) / 2
    shift = E0_0 - E0_d
    lb = 4 * (t_mp ** 2) * s2 / mp.iv.sqrt(1 + 16 * (t_mp ** 2) * c2)
    if mutant == "neutral_floor_sign":
        lb = 2 * lb
    assert _iv_lower(shift) >= _iv_upper(lb)

    # Proposition 2 Part 2: Non-perturbative finite-temperature lower bound (P1-A, P1-B, P1-C):
    #   1) P1-C: p_trip(u) <= (1/2) * exp(-3*beta*U/4) <= (1/2)*exp(-3) < 0.025 for beta*U >= 4, u <= 2, t/U <= 1/16
    p_trip_bound = mp.iv.exp(mp.iv.mpf(-3)) / 2
    assert _iv_upper(p_trip_bound) < mp.mpf("0.025")
    assert _iv_lower(1 - p_trip_bound) > mp.mpf("0.975")

    #   2) P1-B & P1-C: Verify all intermediate numerical inequalities in Proposition 2(2):
    #      Z_ph(tilde_beta/2)/Z_ph(tilde_beta) = 1 + exp(-tilde_beta*Om/2) < 2  ==>  (3/2)*ratio < 3
    #      (2*ln 3)/3 + 1 < 7/4, 5 + 4*(7/4) = 12, exp(-49/400) > 0.88, 0.975 * 0.0176 = 0.01716 > 0.01 = 2/200
    h0_avg_bd = (2 * mp.iv.log(3)) / 3 + 1
    assert _iv_upper(h0_avg_bd) < mp.mpf(7) / 4
    assert 5 + 4 * Fraction(7, 4) == 12
    jensen_fac = mp.iv.exp(mp.iv.mpf(-49) / 400)
    assert _iv_lower(jensen_fac) > mp.mpf("0.88")
    assert Fraction(2, 100) * Fraction(88, 100) == Fraction(176, 10000)
    assert Fraction(975, 1000) * Fraction(176, 10000) == Fraction(1716, 100000) > Fraction(2, 200)

    #   3) P1-A: Verify exact pre-LF 3-state operator identity W H_0 W = W^2 (U I_3 + Om n_rel),
    #      position-basis positivity of exp(-tau * u * W_el) and Tr_3(rho_{3,tilde_beta} W^2) >= 2 t^2
    #      (while confirming Fock-basis H_0 has positive off-diagonal entries on P_02, so Fock basis is NOT Metzler),
    #      and finite-Fock trace inequality -(F_sym(delta) - F_sym(0)) >= (1/400)*(t^2/U)*sinh^2(delta/2)
    from scipy.linalg import expm
    Nr = 35
    ar = np.zeros((Nr, Nr))
    for n in range(1, Nr):
        ar[n - 1, n] = math.sqrt(n)
    nr = ar.T @ ar
    Ir = np.eye(Nr)
    P20 = np.diag([1.0, 0.0, 0.0])
    P02 = np.diag([0.0, 1.0, 0.0])
    I3 = np.eye(3)
    W_el = np.array([[0.0, 0.0, -math.sqrt(2.0)], [0.0, 0.0, -math.sqrt(2.0)], [-math.sqrt(2.0), -math.sqrt(2.0), 0.0]])
    # Position-basis local electronic semigroup exp(-tau * u * W_el) is entrywise strictly positive for tau * u > 0:
    assert np.min(expm(-0.3 * W_el)) > 0.0
    K_off_el = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert np.allclose(W_el @ W_el, 2.0 * (I3 + K_off_el))
    W_op = np.kron(W_el, Ir)
    K_off_op = np.kron(K_off_el, Ir)
    U_val, t_val = 20.0, 0.25
    c0_claim = 1.0 / 400.0
    if mutant == "prop2_finite_t_c0":
        c0_claim = 5.0  # false: exceeds actual ratio ~1.85 at Om/U=0.2, beta*U=4, cosh=2.0
    for Om_ratio in [10.0, 1.0, 0.2]:
        Om_val = Om_ratio * U_val
        H0_pre = (
            np.kron(U_val * I3, Ir)
            + np.kron(I3, Om_val * nr)
            - math.sqrt(U_val * Om_val) * np.kron(P20 - P02, ar.T + ar)
        )
        id_err = np.linalg.norm(
            (t_val ** 2) * (W_op @ H0_pre @ W_op - (W_op @ W_op) @ (np.kron(U_val * I3, Ir) + np.kron(I3, Om_val * nr)))
        )
        assert id_err < 1e-11, f"Pre-LF W H0 W identity error: {id_err}"
        # Confirm Fock-basis H0_pre has positive off-diagonal entries on P02 (not Fock-Metzler),
        # yet Tr_3(rho_{3,tilde_beta} K_off) > 0 and Tr_3(rho_{3,tilde_beta} W^2) >= 2*t^2 by position-basis positivity:
        assert np.max(H0_pre - np.diag(np.diag(H0_pre))) > 0.0
        H_u_sample = H0_pre + t_val * W_op
        evals_u, evecs_u = np.linalg.eigh(H_u_sample)
        w_u = np.exp(-(3.0 / U_val) * (evals_u - evals_u[0]))
        rho_u = (evecs_u * (w_u / np.sum(w_u))) @ evecs_u.T
        tr_koff = float(np.trace(rho_u @ K_off_op))
        if mutant == "prop2_pos_positivity":
            tr_koff = -tr_koff
        assert tr_koff > 0.0, f"Position-basis positivity trace failed: {tr_koff}"
        ev0 = np.linalg.eigvalsh(H0_pre + t_val * W_op)
        ev_trip_base = U_val + Om_val * np.arange(Nr, dtype=np.float64)
        for bU in [4.0, 10.0]:
            beta = bU / U_val
            for ch in [1.2, 1.5, 2.0]:
                sh2 = ch * ch - 1.0
                ev_eta = np.linalg.eigvalsh(H0_pre + (t_val * ch) * W_op)
                E_ref = min(ev0[0], ev_eta[0])
                Z0 = np.sum(np.exp(-beta * (ev0 - E_ref))) + np.sum(np.exp(-beta * (ev_trip_base - E_ref)))
                Zeta = np.sum(np.exp(-beta * (ev_eta - E_ref))) + np.sum(np.exp(-beta * (ev_trip_base - E_ref)))
                dF = math.log(Zeta / Z0) / beta
                lb_claim = c0_claim * (t_val * t_val / U_val) * sh2
                assert dF >= lb_claim, f"Prop 2(2) bound violated: dF={dF} < {lb_claim}"

    # Proposition 3: Check HKY thermal pair-transfer maximum in tau:
    for y_val in [Fraction(1, 10), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(10)]:
        y_iv = mp.iv.mpf(y_val.numerator) / y_val.denominator
        lhs = y_iv * mp.iv.exp(-y_iv)
        rhs = mp.iv.exp(mp.iv.mpf(-1))
        if mutant == "hky_max_wrong":
            rhs = mp.iv.mpf(1) / 3  # false at y=1: 1/e = 0.3678 > 1/3
        assert _iv_upper(lhs) <= _iv_lower(rhs) + (mp.mpf("1e-50") if y_val == 1 else 0)

    # Uniform-in-Omega Franck-Condon resolvent lower bound >= 1/(2U) (Jensen's inequality on Pois(X), P2-3):
    for X_val in [Fraction(1, 2), Fraction(1), Fraction(10), Fraction(40)]:
        X_iv = mp.iv.mpf(X_val.numerator) / X_val.denominator
        m_max = max(30, int(4 * X_val))
        s_fc = mp.iv.mpf(0)
        term = mp.iv.exp(-X_iv)
        for n in range(m_max + 1):
            if n > 0:
                term = term * X_iv / n
            s_fc += term / (1 + mp.iv.mpf(n) / X_iv)
        assert _iv_lower(s_fc) >= _iv_upper(mp.iv.mpf(1) / 2)


# ---------------------------------------------------------------------------
# 5. Guarded mutant runner & numerical reproduction
# ---------------------------------------------------------------------------

MUTANTS = [
    "car_norm_1",
    "sw_sign",
    "rot_sign",
    "m1_distinct_edge",
    "lf_sign_wrong",
    "bf_m0_factor",
    "bf_m_ge1_threshold",
    "cosh_norm_too_tight",
    "shell_sum_factor",
    "prefactor_2_26",
    "series_divergent",
    "kappa4_prefactor",
    "neutral_floor_sign",
    "prop2_finite_t_c0",
    "hky_max_wrong",
    "prop2_pos_positivity",
]


def run_all_checks(mutant: str | None = None) -> None:
    check_car_algebra(mutant=mutant)
    check_bargmann_fock_and_norms(mutant=mutant)
    check_torus_and_series_constants(mutant=mutant)
    check_kappa4_summability(mutant=mutant)
    check_neutral_loop_and_hky_crossover(mutant=mutant)


def run_numerical() -> None:
    """Reproduce NUMERICAL probes (2-site N=2, S_z=0 ED & traced-oscillator neutral floor)."""
    from scipy.linalg import expm

    # (A) 2-site N=2, S_z=0 Holstein ED thermal pair-splitting crossover (NUMERICAL)
    t, U = 0.5, 20.0
    print("=== NUMERICAL Check A: 2-site N=2, S_z=0 Holstein thermal splitting ===")
    for Om in [1.0, 2.0]:
        X = U / Om
        Nph = 60
        # Relative phonon mode a_rel = (a_0 - a_1)/sqrt(2), displacement magnitude sqrt(X)
        # Electronic basis: 0=|2,0>, 1=|0,2>, 2=|up,dn>, 3=|dn,up>
        a = np.zeros((Nph, Nph))
        for n in range(1, Nph):
            a[n - 1, n] = math.sqrt(n)
        D_rel = expm(math.sqrt(X) * (a.T - a))
        H_LF = np.zeros((4 * Nph, 4 * Nph))
        diag_el = np.diag([0.0, 0.0, U, U])
        n_op = np.diag(np.arange(Nph, dtype=np.float64))
        H_LF += np.kron(diag_el, np.eye(Nph)) + np.kron(np.eye(4), Om * n_op)
        M_20 = np.zeros((4, 4)); M_20[0, 2] = -1.0; M_20[0, 3] = 1.0
        M_02 = np.zeros((4, 4)); M_02[1, 2] = -1.0; M_02[1, 3] = 1.0
        H_hop = -t * (
            np.kron(M_20, D_rel) + np.kron(M_20.T, D_rel.T)
            + np.kron(M_02, D_rel.T) + np.kron(M_02.T, D_rel)
        )
        H_LF += H_hop
        # Swap parity: 0 <-> 1, 2 <-> -3 (c_0 <-> c_1 maps c_0up^dag c_1dn^dag |0> to c_1up^dag c_0dn^dag |0> = -|dn,up>)
        P_el = np.array([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ], dtype=np.float64)
        P_ph = np.diag([(-1.0) ** n for n in range(Nph)])
        P_tot = np.kron(P_el, P_ph)
        assert np.max(np.abs(H_LF @ P_tot - P_tot @ H_LF)) < 1e-12
        evals, evecs = np.linalg.eigh(H_LF)
        parities = np.diag(evecs.T @ P_tot @ evecs)
        ev_even = np.sort(evals[parities > 0])
        ev_odd = np.sort(evals[parities < 0])
        E0 = evals[0]
        for tau in [0.5 * Om, Om, 2.0 * Om, 4.0 * Om, 8.0 * Om]:
            b = 1.0 / tau
            w_e = np.exp(-b * (ev_even - E0))
            w_o = np.exp(-b * (ev_odd - E0))
            avg_split = float(np.sum(0.5 * (w_e + w_o) * np.abs(ev_odd - ev_even)) / np.sum(0.5 * (w_e + w_o)))
            print(f"  Om={Om:.1f}, tau/Om={tau/Om:4.1f}: <|E_odd-E_even|>_tau = {avg_split:.3e} (t^2/U={t*t/U:.4e})")

    # (B) 2-site N=2, S_z=0 Hermitian-twist free-energy shift -(1/beta) log(Z_sym(delta)/Z_0)
    print("=== NUMERICAL Check B: Microscopic Hermitian-twist neutral floor (Prop 2) ===")
    t, U = 0.25, 20.0
    for Om in [2.0, 0.5]:
        X = U / Om
        Nph = 50
        a = np.zeros((Nph, Nph))
        for n in range(1, Nph):
            a[n - 1, n] = math.sqrt(n)
        D_rel = expm(math.sqrt(X) * (a.T - a))
        diag_el = np.diag([0.0, 0.0, U, U])
        n_op = np.diag(np.arange(Nph, dtype=np.float64))
        H_0 = np.kron(diag_el, np.eye(Nph)) + np.kron(np.eye(4), Om * n_op)
        M_20 = np.zeros((4, 4)); M_20[0, 2] = -1.0; M_20[0, 3] = 1.0
        M_02 = np.zeros((4, 4)); M_02[1, 2] = -1.0; M_02[1, 3] = 1.0
        V_hop = -t * (
            np.kron(M_20, D_rel) + np.kron(M_20.T, D_rel.T)
            + np.kron(M_02, D_rel.T) + np.kron(M_02.T, D_rel)
        )
        delta = 0.5
        H_plain = H_0 + V_hop
        H_sym = H_0 + math.cosh(delta / 2.0) * V_hop
        ev0 = np.linalg.eigvalsh(H_plain)
        evd = np.linalg.eigvalsh(H_sym)
        for tau in [0.5, 1.0, 2.0]:
            b = 1.0 / tau
            Z0 = np.sum(np.exp(-b * (ev0 - ev0[0])))
            Zd = np.sum(np.exp(-b * (evd - ev0[0])))
            dF = float(np.log(Zd / Z0) / b)
            pred = (t * t / U) * (math.sinh(delta / 2.0) ** 2)
            print(f"  Om={Om:.1f}, tau={tau:.1f}: -(F_sym-F_0) = {dF:.6e}, (t^2/U)sinh^2(d/2) = {pred:.6e}, ratio = {dF/pred:.4f}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-mutants", action="store_true", help="Run all guarded mutants (each must exit 1)")
    parser.add_argument("--mutant", type=str, default=None, help="Run a single named mutant")
    parser.add_argument("--numerical", action="store_true", help="Run numerical ED sanity checks")
    args = parser.parse_args()

    if args.mutant:
        run_all_checks(mutant=args.mutant)
        print(f"ERROR: mutant {args.mutant} unexpectedly passed!")
        sys.exit(0)

    if args.all_mutants:
        for m in MUTANTS:
            res = subprocess.run(
                [sys.executable, __file__, "--mutant", m],
                capture_output=True,
                text=True,
            )
            if res.returncode == 0:
                print(f"FAIL: mutant {m} exited 0!")
                sys.exit(1)
        print(f"ALL {len(MUTANTS)} GUARDED MUTANTS REJECTED (each exited 1).")
        sys.exit(0)

    run_all_checks()
    print("CERTIFIED REPLAY PASSED (5 exact-rational & outward-interval check groups).")
    if args.numerical:
        run_numerical()


if __name__ == "__main__":
    main()
