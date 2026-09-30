"""Numerical cross-check of Theorem D (two-window families); not part of the proof.

(1) For every pair of one-step configurations a, a' and every pair j < k with a_j != a'_j and a_k != a'_k, the coefficient
    w_jk = c_{a_k-a_j} + c_{a'_k-a'_j} - c_{a'_k-a_j} - c_{a_k-a'_j} is >= 2 sec(pi/2d), with equality attained (d = 2..30).
(2) For random two-window families V_k = i^{-a_k}(1-P_k) + i^{-a'_k}P_k with random (non-commuting) projections P_k in M_M,
    F(V) <= F_DKZ (d = 3..9, M = 2..4).
Usage:  python check_two_window.py
"""
import itertools
import time
import numpy as np

t0 = time.time()


def consts(d):
    psi = lambda m: np.pi * m / (2 * d)
    h = {m: 1 / np.cos(psi(m)) - 1j / np.sin(psi(m)) for m in range(1, d)}
    c = lambda t, m: (1 / np.cos(psi(m)), 1 / np.sin(psi(m)), -1 / np.cos(psi(m)), -1 / np.sin(psi(m)))[t % 4]
    FD = sum((d - m) / np.cos(psi(m)) for m in range(1, d))
    return h, c, FD


def one_step(d):
    out = []
    for cc in range(4):
        for r in range(d):
            a = tuple((cc + (1 if k >= r else 0)) % 4 for k in range(d))
            if a not in out:
                out.append(a)
    return out


print('(1) smallest coefficient w_jk over all pairs of one-step configurations and all mixed pairs j < k')
ok1 = True
for d in range(2, 31):
    h, c, FD = consts(d)
    wins = one_step(d)
    Fa = lambda a: sum((h[k - j] * 1j ** ((a[k] - a[j]) % 4)).real for j in range(d) for k in range(j + 1, d))
    assert all(abs(Fa(a) - FD) < 1e-9 for a in wins)
    wmin = np.inf
    for a, ap in itertools.combinations(wins, 2):
        mixed = [k for k in range(d) if a[k] != ap[k]]
        for j, k in itertools.combinations(mixed, 2):
            m = k - j
            w = c(a[k] - a[j], m) + c(ap[k] - ap[j], m) - c(ap[k] - a[j], m) - c(a[k] - ap[j], m)
            wmin = min(wmin, w)
    ref = 2 / np.cos(np.pi / (2 * d))
    ok1 &= abs(wmin - ref) < 1e-9
    print(f'   d = {d:2d}: min w = {wmin:.9f},  2 sec(pi/2d) = {ref:.9f}')

print('(2) random two-window families with random projections: max F(V) - F_DKZ')
rng = np.random.default_rng(20260930)
ok2 = True
for d in range(3, 10):
    h, c, FD = consts(d)
    wins = one_step(d)
    worst = -np.inf
    for trial in range(60):
        M = int(rng.integers(2, 5))
        a, ap = wins[rng.integers(len(wins))], wins[rng.integers(len(wins))]
        V = []
        for k in range(d):
            X = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
            Q, _ = np.linalg.qr(X)
            rank = int(rng.integers(0, M + 1))
            P = Q[:, :rank] @ Q[:, :rank].conj().T
            V.append(1j ** (-a[k]) * (np.eye(M) - P) + 1j ** (-ap[k]) * P)
        F = sum((h[k - j] * np.trace(V[j] @ V[k].conj().T) / M).real for j in range(d) for k in range(j + 1, d))
        worst = max(worst, F - FD)
    ok2 &= worst < 1e-10
    print(f'   d = {d}: max F - F_DKZ over 60 random families = {worst:+.3e}')

print('ALL CHECKS PASSED' if ok1 and ok2 else 'CHECK FAILED')
print(f'# wall time: {time.time() - t0:.0f} s')
