"""local_hessian.py -- second-order analysis of the reduced problem at the twirled DKZ point (direct sum of the 4d
one-step classical configurations).  Claim: F^(2) = -(1/N) sum_{config pairs} sum_{j<k} w_jk |h_j - h_k|^2 with
w_jk = J(s_k-s_j) + J(s'_k-s'_j) - J(s'_k-s_j) - J(s_k-s'_j) >= 0.  Checks (i) the formula against finite
differences on random perturbations, (ii) w >= 0 for all pairs of orbit configurations."""
import itertools
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(5)


def J(m, t, d):
    return 2 * np.sin(np.pi * (m / d + t) / 2) / np.sin(np.pi * m / d)


def orbit(d):
    confs = []
    for a in range(1, d + 1):
        for c in range(4):
            confs.append(tuple((c + (1 if k >= a else 0)) % 4 for k in range(d)))
    return confs


def F_of(E, d):
    N = E[0].shape[0]
    F = 0.0
    for j in range(d):
        for k in range(j + 1, d):
            psi = np.pi * (k - j) / (2 * d)
            h = 1 / np.cos(psi) - 1j / np.sin(psi)
            F += (h * np.trace(E[j].conj().T @ E[k]) / N).real
    return F


for d in range(3, 9):
    confs = orbit(d)
    # (ii) weights
    minw = np.inf
    for s, sp in itertools.combinations(confs, 2):
        for j in range(d):
            for k in range(j + 1, d):
                m = k - j
                w = J(m, s[k] - s[j], d) + J(m, sp[k] - sp[j], d) - J(m, sp[k] - s[j], d) - J(m, s[k] - sp[j], d)
                dj = 1j ** s[j] - 1j ** sp[j]
                dk = 1j ** s[k] - 1j ** sp[k]
                if abs(dj) > 1e-12 and abs(dk) > 1e-12:
                    minw = min(minw, w)
    # (i) finite-difference check of F^(2) on the direct sum
    N = len(confs)
    E0 = [np.diag([1j ** c[k] for c in confs]) for k in range(d)]
    F0 = F_of(E0, d)
    H = []
    for k in range(d):
        A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
        H.append((A + A.conj().T) / 2)
    eps = 1e-3
    def Fe(e):
        return F_of([expm(1j * e * H[k]) @ E0[k] @ expm(-1j * e * H[k]) for k in range(d)], d)
    F2_fd = (Fe(eps) + Fe(-eps) - 2 * F0) / (2 * eps ** 2)
    F2 = 0.0
    for a_, b_ in itertools.combinations(range(N), 2):
        s, sp = confs[a_], confs[b_]
        hk = np.array([H[k][a_, b_] for k in range(d)])
        for j in range(d):
            for k in range(j + 1, d):
                m = k - j
                w = J(m, s[k] - s[j], d) + J(m, sp[k] - sp[j], d) - J(m, sp[k] - s[j], d) - J(m, s[k] - sp[j], d)
                dj = 1j ** s[j] - 1j ** sp[j]
                dk = 1j ** s[k] - 1j ** sp[k]
                # weight 2Re[kappa conj(dj) dk] equals w only if both nonzero; general formula:
                psi = np.pi * m / (2 * d)
                kap = 0.5 * (1 / np.cos(psi) - 1j / np.sin(psi))
                wt = 2 * (kap * np.conj(dj) * dk).real
                F2 -= wt * abs(hk[j] - hk[k]) ** 2
    F2 /= N
    # F(eps) = F0 + eps^2 F2 + ... -> second derivative/2 = F2
    print(f"d={d}: F0={F0:.10f} (F_DKZ={np.sum((d-np.arange(1,d))/np.cos(np.pi*np.arange(1,d)/(2*d))):.10f}); "
          f"min weight over config pairs (on supp) = {minw:+.4f};  F2 formula {F2:.6f} vs finite diff {F2_fd:.6f}")
