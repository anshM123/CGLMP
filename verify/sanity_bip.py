"""
sanity_bip.py -- end-to-end numerical sanity check of an all-states certificate cert_bip_d{d}.pkl on RANDOM
bipartite strategies: average the moment functional over the images of the strategy under rho, sigma
(party swap) and complex conjugation; compare S - lam with sum_ij X_ij L(e_j^+ e_i).
"""
import sys
import pickle
import numpy as np
import mpmath as mp
from cyclo import Field
from extfield import Ext


def rand_u(d, D, rng):
    Z = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    Q, _ = np.linalg.qr(Z)
    w = np.exp(2j * np.pi * rng.integers(0, d, size=D) / d)
    return (Q * w) @ Q.conj().T


def images(ops, psi, d):
    """ops = [A1,B1,A2,B2] (A on Alice space, B on Bob space); psi as DA x DB matrix."""
    w = np.exp(2j * np.pi / d)
    inv = lambda M: M.conj().T
    out = []
    stack = [(tuple(ops), psi)]
    seen = []
    while stack:
        T, ps = stack.pop()
        if any(all(np.allclose(a, b) for a, b in zip(T, U)) and np.allclose(ps, p2) for (U, p2) in seen):
            continue
        seen.append((T, ps))
        # rho: A1'=B1, B1'=A2, A2'=B2, B2'=w A1 ; parties swapped -> psi transposed
        stack.append(((T[1], T[2], T[3], w * T[0]), ps.T))
        stack.append(((inv(T[3]), inv(T[2]), inv(T[1]), inv(T[0])), ps.T))
        stack.append((tuple(M.T for M in T), ps.conj()))
        if len(seen) > 300:
            break
    return seen


def moment(T, ps, a, b, d):
    MA = np.eye(ps.shape[0], dtype=complex)
    for (g, n) in a:
        MA = MA @ np.linalg.matrix_power(T[g], n % d)
    MB = np.eye(ps.shape[1], dtype=complex)
    for (g, n) in b:
        MB = MB @ np.linalg.matrix_power(T[g], n % d)
    # <psi| MA (x) MB |psi> with psi as matrix: tr(psi^dag MA psi MB^T)
    return np.trace(ps.conj().T @ MA @ ps @ MB.T)


if __name__ == "__main__":
    d = int(sys.argv[1])
    D = int(sys.argv[2]) if len(sys.argv) > 2 else d
    cert = pickle.load(open(f"cert_bip_d{d}.pkl", "rb"))
    F0 = Field(4 * d)
    mp.mp.dps = 40
    L = Ext(F0, cert["p"], mp.mpf(cert["lamstr"]))
    lam = L.to_cfloat(cert["lam"]).real
    w = np.exp(2j * np.pi / d)
    rng = np.random.default_rng(9)
    for trial in range(2):
        ops = [rand_u(d, D, rng) for _ in range(4)]
        ps = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
        ps /= np.linalg.norm(ps)
        imgs = images(ops, ps, d)
        def Lbar(a, b):
            return np.mean([moment(T, p, a, b, d) for (T, p) in imgs])
        S = 2 * (d - 1)
        for n in range(1, d):
            c = -1 / (1 - w ** (-n))
            S += c * (Lbar(((0, n),), ((1, -n % d),)) + Lbar(((2, -n % d),), ((1, n),))
                      + Lbar(((2, n),), ((3, -n % d),)) + w ** (-n) * Lbar(((0, -n % d),), ((3, n),)))
        tot = 0.0
        for bidx in cert["keep"]:
            lab, els, _ = cert["bdata"][bidx]
            N = np.array([[L.to_cfloat(x) for x in col] for col in cert["kernels"][bidx]]).T
            Y = np.array([[L.to_cfloat(x) for x in row] for row in cert["Y"][bidx]])
            X = N @ Y @ N.conj().T
            els_c = [[(F0.to_cfloat(c), wd) for (c, wd) in el] for el in els]
            n = len(els)
            G = np.zeros((n, n), complex)
            for j in range(n):
                for i in range(n):
                    val = 0
                    for (cj, (aj, bj)) in els_c[j]:
                        for (ci, (ai, bi)) in els_c[i]:
                            adj = lambda wd: tuple((g, (-k) % d) for (g, k) in reversed(wd))
                            val += np.conj(cj) * ci * Lbar(adj(aj) + ai, adj(bj) + bi)
                    G[j, i] = val
            tot += np.real(np.sum(X * G.T))
        print(f"d={d} D={D} trial {trial}: #images {len(imgs)}  S - lam = {S.real - lam:.10f}  SOS = {tot:.10f}  diff = {S.real - lam - tot:.2e}", flush=True)
