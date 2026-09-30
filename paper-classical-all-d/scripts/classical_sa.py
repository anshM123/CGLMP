"""classical_sa.py -- Sherali-Adams-type LP relaxations of the CLASSICAL reduced problem (Z_4 clock model, s in Z_4^d):
variables: for each triple j<k<l a distribution q_jkl(a,b) on (s_k - s_j, s_l - s_k) in Z_4^2; pair marginals
(differences) must agree across triples.  Objective F = sum_{j<k} sum_t J_{k-j}(t) pi_jk(t).  Level 2 (pairs only) and
level 3 (triples).  Solved with scipy HiGHS."""
import sys, itertools
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog


def J(m, t, d):
    return 2 * np.sin(np.pi * (m / d + t) / 2) / np.sin(np.pi * m / d)


def solve(d, level=3):
    pairs = [(j, k) for j in range(d) for k in range(j + 1, d)]
    pidx = {p: i for i, p in enumerate(pairs)}
    npv = 4 * len(pairs)          # pi_jk(t)
    triples = list(itertools.combinations(range(d), 3)) if level >= 3 else []
    ntv = 16 * len(triples)
    nv = npv + ntv
    rows, cols, vals, rhs = [], [], [], []
    r = 0
    # normalisation of pairs
    for p in pairs:
        for t in range(4):
            rows.append(r); cols.append(4 * pidx[p] + t); vals.append(1.0)
        rhs.append(1.0); r += 1
    # triple consistency: marginals
    for ti, (j, k, l) in enumerate(triples):
        base = npv + 16 * ti
        # (j,k) marginal: sum_b q(a,b) = pi_jk(a)
        for a in range(4):
            for b in range(4):
                rows.append(r); cols.append(base + 4 * a + b); vals.append(1.0)
            rows.append(r); cols.append(4 * pidx[(j, k)] + a); vals.append(-1.0)
            rhs.append(0.0); r += 1
        # (k,l) marginal
        for b in range(4):
            for a in range(4):
                rows.append(r); cols.append(base + 4 * a + b); vals.append(1.0)
            rows.append(r); cols.append(4 * pidx[(k, l)] + b); vals.append(-1.0)
            rhs.append(0.0); r += 1
        # (j,l) marginal: a+b = c
        for c in range(4):
            for a in range(4):
                b = (c - a) % 4
                rows.append(r); cols.append(base + 4 * a + b); vals.append(1.0)
            rows.append(r); cols.append(4 * pidx[(j, l)] + c); vals.append(-1.0)
            rhs.append(0.0); r += 1
    A = sp.csr_matrix((vals, (rows, cols)), shape=(r, nv))
    c = np.zeros(nv)
    for (j, k) in pairs:
        for t in range(4):
            c[4 * pidx[(j, k)] + t] = -J(k - j, t, d)
    res = linprog(c, A_eq=A, b_eq=np.array(rhs), bounds=(0, None), method="highs")
    return -res.fun, res


if __name__ == "__main__":
    for d in range(3, int(sys.argv[1]) + 1):
        m = np.arange(1, d)
        Fd = np.sum((d - m) / np.cos(np.pi * m / (2 * d)))
        v2, _ = solve(d, 2)
        v3, _ = solve(d, 3)
        print(f"d={d}: F_DKZ={Fd:.8f}  SA-pairs={v2:.6f}  SA-triples={v3:.8f}  gap3={v3 - Fd:+.2e}", flush=True)
