"""classical_z4.py -- brute force of the classical (commuting) reduced problem: s in Z_4^d,
F(s) = sum_{j<k} J_{k-j}(s_k - s_j),  J_m(t) = 2 sin(pi(m/d + t)/2)/sin(pi m/d)."""
import itertools, numpy as np, sys
for d in range(2, 11):
    m = np.arange(1, d)
    Fdkz = np.sum((d - m) / np.cos(np.pi * m / (2 * d)))
    J = np.zeros((d, 4))
    for mm in range(1, d):
        for t in range(4):
            J[mm, t] = 2 * np.sin(np.pi * (mm / d + t) / 2) / np.sin(np.pi * mm / d)
    vals = {}
    best = []
    for s in itertools.product(range(4), repeat=d - 1):
        s = (0,) + s
        F = 0.0
        for j in range(d):
            for k in range(j + 1, d):
                F += J[k - j, (s[k] - s[j]) % 4]
        best.append((F, s))
    best.sort(reverse=True)
    top = [b for b in best if b[0] > Fdkz - 1e-9]
    nxt = [b for b in best if b[0] <= Fdkz - 1e-9][:3]
    print(f"d={d}: F_DKZ={Fdkz:.6f} max={best[0][0]:.6f} #optimal={len(top)}: {[b[1] for b in top][:12]}  next: {[(round(b[0],4), b[1]) for b in nxt]}")
