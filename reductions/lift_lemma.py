"""lift_lemma.py -- classical 'telescoping lift' lemma check.
For s in Z_4^d: F(s) = sum_{j in Z_d} sum_{m=1}^{d-1} csc(pi m/2d) chi(n_{j+m} - n_j) (chi = [.=1]-[.=3] mod 4), and for ANY lift
n (n_{k+d} = n_k + 1) sum_j (n_{j+m}-n_j) = m.  If a lift has all D_m(j) >= -1 then F <= F_DKZ.  Count such configs and verify."""
import itertools, sys
import numpy as np


def F_arc(s, d):
    tot = 0.0
    for j in range(d):
        for m in range(1, d):
            k = j + m
            t = (s[k % d] + (k // d) - s[j]) % 4
            tot += ((t == 1) - (t == 3)) / np.sin(np.pi * m / (2 * d))
    return tot


for d in range(3, int(sys.argv[1]) + 1):
    Fd = sum(m / np.sin(np.pi * m / (2 * d)) for m in range(1, d))
    ok = 0; tot = 0; viol = 0; max_nonlift = -np.inf
    for s in itertools.product(range(4), repeat=d - 1):
        s = (0,) + s
        tot += 1
        F = F_arc(s, d)
        # try to find an almost-monotone lift: pairwise delta in [-1,2] for j<k, consistency
        n = [0]
        good = True
        for k in range(1, d):
            # n_k must satisfy n_k - n_j in [-1,2] for all j<k and n_k = s_k mod 4
            lo = max(n[j] - 1 for j in range(k)); hi = min(n[j] + 2 for j in range(k))
            cand = [x for x in range(lo, hi + 1) if (x - s[k]) % 4 == 0]
            if not cand:
                good = False; break
            n.append(cand[0])
        if good:
            ok += 1
            if F > Fd + 1e-9:
                viol += 1
        else:
            max_nonlift = max(max_nonlift, F)
    print(f"d={d}: configs {tot}, almost-monotone liftable {ok} ({100*ok/tot:.1f}%), violations {viol}; "
          f"max F among non-liftable = {max_nonlift:.4f} (F_DKZ = {Fd:.4f}, and F_DKZ-2csc/.. gap)")
