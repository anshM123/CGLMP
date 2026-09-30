"""full_check_r5.py -- independent (non-reduced) check of the R5 KL certificate on |Phi_5>:
full correlation table q(a,b|x,y) from explicit R5 bases (numpy), test factor r(a,b|x,y) rational,
exact feasibility over ALL 5^4 deterministic strategies, lower bound in interval arithmetic with q
recomputed from the exact-phase formula."""
import itertools
from fractions import Fraction as Fr
import numpy as np
from mpmath import iv, mpf
from bell22d import probs_pure, maxent, kl_strength
from state_opt import r5_bases
from kl_certificate import r5_phases_q, P_interval

d = 5
A, B = r5_bases()
q = probs_pure(maxent(d), A, B)
res = kl_strength(q)
R = {}
for x, y, a, b in itertools.product(range(2), range(2), range(d), range(d)):
    R[(x, y, a, b)] = Fr(float(q[x, y, a, b] / res['p'][x, y, a, b])).limit_denominator(10 ** 12)
Mx = max(sum(R[(x, y, lam[x], lam[2 + y])] for x in range(2) for y in range(2)) / 4
         for lam in itertools.product(range(d), repeat=4))
th, et = r5_phases_q()
P = P_interval(th, et, d)          # link distributions P_xy(m), m = b - a
# consistency of the numpy table with the exact formula
dev = max(abs(q[x, y, a, b] - float(mpf(P[x][y][(b - a) % d].mid)) / d)
          for x, y, a, b in itertools.product(range(2), range(2), range(d), range(d)))
L = iv.mpf(0)
for (x, y, a, b), rr in R.items():
    rn = rr / Mx * (1 - Fr(1, 10 ** 9))
    L += (P[x][y][(b - a) % d] / d) * iv.log(iv.mpf(rn.numerator) / rn.denominator) / 4
L = L / iv.log(iv.mpf(2))
print(f"R5 full check: max |q_numpy - q_exact| = {dev:.1e}; exact normaliser over all 625 strategies = {float(Mx):.15f}; "
      f"certified KL >= {float(mpf(L.a)):.12f} bits (numerical {res['kl']/np.log(2):.12f})")
