"""sym_check.py -- DKZ correlations on |Phi_d> are invariant under
  G1: q'(a,b|x,y) = q(-a, -b + c_y | 1-x, y),  c = (0, 1)
  G2: q'(a,b|x,y) = q(-a + c'_x, -b | x, 1-y), c' = (-1, 0)
  Zd: q'(a,b|x,y) = q(a+s, b+s | x, y)
(closed form p = 1/(2 d^3 sin^2(pi (b-a-delta_xy)/d)), delta = (1/4, 3/4; -1/4, 1/4)).  G1, G2 act on
setting pairs as (x,y) -> (1-x,y), (x,1-y): transitive  =>  uniform settings optimal (Lemma S)."""
import numpy as np
from bell22d import dkz_bases, maxent, probs_pure

for d in range(2, 10):
    A, B = dkz_bases(d)
    q = probs_pure(maxent(d), A, B)
    a = np.arange(d)
    e1 = e2 = e3 = e4 = 0.0
    for x in range(2):
        for y in range(2):
            c = (0, 1)[y]
            q1 = q[1 - x, y][np.ix_((-a) % d, (-a + c) % d)]
            e1 = max(e1, np.abs(q1 - q[x, y]).max())
            cp = (-1, 0)[x]
            q2 = q[x, 1 - y][np.ix_((-a + cp) % d, (-a) % d)]
            e2 = max(e2, np.abs(q2 - q[x, y]).max())
            q3 = q[x, y][np.ix_((a + 1) % d, (a + 1) % d)]
            e3 = max(e3, np.abs(q3 - q[x, y]).max())
            m = (a[None, :] - a[:, None]) % d
            delta = {(0, 0): 0.25, (0, 1): 0.75, (1, 0): -0.25, (1, 1): 0.25}[(x, y)]
            closed = 1 / (2 * d ** 3 * np.sin(np.pi * (m - delta) / d) ** 2)
            e4 = max(e4, np.abs(closed - q[x, y]).max())
    print(f"d={d}: max deviation  G1 {e1:.1e}  G2 {e2:.1e}  Z_d {e3:.1e}  closed form {e4:.1e}")
