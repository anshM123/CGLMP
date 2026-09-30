"""adgl_check.py -- sanity check of B-T3 on the ADGL optimum (DKZ measurements + Perron state)."""
import sys
import numpy as np
sys.path.insert(0, r"<project>/iqoqi\programs\oqp27B")
from dkz_tracial import unitaries
from twirl_allstates import S_beta
for d in (3, 4, 5, 6):
    R = unitaries(d)
    k = np.arange(d)
    K = 1 / np.cos(np.pi * (k[:, None] - k[None, :]) / (2 * d))
    ev, U = np.linalg.eigh(K)
    v = np.abs(U[:, -1])
    best = None
    # the Schmidt basis of the DKZ convention may be permuted/reflected; try v, reversed v
    for vv in (v, v[::-1]):
        L = np.diag(vv ** 2)
        S = S_beta(R, L, d).real
        I = 4 - 2 * S / (d - 1)
        best = I if best is None else max(best, I)
    print(f"d={d}: I(DKZ, Perron state) = {best:.10f}   2(lmax-1)/(d-1) = {2*(ev[-1]-1)/(d-1):.10f}")
