"""Full chain on random configurations:
  (a) random 3-colourings -> swap flow -> monotone; Q increases at each swap by exactly sum_{t in C} D(p-t).
  (b) for the final monotone colouring: Q_c <= Q_c,block (continuous Stein-Weiss), E <= E_block (junction inequality),
      hence Q <= Q_block; report margins.  Q_c computed exactly from kbar (tent averages)."""
import sys
import numpy as np
from verify_junction import kbar
rng = np.random.default_rng(7)
def run(N, a, b, trials):
    k = np.zeros(N); k[1:] = 1/np.tan(np.pi*np.arange(1, N)/N)
    kb = np.zeros(N); kb[1:] = [kbar(m, N) for m in range(1, N)]
    def pair(A, B, ker):
        D = (A[:, None] - B[None, :]) % N
        return float(np.sum(ker[D]))
    A0 = np.arange(a); B0 = (np.arange(-b, 0)) % N
    Qb, Qcb = pair(A0, B0, k), pair(A0, B0, kb); Eb = Qb - Qcb
    worst = {'swapid': 0, 'Qc': -1e9, 'E': -1e9, 'Q': -1e9}
    for t in range(trials):
        col = np.full(N, 2); perm = rng.permutation(N); col[perm[:a]] = 0; col[perm[a:a+b]] = 1
        while True:
            oo = [x for x in range(N) if (col[x], col[(x+1) % N]) in {(0,1),(1,2),(2,0)}]
            if not oo: break
            x = oo[rng.integers(len(oo))]; y = (x+1) % N
            A = np.where(col == 0)[0]; B = np.where(col == 1)[0]; Q0 = pair(A, B, k)
            third = {(0,1): 2, (1,2): 0, (2,0): 1}[(col[x], col[y])]
            T = np.where(col == third)[0]
            pred = float(np.sum(k[(x - T) % N] - k[(x - T + 1) % N]))
            col[x], col[y] = col[y], col[x]
            A = np.where(col == 0)[0]; B = np.where(col == 1)[0]; Q1 = pair(A, B, k)
            worst['swapid'] = max(worst['swapid'], abs((Q1 - Q0) - pred))
        A = np.where(col == 0)[0]; B = np.where(col == 1)[0]
        Q, Qc = pair(A, B, k), pair(A, B, kb); E = Q - Qc
        worst['Qc'] = max(worst['Qc'], Qc - Qcb); worst['E'] = max(worst['E'], E - Eb); worst['Q'] = max(worst['Q'], Q - Qb)
    c = N/np.pi
    print(f"N={N} a={a} b={b}: swap identity max err {worst['swapid']:.1e}; max(Qc-Qc_block)={worst['Qc']:+.2e}; "
          f"max(E-E_block)/(N/pi)={worst['E']/c:+.4f}; max(Q-Q_block)={worst['Q']:+.3e}", flush=True)
for spec in sys.argv[1:]:
    N, a, b, tr = map(int, spec.split(','))
    run(N, a, b, tr)
