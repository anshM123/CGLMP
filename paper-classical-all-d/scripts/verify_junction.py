"""Verification of the junction inequality ingredients (lemmas "junction function" and "junction inequality").
e(m) = cot(pi m/N) - kbar(m),  kbar(m) = int_{-1}^{1} (1-|v|) cot(pi(m+v)/N) dv.
Claims: (i) e(m) <= 0 on [1, N/2], (N/pi) eps(m) <= e(m) with eps(1) = 1 - 2 log 2, eps(u) = 1/u - D2(x log x)(u) (u>=2);
(ii) |e(1)| >= (N/pi)(|eps(1)| - |eps(N-1)|);  (iii) sc = sum_{m=2}^{N/2} (m-1)|e(m)| <= 0.07785 N/pi,
fw = sum_{m=2}^{N/2} m|e(m)| <= 0.11481 N/pi;  (iv) the chain Q = Q_c + E, E(A,B) - E_block <= -(w-1)|e1| + w sc + fw < 0."""
import sys
import numpy as np
from scipy.integrate import quad
def kbar(m, N):
    f = lambda v: (1 - abs(v)) / np.tan(np.pi * (m + v) / N)
    pts = [p for p in (-m, N - m) if -1 < p < 1]
    if pts:   # singular point inside: split and use the log-integrable behaviour
        p = pts[0]
        a1, _ = quad(f, -1, p, limit=400); a2, _ = quad(f, p, 1, limit=400)
        return a1 + a2
    val, _ = quad(f, -1, 1, limit=200)
    return val
def eps(u):
    if u == 1: return 1 - 2*np.log(2)
    x = float(u)
    return 1/x - ((x+1)*np.log(x+1) - 2*x*np.log(x) + (x-1)*np.log(x-1))
if __name__ == "__main__":
  for N in [int(a) for a in sys.argv[1:]]:
      e = {m: 1/np.tan(np.pi*m/N) - kbar(m, N) for m in range(1, N)}
      ok_sign = all(e[m] <= 1e-9 for m in range(1, N//2 + 1))
      ok_low = all(e[m] >= (N/np.pi)*eps(m) - 1e-7 for m in range(1, N//2 + 1))
      e1 = -e[1]; e1_low = (N/np.pi)*(abs(eps(1)) - abs(eps(N-1)))
      sc = sum((m-1)*abs(e[m]) for m in range(2, N//2 + 1)); fw = sum(m*abs(e[m]) for m in range(2, N//2 + 1))
      c = N/np.pi
      print(f"N={N}: sign ok {ok_sign}, lower ok {ok_low}; |e1|/(N/pi)={e1/c:.5f} (bound {e1_low/c:.5f}); sc/(N/pi)={sc/c:.5f} (<=0.07785); fw/(N/pi)={fw/c:.5f} (<=0.11481); w=2 margin {( -e1 + 2*sc + fw)/c:+.4f}", flush=True)
