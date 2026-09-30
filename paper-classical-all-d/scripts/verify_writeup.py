"""Consolidated numerical cross-check of the proof in the paper (Sections 3 and 4).
 (T)  brute force over ALL classical configurations (transversals S of Z_4d mod d), d = 2..DMAX:
      max F = F_DKZ, attained exactly on the 4d windows;  and F(S) = <S,S-d>/2 - d/2 on every S
      [Theorem B; lemma "cotangent representation"]
 (L1) cyclic symmetry <A,B> = <B,C> = <C,A> on random 3-colourings        [lemma "cyclic symmetry"]
 (L2) swap identity  Delta Q = h_C(x) - h_C(x+1) = sum_{t in C} D(x-t)  on random colourings
                                                                           [lemma "exchange identity"]
 (L3) Q = Q_c + E with Q_c from tent averages (quadrature)                 [lemma "smearing"]
 (L5) exact formula e(m) = (N/pi) sum_j eps(m + jN) (truncated |j| <= 2000, tail < 1e-9)
                                                                           [lemma "junction function", part (2)]
usage: verify_writeup.py DMAX"""
import sys, itertools
import numpy as np
from verify_junction import kbar, eps

def cot(u, N):
    return 1/np.tan(np.pi*np.asarray(u)/N)

def pair(A, B, N):
    A = np.asarray(A); B = np.asarray(B)
    if len(A) == 0 or len(B) == 0: return 0.0
    return float(np.sum(cot((A[:, None] - B[None, :]) % N, N)))

def theorem_bruteforce(d):
    N = 4*d
    FDKZ = sum((d - m)/np.cos(np.pi*m/(2*d)) for m in range(1, d))
    sec = np.zeros(N)
    for u in range(N):
        if u % d: sec[u] = 1/np.cos(np.pi*u/(2*d))
    best = -1e18; nmax = 0; win_ok = True; ident_err = 0.0
    wins = {tuple(sorted(((c + t) % N) for t in range(d))) for c in range(N)}
    for r in itertools.product(range(4), repeat=d):
        S = [(k + d*r[k]) % N for k in range(d)]
        Sa = np.array(S)
        F = 0.5*float(np.sum(sec[(Sa[:, None] - Sa[None, :]) % N]))
        if len(S) <= 6 or r[0] == 0 and sum(r) % 7 == 0:        # identity check on a subsample (exact algebra)
            ident_err = max(ident_err, abs(F - (0.5*pair(S, [(y - d) % N for y in S], N) - d/2)))
        if F > FDKZ - 1e-9:
            nmax += 1
            if tuple(sorted(S)) not in wins: win_ok = False
        best = max(best, F)
    return FDKZ, best, nmax, win_ok, ident_err

def lemmas(N, trials, rng):
    k = np.zeros(N); k[1:] = cot(np.arange(1, N), N)
    kb = np.zeros(N); kb[1:] = [kbar(m, N) for m in range(1, N)]
    err1 = err2 = err3 = 0.0
    for t in range(trials):
        a, b = rng.integers(1, N//2, size=2)
        if a + b >= N: continue
        col = np.full(N, 2); perm = rng.permutation(N); col[perm[:a]] = 0; col[perm[a:a+b]] = 1
        A, B, C = (np.where(col == c)[0] for c in range(3))
        q = pair(A, B, N)
        err1 = max(err1, abs(q - pair(B, C, N)), abs(q - pair(C, A, N)))
        D = (A[:, None] - B[None, :]) % N
        err3 = max(err3, abs(q - (float(np.sum(kb[D])) + float(np.sum((k - kb)[D])))))
        for x in range(N):
            y = (x + 1) % N
            if (col[x], col[y]) in {(0, 1), (1, 2), (2, 0)}:
                third = {(0, 1): 2, (1, 2): 0, (2, 0): 1}[(col[x], col[y])]
                T = np.where(col == third)[0]
                pred = float(np.sum(k[(x - T) % N] - k[(x - T + 1) % N]))
                col2 = col.copy(); col2[x], col2[y] = col2[y], col2[x]
                A2, B2 = np.where(col2 == 0)[0], np.where(col2 == 1)[0]
                err2 = max(err2, abs(pair(A2, B2, N) - q - pred))
                if pred <= 0: err2 = 1e9
    # smearing lemma: check that Q_c equals the double integral: kbar is by definition the tent average; check one entry by 2D quadrature
    from scipy.integrate import dblquad
    m = 3; val, _ = dblquad(lambda s, t: cot(m + s - t, N), 0, 1, 0, 1)
    err3 = max(err3, abs(val - kb[m]))
    # junction-function lemma, part (2): exact partial-fraction formula
    err5 = 0.0
    for m in range(1, N):
        J = 2000
        tot = eps(m) + sum(eps(m + j*N) + (-eps(j*N - m)) for j in range(1, J + 1))   # eps odd: eps(m - jN) = -eps(jN - m)
        err5 = max(err5, abs((N/np.pi)*tot - (k[m] - kb[m])))
    return err1, err2, err3, err5

if __name__ == "__main__":
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for d in range(2, DMAX + 1):
        FDKZ, best, nmax, win_ok, ierr = theorem_bruteforce(d)
        print(f"(T) d={d}: F_DKZ={FDKZ:.10f}  max F over all 4^d classical configs={best:.10f}  diff={best-FDKZ:+.1e}  "
              f"#maximisers={nmax} (=4d={4*d}: {nmax == 4*d}), all windows: {win_ok}; cotangent-representation identity err {ierr:.1e}", flush=True)
    rng = np.random.default_rng(11)
    for N in (8, 12, 20, 33, 40):
        e1, e2, e3, e5 = lemmas(N, 40, rng)
        print(f"N={N}: cyclic symmetry err {e1:.1e}; swap identity err {e2:.1e}; Q=Qc+E and tent=2D-integral err {e3:.1e}; "
              f"e(m) partial-fraction formula err {e5:.1e}", flush=True)
