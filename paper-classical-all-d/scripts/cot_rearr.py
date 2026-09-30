"""Discrete cot-rearrangement test: for disjoint A, B in Z_N with |A| = a, |B| = b, is
<A,B> := sum_{x in A, y in B} cot(pi (x-y)/N) maximised by adjacent intervals B = [-b,-1], A = [0,a-1] ?
Also check the classical identity F(S) = (1/2)<S, S-d> - d/2 for transversals S of Z_4d mod d."""
import sys, itertools
import numpy as np

def pair(A, B, N):
    A = np.array(A); B = np.array(B)
    D = (A[:, None] - B[None, :]) % N
    return np.sum(1/np.tan(np.pi*D/N))

def check_identity(d):
    N = 4*d
    rng = np.random.default_rng(0)
    for t in range(5):
        r = rng.integers(0, 4, size=d)
        S = [(k + d*r[k]) % N for k in range(d)]
        Fdir = 0.5*sum(1/np.cos(np.pi*(x-y)/(2*d)) for x in S for y in S if x != y)
        Fcot = 0.5*pair(S, [(y - d) % N for y in S], N) - d/2
        assert abs(Fdir - Fcot) < 1e-8, (Fdir, Fcot)
    return True

def brute(N, a, b):
    best = -1e18; arg = None
    adj = pair(list(range(a)), [(-j) % N for j in range(1, b+1)], N)
    for B in itertools.combinations(range(1, N), b-1):
        B = (0,) + B
        rest = [x for x in range(N) if x not in B]
        # for fixed B, the best A = top-a values of h_B on the complement (bathtub): exact max over A
        h = np.array([sum(1/np.tan(np.pi*((x-y) % N)/N) for y in B) for x in rest])
        val = np.sort(h)[::-1][:a].sum()
        if val > best + 1e-9: best = val; arg = (B, [rest[i] for i in np.argsort(-h)[:a]])
    return adj, best, arg

if __name__ == "__main__":
    for d in [3, 4, 5, 6]:
        check_identity(d)
    print("identity F(S) = <S,S-d>/2 - d/2 verified on random transversals")
    for spec in sys.argv[1:]:
        N, a, b = map(int, spec.split(','))
        adj, best, arg = brute(N, a, b)
        print(f"N={N} a={a} b={b}: adjacent {adj:.6f}  max over (A,B) {best:.6f}  diff {best-adj:+.2e}  argmax B={arg[0]} A={sorted(arg[1])}", flush=True)
