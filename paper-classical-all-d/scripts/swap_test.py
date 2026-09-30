"""3-colorings of Z_N (A,B,C with counts a,b,c).  Q = <A,B> = sum_{x in A, y in B} cot(pi(x-y)/N).
(1) Does swapping every out-of-order adjacent pair (A,B)->(B,A), (B,C)->(C,B), (C,A)->(A,C) increase Q?
(2) In-order (monotone) multi-winding colorings: Q vs block."""
import sys, itertools
import numpy as np

def Q(col, N):
    A = [i for i in range(N) if col[i] == 0]; B = [i for i in range(N) if col[i] == 1]
    if not A or not B: return 0.0
    D = (np.array(A)[:, None] - np.array(B)[None, :]) % N
    return float(np.sum(1/np.tan(np.pi*D/N)))

out_of_order = {(0, 1), (1, 2), (2, 0)}   # (A,B), (B,C), (C,A) with A=0,B=1,C=2
for spec in sys.argv[1:]:
    N, a, b = map(int, spec.split(',')); c = N - a - b
    block = [2]*c + [1]*b + [0]*a
    Qb = Q(block, N)
    nbad = 0; ntot = 0; maxQ = -1e9
    for pos in itertools.combinations(range(N), a):
        rest = [i for i in range(N) if i not in pos]
        for posB in itertools.combinations(rest, b):
            col = [2]*N
            for i in pos: col[i] = 0
            for i in posB: col[i] = 1
            q = Q(col, N); maxQ = max(maxQ, q)
            for x in range(N):
                pr = (col[x], col[(x+1) % N])
                if pr in out_of_order:
                    ntot += 1
                    col2 = col.copy(); col2[x], col2[(x+1) % N] = col2[(x+1) % N], col2[x]
                    if Q(col2, N) < q - 1e-9: nbad += 1
    print(f"N={N} a={a} b={b}: Q_block={Qb:.4f} maxQ={maxQ:.4f}; out-of-order swaps tested {ntot}, decreasing: {nbad}", flush=True)
