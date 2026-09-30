"""
modsolve3.py -- parallel version of modsolve2 (same mathematics): the per-prime work (evaluation under all
embeddings, vectorised Gauss-Jordan in A_j = F_P[t]/(p_j), Vandermonde reconstruction) runs in worker
processes; the parent does CRT, rational reconstruction and the caller's exact check.
"""
import os
import math
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import modsolve2 as M2
from modsolve import is_prime, mat_inv_mod, try_recon, to_elements

_G = {}


def _init(nums, dens, m, k, n, N, p_ext):
    _G.update(nums=nums, dens=dens, m=m, k=k, n=n, N=N, p_ext=p_ext)


def _residues(P):
    nums, dens = _G["nums"], _G["dens"]
    if isinstance(nums, np.ndarray):
        return nums % P, dens % P
    num = np.array([[ai % P for ai in a] for a in nums], dtype=np.int64)
    den = np.array([dn % P for dn in dens], dtype=np.int64)
    return num, den


def _one_prime(args):
    P, pivots = args
    g = _G
    m, k, n, N, p_ext = g["m"], g["k"], g["n"], g["N"], g["p_ext"]
    try:
        emb = M2.Emb(N, n, p_ext, P)
        num, den = _residues(P)
        if np.any(den == 0):
            return None
        inv = M2._modinv_vec(den, P)
        vals = np.zeros((len(den), n), dtype=np.int64)
        for i in range(n):
            vals = (vals + (num[:, i:i + 1] * emb.T[None, :, i]) % P) % P
        vals = (vals * inv[:, None]) % P
        A = vals.reshape(m, m + 1, k, n).transpose(3, 0, 1, 2).copy()
    except (ValueError, ZeroDivisionError):
        return None
    pc = emb.pc if k > 1 else None
    order = []
    rows_used = np.zeros(m, dtype=bool)
    try:
        cols_iter = range(m) if pivots is None else [c for (_, c) in pivots]
        pivmap = None if pivots is None else {c: r for (r, c) in pivots}
        for c in cols_iter:
            if pivots is None:
                r = None
                for rr in range(m):
                    if not rows_used[rr] and np.all(np.any(A[:, rr, c, :] % P != 0, axis=-1)):
                        r = rr
                        break
                if r is None:
                    continue
            else:
                r = pivmap[c]
            rows_used[r] = True
            order.append((r, c))
            pv = A[:, r, c, :]
            if k == 1:
                if np.any(pv[:, 0] % P == 0):
                    return None
                ivv = M2._modinv_vec(pv[:, 0], P)[:, None]
                A[:, r, :, :] = (A[:, r, :, :] * ivv[:, None, :]) % P
                f = A[:, :, c, :].copy()
                f[:, r, :] = 0
                A = (A - (f[:, :, None, :] * A[:, r, None, :, :]) % P) % P
            else:
                ivv = M2.ainv_scalar(pv, pc, P)
                rowr = A[:, r, :, :]
                A[:, r, :, :] = M2.amul(rowr.transpose(1, 0, 2), np.broadcast_to(ivv, rowr.transpose(1, 0, 2).shape), pc, P).transpose(1, 0, 2)
                f = A[:, :, c, :].copy()
                f[:, r, :] = 0
                prod = M2.amul(np.broadcast_to(f[:, :, None, :], A.shape).transpose(1, 2, 0, 3),
                               np.broadcast_to(A[:, r, None, :, :], A.shape).transpose(1, 2, 0, 3), pc, P).transpose(2, 0, 1, 3)
                A = (A - prod) % P
    except ZeroDivisionError:
        return None
    piv = order if pivots is None else pivots
    Z = np.zeros((n, m, k), dtype=np.int64)
    for (r, c) in piv:
        Z[:, c, :] = A[:, r, m, :]
    V = [[pow(emb.rj[j], i, P) for i in range(n)] for j in range(n)]
    Vinv = np.array(mat_inv_mod(V, P), dtype=object)
    C = (Vinv.dot(Z.astype(object).reshape(n, m * k))) % P
    C2 = np.empty((n * k, m), dtype=object)
    for i in range(n):
        for c in range(m):
            for b in range(k):
                C2[i * k + b, c] = C[i, c * k + b]
    return P, order, C2


def solve(Mat, rhs, N, n_phi, p_ext=None, verbose=False, max_primes=6000, check=None, workers=None, batch=None):
    m = len(Mat)
    k = 1 if p_ext is None else len(p_ext) - 1
    n = n_phi
    rows_list = [list(Mat[i]) + [rhs[i]] for i in range(m)]
    nums, dens = M2._flatten(rows_list, k)
    # compact int64 packing when all numerators/denominators fit (saves worker memory)
    lim = 1 << 62
    if all(abs(v) < lim for a in nums for v in a) and all(dn < lim for dn in dens):
        nums = np.array(nums, dtype=np.int64)
        dens = np.array(dens, dtype=np.int64)
    workers = workers or max(1, min(8, (os.cpu_count() or 2) - 4))
    batch = batch or workers
    # prime generator
    def primes():
        P = (1 << 31) - 1
        P = P - (P % N) + 1
        while P >= (1 << 31):
            P -= N
        while True:
            P -= N
            if is_prime(P):
                yield P
    gen = primes()
    pivots = None
    acc, M_acc, nprimes = None, 1, 0
    with ProcessPoolExecutor(max_workers=workers, initializer=_init,
                             initargs=(nums, dens, m, k, n, N, p_ext)) as ex:
        # first prime serially (in a worker) to fix the pivot pattern
        while pivots is None:
            res = ex.submit(_one_prime, (next(gen), None)).result()
            if res is not None:
                P, order, C2 = res
                pivots = order
                acc, M_acc, nprimes = C2, P, 1
                if verbose:
                    print(f"   modsolve3: rank {len(pivots)} of {m}; {workers} workers", flush=True)
        while nprimes < max_primes:
            tasks = [(next(gen), pivots) for _ in range(batch)]
            for res in ex.map(_one_prime, tasks):
                if res is None:
                    continue
                P, order, C2 = res
                invM = pow(M_acc % P, P - 2, P)
                diff = ((C2 - acc) % P * invM) % P
                acc = acc + M_acc * diff
                M_acc *= P
                nprimes += 1
            if verbose and nprimes % 50 < batch:
                print(f"   modsolve3: {nprimes} primes, {M_acc.bit_length()} bits", flush=True)
            rec = try_recon(acc, M_acc)
            if rec is not None:
                sol = to_elements(rec, m, n, k)
                if check is None or check(sol):
                    if verbose:
                        print(f"   modsolve3: success with {nprimes} primes ({M_acc.bit_length()} bits)", flush=True)
                    return sol
    raise RuntimeError("modsolve3: too many primes")
