"""
exact_tracial.py -- exact (cyclotomic) data for the symmetry-reduced level-2 tracial relaxation of the
maximally-entangled CGLMP problem, and exact DKZ evaluation / kernels.

Field: Q(zeta), zeta = exp(2 pi i/(4d)); w = zeta^4.
Generators R1..R4 -> g = 0..3, letters (g, n), n in 1..d-1.
Symmetries on moments (tracial functional L):
   cyclic:    L(uv) = L(vu)
   reversal:  L(w) = L(rev w)                    (strategy R -> R^T)
   rho:       L(w) = L(rho w), rho: R_g -> R_{g+1}, R_4 -> w R_1   (phase w^n for letters (3,n))
   sigma:     L(w) = L(sigma w), sigma: R_g -> R_{3-g}^{-1}
   hermitian: L(w)^* = L(w^dagger)
"""
import pickle
import sys
import itertools
from fractions import Fraction
from ncwords import reduce_word, dagger, cyclic_canon
from cyclo import Field


class Setup:
    def __init__(self, d, words):
        self.d = d
        self.N = 4 * d
        self.F = Field(self.N)
        self.words = words
        self.classes = {}      # canonical rep -> info
        self.cache = {}
        self.nvar = 0
        self.var_info = []     # var index -> (rep, 'real'/'re'/'im')

    # ---------- symmetry actions: return (exponent of zeta, word)
    def rho(self, w):
        e = 0
        out = []
        for (g, n) in w:
            if g == 3:
                e += 4 * n
                out.append((0, n))
            else:
                out.append((g + 1, n))
        return e % self.N, tuple(out)

    def sigma(self, w):
        d = self.d
        return 0, tuple((3 - g, (-n) % d) for (g, n) in w)

    def orbit(self, w):
        d, N = self.d, self.N
        start = cyclic_canon(w, d)
        seen = {start: 0}
        stack = [start]
        bad = False
        while stack:
            u = stack.pop()
            eu = seen[u]
            cands = [(0, cyclic_canon(tuple(reversed(u)), d))]
            e2, u2 = self.rho(u)
            cands.append((e2, cyclic_canon(u2, d)))
            e3, u3 = self.sigma(u)
            cands.append((e3, cyclic_canon(u3, d)))
            for (e, v) in cands:
                # L(u) = zeta^e L(v)  => L(w) = zeta^{eu} L(u) = zeta^{eu+e} L(v)
                ev = (eu + e) % N
                if v in seen:
                    if seen[v] != ev:
                        bad = True
                else:
                    seen[v] = ev
                    stack.append(v)
        return seen, bad

    def moment(self, w):
        """L(w) as dict var -> F coefficient ; var -1 = constant."""
        F, d, N = self.F, self.d, self.N
        w = reduce_word(w, d)
        key = cyclic_canon(w, d)
        if key in self.cache:
            return self.cache[key]
        orb, bad = self.orbit(key)
        kd = cyclic_canon(dagger(key, d), d)
        orb_d, bad_d = self.orbit(kd)
        if bad or bad_d:
            self.cache[key] = {}
            return {}
        rep = min(set(orb) | set(orb_d))
        if rep == ():
            if () in orb:
                res = {-1: F.z(orb[()])}
            else:
                res = {-1: F.z(-orb_d[()])}
            self.cache[key] = res
            return res
        cls = self.classes.get(rep)
        if cls is None:
            rd = cyclic_canon(dagger(rep, d), d)
            orb_rd, _ = self.orbit(rd)
            if rep in orb_rd:
                k = orb_rd[rep] % N          # L(rep)^* = L(rd) = zeta^k L(rep)
                if k == (2 * d) % N:
                    f = F.i_unit()
                else:
                    f = F.add(F.one, F.z(-k))
                v = self.nvar
                self.nvar += 1
                self.var_info.append((rep, 'real'))
                cls = ('real', v, f)
            else:
                v1, v2 = self.nvar, self.nvar + 1
                self.nvar += 2
                self.var_info.append((rep, 're'))
                self.var_info.append((rep, 'im'))
                cls = ('cplx', v1, v2)
            self.classes[rep] = cls
        # value of L(rep) as dict
        if cls[0] == 'real':
            Lrep = {cls[1]: cls[2]}
        else:
            Lrep = {cls[1]: F.one, cls[2]: F.i_unit()}
        if rep in orb:
            ph = F.z(orb[rep])
            res = {v: F.mul(ph, c) for v, c in Lrep.items()}
        else:
            ph = F.z(-orb_d[rep])
            res = {v: F.mul(ph, F.conj(c)) for v, c in Lrep.items()}
        self.cache[key] = res
        return res

    # ---------- symmetry-adapted SOS basis
    def blocks(self):
        F, d, N = self.F, self.d, self.N
        bysector = {}
        for u in self.words:
            q = sum(n for g, n in u) % d
            bysector.setdefault(q, []).append(u)
        out = []
        for q, ws in sorted(bysector.items()):
            orbits, done = [], set()
            for u in ws:
                if u in done:
                    continue
                orb = []
                e, cur = 0, u
                for k in range(4):
                    orb.append((e, cur))
                    done.add(cur)
                    e2, cur = self.rho(cur)
                    e = (e + e2) % N
                orbits.append(orb)
            for t in range(4):
                mu_exp = (q + d * t) % N          # mu = zeta^{q + d t}
                els = []
                for orb in orbits:
                    if orb[0][1] == ():
                        if mu_exp == 0:
                            els.append([(F.one, ())])
                        continue
                    el = [(F.z(e - k * mu_exp), wd) for k, (e, wd) in enumerate(orb)]
                    els.append(el)
                # drop structurally zero elements
                els = [el for el in els if self.inner(el, el)]
                if els:
                    out.append(((q, t), els))
        return out

    def inner(self, a, b):
        """<a|b> = L(a^dag b) as dict var -> coef"""
        F, d = self.F, self.d
        acc = {}
        for (ca, wa) in a:
            cac = F.conj(ca)
            for (cb, wb) in b:
                c0 = F.mul(cac, cb)
                for v, c in self.moment(reduce_word(dagger(wa, d) + wb, d)).items():
                    acc[v] = F.add(acc.get(v, F.zero), F.mul(c0, c))
        return {v: c for v, c in acc.items() if not F.is_zero(c)}

    def objective(self):
        """S = const + sum_v c_v y_v ; returns (dict var->F coef (real), const F)."""
        F, d = self.F, self.d
        acc = {-1: F.from_int(2 * (d - 1))}
        for n in range(1, d):
            # c_n = -1/(1 - w^{-n})
            cn = F.neg(F.inv(F.sub(F.one, F.z(-4 * n))))
            terms = [(cn, ((0, n), (1, (-n) % d))), (cn, ((1, n), (2, (-n) % d))),
                     (cn, ((2, n), (3, (-n) % d))), (F.mul(cn, F.z(-4 * n)), ((3, n), (0, (-n) % d)))]
            for (c, w) in terms:
                for v, cc in self.moment(w).items():
                    acc[v] = F.add(acc.get(v, F.zero), F.mul(c, cc))
        for v, c in acc.items():
            assert F.is_real(c), ("objective coefficient not real", v, c)
        const = acc.pop(-1)
        return {v: c for v, c in acc.items() if not F.is_zero(c)}, const


# ---------- exact DKZ strategy as generalised permutation matrices: M|k> = zeta^{e[k]} |p[k]>
def gp_mul(A, B, N):
    pA, eA = A
    pB, eB = B
    p = [pA[pB[k]] for k in range(len(pB))]
    e = [(eB[k] + eA[pB[k]]) % N for k in range(len(pB))]
    return (p, e)


def gp_pow(A, n, d, N):
    n %= d
    I = (list(range(d)), [0] * d)
    R = I
    for _ in range(n):
        R = gp_mul(A, R, N)
    return R


def gp_transpose(A, N):
    p, e = A
    d = len(p)
    pt = [0] * d
    et = [0] * d
    for k in range(d):
        pt[p[k]] = k
        et[p[k]] = e[k]
    return (pt, et)


def dkz(d):
    N = 4 * d
    R = []
    for g in range(4):
        a = 2 - g                      # R_{g+1} = G^a X^{-1} G^{-a}
        p = [(k - 1) % d for k in range(d)]
        e = [(-a) % N if k >= 1 else (a * (d - 1)) % N for k in range(d)]
        R.append((p, e))
    return R


def strategy_images(R, d):
    N = 4 * d
    def key(S):
        return tuple((tuple(p), tuple(e)) for (p, e) in S)
    def rho_s(S):
        p0, e0 = S[0]
        return [S[1], S[2], S[3], (p0, [(x + 4) % N for x in e0])]
    def sigma_s(S):
        inv = lambda A: gp_pow(A, d - 1, d, N)
        return [inv(S[3]), inv(S[2]), inv(S[1]), inv(S[0])]
    def rev_s(S):
        return [gp_transpose(A, N) for A in S]
    seen = {key(R): R}
    stack = [R]
    while stack:
        S = stack.pop()
        for T in (rho_s(S), sigma_s(S), rev_s(S)):
            k = key(T)
            if k not in seen:
                seen[k] = T
                stack.append(T)
    return list(seen.values())


def eval_element(el, S, d, F):
    """element = list of (coef, word) ; returns dict (row,col) -> F value of the d x d matrix."""
    N = 4 * d
    out = {}
    for (c, w) in el:
        M = (list(range(d)), [0] * d)
        for (g, n) in w:
            M = gp_mul(M, gp_pow(S[g], n, d, N), N)
        p, e = M
        for k in range(d):
            key = (p[k], k)
            out[key] = F.add(out.get(key, F.zero), F.mul(c, F.z(e[k])))
    return {k: v for k, v in out.items() if not F.is_zero(v)}


def nullspace(F, rows, ncols):
    """rows: list of dict col->F. returns list of basis vectors (list of F) of {x : row.x = 0}."""
    # RREF
    R = [dict(r) for r in rows if r]
    pivots = []
    prow = []
    for c in range(ncols):
        idx = None
        for i, r in enumerate(R):
            if c in r and not F.is_zero(r[c]):
                idx = i
                break
        if idx is None:
            continue
        r = R.pop(idx)
        inv = F.inv(r[c])
        r = {k: F.mul(v, inv) for k, v in r.items()}
        # eliminate from others
        newR = []
        for s in R:
            if c in s:
                f = s[c]
                t = dict(s)
                for k, v in r.items():
                    t[k] = F.sub(t.get(k, F.zero), F.mul(f, v))
                t = {k: v for k, v in t.items() if not F.is_zero(v)}
                if t:
                    newR.append(t)
            else:
                newR.append(s)
        R = newR
        for j in range(len(prow)):
            s = prow[j]
            if c in s:
                f = s[c]
                t = dict(s)
                for k, v in r.items():
                    t[k] = F.sub(t.get(k, F.zero), F.mul(f, v))
                prow[j] = {k: v for k, v in t.items() if not F.is_zero(v)}
        pivots.append(c)
        prow.append(r)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        x = [F.zero] * ncols
        x[fc] = F.one
        for pc, r in zip(pivots, prow):
            if fc in r:
                x[pc] = F.neg(r[fc])
        basis.append(x)
    return basis


def level2_words(d, filt=None):
    ws = [()]
    for i in range(4):
        for n in range(1, d):
            ws.append(((i, n),))
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            for n in range(1, d):
                for m in range(1, d):
                    u = ((i, n), (j, m))
                    if filt is None or filt(u, d):
                        ws.append(u)
    return ws


ADJ = lambda u, d: (u[1][0] - u[0][0]) % 4 in (1, 3)


def build_all(d, verbose=True, filt=None):
    import time
    t0 = time.time()
    st = Setup(d, level2_words(d, filt))
    F = st.F
    blocks = st.blocks()
    # entries
    bdata = []
    for (lab, els) in blocks:
        n = len(els)
        ent = {}
        for i in range(n):
            for j in range(i, n):
                ent[(i, j)] = st.inner(els[i], els[j])
        bdata.append((lab, els, ent))
    obj, const = st.objective()
    if verbose:
        print(f"d={d}: {len(blocks)} blocks, sizes {[len(e) for _, e in blocks]}, nvar={st.nvar} [{time.time()-t0:.1f}s]", flush=True)
    # exact DKZ value
    R = dkz(d)
    val = F.from_int(2 * (d - 1))
    tau = lambda el: F.scal(eval_trace(el, R, d, F), 1, d)
    for n in range(1, d):
        cn = F.neg(F.inv(F.sub(F.one, F.z(-4 * n))))
        terms = [(cn, ((0, n), (1, (-n) % d))), (cn, ((1, n), (2, (-n) % d))),
                 (cn, ((2, n), (3, (-n) % d))), (F.mul(cn, F.z(-4 * n)), ((3, n), (0, (-n) % d)))]
        for (c, w) in terms:
            val = F.add(val, F.mul(c, tau([(F.one, w)])))
    assert F.is_real(val)
    # kernels
    imgs = strategy_images(R, d)
    kernels = []
    for (lab, els, ent) in bdata:
        rows = []
        for S in imgs:
            mats = [eval_element(el, S, d, F) for el in els]
            keys = set()
            for m in mats:
                keys |= set(m.keys())
            for kk in keys:
                row = {j: mats[j][kk] for j in range(len(els)) if kk in mats[j]}
                rows.append(row)
        ker = nullspace(F, rows, len(els))
        kernels.append(ker)
    if verbose:
        print(f"   images={len(imgs)}  kernel dims {[len(k) for k in kernels]}  S_DKZ={F.to_complex(val).real} [{time.time()-t0:.1f}s]", flush=True)
    return dict(d=d, st_nvar=st.nvar, var_info=st.var_info, bdata=bdata, obj=obj, const=const,
                S_dkz=val, kernels=kernels, N=st.N)


def eval_trace(el, S, d, F):
    m = eval_element(el, S, d, F)
    tot = F.zero
    for (r, c), v in m.items():
        if r == c:
            tot = F.add(tot, v)
    return tot


if __name__ == "__main__":
    args = sys.argv[1:]
    filt = None
    if args and args[0] == "adj":
        filt = ADJ
        args = args[1:]
    for d in [int(t) for t in args]:
        data = build_all(d, filt=filt)
        data["wordset"] = "adjacent pairs" if filt else "full level 2"
        with open(f"exact_tracial_d{d}.pkl", "wb") as fh:
            pickle.dump(data, fh)
