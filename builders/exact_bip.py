"""
exact_bip.py -- exact data for the symmetry-reduced NPA level-1+AB relaxation of CGLMP_d (all states).

Chain generators g = 0: A1, 1: B1, 2: A2, 3: B2 (Alice = {0,2}, Bob = {1,3}); letters (g, n), n in 1..d-1.
Moments L(alpha, beta) = <psi| alpha(A) (x) beta(B) |psi>.
Symmetries imposed on moments (valid after symmetrising a strategy; S is invariant):
   reversal (complex conjugation):  L(a, b) = L(rev a, rev b)
   rho (chain rotation, swaps parties): A1->B1, B1->A2, A2->B2, B2-> w A1
   sigma (chain reflection, swaps parties): g -> 3-g, n -> -n
   hermiticity: L(a,b)^* = L(a^dag, b^dag)
Objective (Gill form): S = 2(d-1) + sum_n c_n [L(A1^n,B1^-n) + L(A2^-n,B1^n) + L(A2^n,B2^-n) + w^-n L(A1^-n,B2^n)].
Optimal strategy (conjectured, Acin et al.): DKZ operators A1 = R1, B1 = R2^T, A2 = R3, B2 = R4^T
(R from exact_tracial.dkz) and psi = sum_k v_k |kk>, v = Perron vector of K_d; S* = 2d-1-lambda_max(K_d).
"""
import sys
import pickle
import time
import mpmath as mp
from ncwords import reduce_word, dagger
from cyclo import Field
from extfield import Ext
from exact_tracial import dkz, gp_mul, gp_pow, gp_transpose, nullspace


def is_alice(g):
    return g in (0, 2)


class SetupB:
    def __init__(self, d, words):
        self.d, self.N = d, 4 * d
        self.F = Field(self.N)
        self.words = words
        self.classes, self.cache = {}, {}
        self.nvar = 0
        self.var_info = []

    def norm_pair(self, a, b):
        return (reduce_word(a, self.d), reduce_word(b, self.d))

    def subst(self, w, kind):
        """apply rho or sigma letterwise; returns (exp, word)"""
        d, N = self.d, self.N
        e, out = 0, []
        for (g, n) in w:
            if kind == 'rho':
                if g == 3:
                    e += 4 * n
                    out.append((0, n))
                else:
                    out.append((g + 1, n))
            else:
                out.append((3 - g, (-n) % d))
        return e % N, tuple(out)

    def gens(self, key):
        a, b = key
        out = [(0, self.norm_pair(tuple(reversed(a)), tuple(reversed(b))))]
        for kind in ('rho', 'sigma'):
            e1, a1 = self.subst(a, kind)
            e2, b1 = self.subst(b, kind)
            # a1 is now a Bob word, b1 an Alice word
            out.append(((e1 + e2) % self.N, self.norm_pair(b1, a1)))
        return out

    def orbit(self, key):
        N = self.N
        seen = {key: 0}
        stack = [key]
        bad = False
        while stack:
            u = stack.pop()
            for (e, v) in self.gens(u):
                ev = (seen[u] + e) % N
                if v in seen:
                    if seen[v] != ev:
                        bad = True
                else:
                    seen[v] = ev
                    stack.append(v)
        return seen, bad

    def moment(self, a, b):
        F, d, N = self.F, self.d, self.N
        key = self.norm_pair(a, b)
        if key in self.cache:
            return self.cache[key]
        orb, bad = self.orbit(key)
        kd = self.norm_pair(dagger(key[0], d), dagger(key[1], d))
        orb_d, bad_d = self.orbit(kd)
        if bad or bad_d:
            self.cache[key] = {}
            return {}
        rep = min(set(orb) | set(orb_d))
        if rep == ((), ()):
            res = {-1: F.z(orb[rep])} if rep in orb else {-1: F.z(-orb_d[rep])}
            self.cache[key] = res
            return res
        cls = self.classes.get(rep)
        if cls is None:
            rd = self.norm_pair(dagger(rep[0], d), dagger(rep[1], d))
            orb_rd, _ = self.orbit(rd)
            if rep in orb_rd:
                k = orb_rd[rep] % N
                f = F.i_unit() if k == (2 * d) % N else F.add(F.one, F.z(-k))
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
        Lrep = {cls[1]: cls[2]} if cls[0] == 'real' else {cls[1]: F.one, cls[2]: F.i_unit()}
        if rep in orb:
            ph = F.z(orb[rep])
            res = {v: F.mul(ph, c) for v, c in Lrep.items()}
        else:
            ph = F.z(-orb_d[rep])
            res = {v: F.mul(ph, F.conj(c)) for v, c in Lrep.items()}
        self.cache[key] = res
        return res

    def inner(self, x, y):
        """<x|y> = L(x^dag y) for elements = list of (coef, (a, b))"""
        F, d = self.F, self.d
        acc = {}
        for (cx, (ax, bx)) in x:
            cxc = F.conj(cx)
            for (cy, (ay, by)) in y:
                c0 = F.mul(cxc, cy)
                for v, c in self.moment(dagger(ax, d) + ay, dagger(bx, d) + by).items():
                    acc[v] = F.add(acc.get(v, F.zero), F.mul(c0, c))
        return {v: c for v, c in acc.items() if not F.is_zero(c)}

    def rho_word(self, key):
        a, b = key
        e1, a1 = self.subst(a, 'rho')
        e2, b1 = self.subst(b, 'rho')
        return (e1 + e2) % self.N, self.norm_pair(b1, a1)

    def blocks(self):
        F, d, N = self.F, self.d, self.N
        bysector = {}
        for u in self.words:
            q = (sum(n for g, n in u[0]) + sum(n for g, n in u[1])) % d
            bysector.setdefault(q, []).append(u)
        out = []
        for q, ws in sorted(bysector.items()):
            wset = set(ws)
            orbits, done = [], set()
            for u in ws:
                if u in done:
                    continue
                orb = []
                e, cur = 0, u
                for k in range(4):
                    orb.append((e, cur))
                    done.add(cur)
                    e2, cur = self.rho_word(cur)
                    e = (e + e2) % N
                    assert cur in wset or cur == u
                orbits.append(orb)
            for t in range(4):
                mu_exp = (q + d * t) % N
                els = []
                for orb in orbits:
                    if orb[0][1] == ((), ()):
                        if mu_exp == 0:
                            els.append([(F.one, ((), ()))])
                        continue
                    els.append([(F.z(e - k * mu_exp), wd) for k, (e, wd) in enumerate(orb)])
                els = [el for el in els if self.inner(el, el)]
                if els:
                    out.append(((q, t), els))
        return out

    def objective(self):
        F, d = self.F, self.d
        acc = {-1: F.from_int(2 * (d - 1))}
        for n in range(1, d):
            cn = F.neg(F.inv(F.sub(F.one, F.z(-4 * n))))
            terms = [(cn, (((0, n),), ((1, (-n) % d),))), (cn, (((2, (-n) % d),), ((1, n),))),
                     (cn, (((2, n),), ((3, (-n) % d),))), (F.mul(cn, F.z(-4 * n)), (((0, (-n) % d),), ((3, n),)))]
            for (c, (a, b)) in terms:
                for v, cc in self.moment(a, b).items():
                    acc[v] = F.add(acc.get(v, F.zero), F.mul(c, cc))
        for v, c in acc.items():
            assert F.is_real(c), ("objective coefficient not real", v)
        const = acc.pop(-1)
        return {v: c for v, c in acc.items() if not F.is_zero(c)}, const


def level_1AB_words(d):
    A = [()] + [((g, n),) for g in (0, 2) for n in range(1, d)]
    B = [()] + [((g, n),) for g in (1, 3) for n in range(1, d)]
    return [(a, b) for a in A for b in B]


def perron_field(d, F):
    """return Ext L = F(lambda*), lambda* = lambda_max(K_d), and the Perron vector v (v_0 = 1) over L.
    uses the symmetric block of K_d under k -> d-1-k (sizes h = ceil(d/2))."""
    import sympy as sp
    x = sp.symbols('x')
    N = F.N
    # sec(pi j/(2d)) as F element: 1/cos = 2/(z^j + z^-j) with z = zeta_{4d}
    def secF(j):
        return F.div(F.from_int(2), F.add(F.z(j), F.z(-j)))
    h = (d + 1) // 2
    Ks = [[None] * h for _ in range(h)]
    for i in range(h):
        for j in range(h):
            v = secF(abs(i - j))
            if j != d - 1 - j:
                v = F.add(v, secF(abs(i - (d - 1 - j))))
            Ks[i][j] = v
    # characteristic polynomial det(t I - Ks) via Faddeev-LeVerrier over F (monic)
    n = h
    I = [[F.one if i == j else F.zero for j in range(n)] for i in range(n)]
    def matmul(A, B):
        return [[sum_F(F, [F.mul(A[i][k], B[k][j]) for k in range(n)]) for j in range(n)] for i in range(n)]
    coeffs = [F.one]
    M = [[F.zero] * n for _ in range(n)]
    Ak = None
    c_prev = F.one
    Mk = [[F.zero] * n for _ in range(n)]
    cs = [F.one]
    for k in range(1, n + 1):
        # M_k = A M_{k-1} + c_{n-k+1} I
        AM = matmul(Ks, Mk)
        Mk = [[F.add(AM[i][j], cs[-1] if i == j else F.zero) for j in range(n)] for i in range(n)]
        AMk = matmul(Ks, Mk)
        tr = sum_F(F, [AMk[i][i] for i in range(n)])
        ck = F.neg(F.div(tr, F.from_int(k)))
        cs.append(ck)
    # polynomial t^n + c1 t^{n-1} + ... + cn ; low->high list
    p = [cs[n - i] for i in range(n)] + [F.one]
    # numeric largest root
    mp.mp.dps = 80
    pc = [F.to_complex(c, 80).real for c in p]
    roots = mp.polyroots(list(reversed(pc)), maxsteps=200, extraprec=200)
    lam = max(mp.re(r) for r in roots)
    L = Ext(F, p, lam)
    # Perron vector on symmetric block: solve (Ks - t I) u = 0 with u_0 = 1 over L
    t = L.t() if L.m >= 2 else None
    Aml = [[L.sub(L.from_F(Ks[i][j]), t if i == j else L.zero) for j in range(n)] for i in range(n)]
    # set u_0 = 1: solve rows 1..n-1 for u_1..u_{n-1}:  sum_{j>=1} A_ij u_j = -A_i0
    m = n - 1
    if m == 0:
        u = [L.one]
    else:
        Sys = [[Aml[i][j] for j in range(1, n)] + [L.neg(Aml[i][0])] for i in range(1, n)]
        for c in range(m):
            pr = next(r for r in range(c, m) if not L.is_zero(Sys[r][c]))
            Sys[c], Sys[pr] = Sys[pr], Sys[c]
            iv = L.inv(Sys[c][c])
            Sys[c] = [L.mul(v, iv) for v in Sys[c]]
            for r in range(m):
                if r != c and not L.is_zero(Sys[r][c]):
                    f = Sys[r][c]
                    Sys[r] = [L.sub(vr, L.mul(f, vc)) for vr, vc in zip(Sys[r], Sys[c])]
        u = [L.one] + [Sys[i][m] for i in range(m)]
        # check row 0
        r0 = L.zero
        for j in range(n):
            r0 = L.add(r0, L.mul(Aml[0][j], u[j]))
        assert L.is_zero(r0), "Perron vector equation fails"
    v = [u[k] if k < h else u[d - 1 - k] for k in range(d)]
    return L, v, p, lam


def sum_F(F, lst):
    acc = F.zero
    for x in lst:
        acc = F.add(acc, x)
    return acc


def strategy_images_bip(d, F, L, v):
    """strategies as (ops list [A1,B1,A2,B2] gp-matrices, psi dict (i,j)->L)."""
    N = 4 * d
    R = dkz(d)
    ops = [R[0], gp_transpose(R[1], N), R[2], gp_transpose(R[3], N)]
    psi = {(k, k): v[k] for k in range(d)}

    def key(S):
        ops, psi = S
        return (tuple((tuple(p), tuple(e)) for (p, e) in ops),
                tuple(sorted((k, tuple(tuple(c[0]) + (c[1],) for c in val)) for k, val in psi.items())))

    def swap(psi):
        return {(j, i): val for (i, j), val in psi.items()}

    def rho_s(S):
        ops, psi = S
        p0, e0 = ops[0]
        return ([ops[1], ops[2], ops[3], (p0, [(x + 4) % N for x in e0])], swap(psi))

    def sigma_s(S):
        ops, psi = S
        inv = lambda A: gp_pow(A, d - 1, d, N)
        return ([inv(ops[3]), inv(ops[2]), inv(ops[1]), inv(ops[0])], swap(psi))

    def rev_s(S):
        ops, psi = S
        return ([gp_transpose(A, N) for A in ops], {k: L.conj(val) for k, val in psi.items()})

    S0 = (ops, psi)
    seen = {key(S0): S0}
    stack = [S0]
    while stack:
        S = stack.pop()
        for T in (rho_s(S), sigma_s(S), rev_s(S)):
            k = key(T)
            if k not in seen:
                seen[k] = T
                stack.append(T)
    return list(seen.values())


def apply_word(ops, w, d, N):
    M = (list(range(d)), [0] * d)
    for (g, n) in w:
        M = gp_mul(M, gp_pow(ops[g], n, d, N), N)
    return M


def eval_element_bip(el, S, d, F, L):
    """vector (dict (i,j)->L) of element applied to the state."""
    N = 4 * d
    ops, psi = S
    out = {}
    for (c, (a, b)) in el:
        Ma = apply_word(ops, a, d, N)
        Mb = apply_word(ops, b, d, N)
        for (i, j), val in psi.items():
            ph = F.mul(c, F.z(Ma[1][i] + Mb[1][j]))
            key = (Ma[0][i], Mb[0][j])
            out[key] = L.add(out.get(key, L.zero), L.mulF(val, ph))
    return {k: v for k, v in out.items() if not L.is_zero(v)}


def nullspace_L(L, rows, ncols):
    R = [dict(r) for r in rows if r]
    pivots, prow = [], []
    for c in range(ncols):
        idx = next((i for i, r in enumerate(R) if c in r and not L.is_zero(r[c])), None)
        if idx is None:
            continue
        r = R.pop(idx)
        iv = L.inv(r[c])
        r = {k: L.mul(v, iv) for k, v in r.items()}
        def elim(s):
            if c not in s:
                return s
            f = s[c]
            t = dict(s)
            for k, v in r.items():
                t[k] = L.sub(t.get(k, L.zero), L.mul(f, v))
            return {k: v for k, v in t.items() if not L.is_zero(v)}
        R = [x for x in (elim(s) for s in R) if x]
        prow = [elim(s) for s in prow]
        pivots.append(c)
        prow.append(r)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        x = [L.zero] * ncols
        x[fc] = L.one
        for pc, r in zip(pivots, prow):
            if fc in r:
                x[pc] = L.neg(r[fc])
        basis.append(x)
    return basis


def build_all(d, verbose=True):
    t0 = time.time()
    st = SetupB(d, level_1AB_words(d))
    F = st.F
    blocks = st.blocks()
    bdata = []
    for (lab, els) in blocks:
        n = len(els)
        ent = {(i, j): st.inner(els[i], els[j]) for i in range(n) for j in range(i, n)}
        bdata.append((lab, els, ent))
    obj, const = st.objective()
    if verbose:
        print(f"d={d}: {len(blocks)} blocks, sizes {[len(e) for _, e in blocks]}, nvar={st.nvar} [{time.time()-t0:.1f}s]", flush=True)
    L, v, p, lam = perron_field(d, F)
    if verbose:
        print(f"   lambda* = {mp.nstr(lam, 25)}, extension degree {L.m} over Q(zeta_{4*d}) [{time.time()-t0:.1f}s]", flush=True)
    # exact objective value at the conjectured optimum: S* = 2d-1-lambda*  (checked numerically + via eigen-equation)
    imgs = strategy_images_bip(d, F, L, v)
    kernels = []
    for (lab, els, ent) in bdata:
        rows = []
        for S in imgs:
            vecs = [eval_element_bip(el, S, d, F, L) for el in els]
            keys = set()
            for vv in vecs:
                keys |= set(vv.keys())
            for kk in keys:
                rows.append({j: vecs[j][kk] for j in range(len(els)) if kk in vecs[j]})
        kernels.append(nullspace_L(L, rows, len(els)))
    if verbose:
        print(f"   images={len(imgs)} kernel dims {[len(k) for k in kernels]} [{time.time()-t0:.1f}s]", flush=True)
    # S* as element of L: 2d-1 - t
    Sstar = L.sub(L.from_F(F.from_int(2 * d - 1)), L.t() if L.m >= 2 else L.from_F(F.zero))
    return dict(d=d, N=st.N, st_nvar=st.nvar, var_info=st.var_info, bdata=bdata, obj=obj, const=const,
                p=p, lam=str(lam), v=v, kernels=kernels, Sstar=Sstar)


if __name__ == "__main__":
    for d in [int(t) for t in sys.argv[1:]]:
        data = build_all(d)
        with open(f"exact_bip_d{d}.pkl", "wb") as fh:
            pickle.dump(data, fh)
