"""
verify_rigidity.py -- checks, for a maximally-entangled certificate cert_tracial_d{d}.pkl, the support
properties used in the rigidity (self-testing) argument of LOG.md section 6:

 (a) sector-0 blocks with kappa != 1 (labels (0,1),(0,2),(0,3)) are entirely in the SOS support
     (kernel = whole block, Y_b > 0 checked by verify_tracial.py)  => all chain modes
     sum_k kappa^-k rho^k(R1^n R2^-n) vanish on optimal strategies  => equal links (E_n);
 (b) in the kappa = 1 sector-0 block the SOS support contains, for n = 1..d-1,
        T_n = sum_{k=0}^{3} rho^k( U_n^{-1} - (1-i) z^n 1 - i z^{2n} U_n ),   U_n = R1^n R2^-n,
     i.e. (with (E_n)) the two-eigenvalue relation (U_n - z^-n)(U_n - i z^-n) = 0, z = zeta_{4d}.
Together with the analytic argument in LOG.md this proves: every maximally entangled strategy attaining
the bound is unitarily equivalent to DKZ (x) 1_M.
"""
import sys
import pickle
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "builders"))  # cyclo.py lives in builders/
from cyclo import Field
from ncwords import reduce_word


def rho(w, d):
    e, out = 0, []
    for (g, n) in w:
        if g == 3:
            e += 4 * n
            out.append((0, n))
        else:
            out.append((g + 1, n))
    return e, tuple(out)


def poly_add(F, P, w, c):
    P[w] = F.add(P.get(w, F.zero), c)


def elem_poly(F, el, d):
    P = {}
    for (c, w) in el:
        poly_add(F, P, reduce_word(w, d), c)
    return {w: c for w, c in P.items() if not F.is_zero(c)}


def rank(F, rows, ncols):
    """rank of list of row dicts col->F (exact)."""
    R = [dict(r) for r in rows if r]
    rk = 0
    for c in range(ncols):
        idx = next((i for i, r in enumerate(R) if c in r and not F.is_zero(r[c])), None)
        if idx is None:
            continue
        r = R.pop(idx)
        iv = F.inv(r[c])
        r = {k: F.mul(v, iv) for k, v in r.items()}
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
        rk += 1
    return rk


def check(d):
    cert = pickle.load(open(f"cert_tracial_d{d}.pkl", "rb"))
    F = Field(4 * d)
    z = lambda k: F.z(k)
    iu = F.i_unit()
    labels = [lab for (lab, els, ent) in cert["bdata"]]
    # (a)
    for t in (1, 2, 3):
        b = labels.index((0, t))
        n_el = len(cert["bdata"][b][1])
        k = len(cert["kernels"][b])
        assert k == n_el, f"block (0,{t}) not fully in the support"
        assert b in cert["keep"]
    print(f"d={d}: (a) sector-0 blocks kappa=i,-1,-i: kernel = whole block (sizes "
          f"{[len(cert['bdata'][labels.index((0, t))][1]) for t in (1, 2, 3)]})")
    # (b)
    b = labels.index((0, 0))
    els = cert["bdata"][b][1]
    ker = cert["kernels"][b]
    polys = [elem_poly(F, el, d) for el in els]
    words = sorted(set(w for P in polys for w in P))
    widx = {w: i for i, w in enumerate(words)}
    # express kernel vectors as polynomials
    kpolys = []
    for vec in ker:
        P = {}
        for x, Pe in zip(vec, polys):
            if F.is_zero(x):
                continue
            for w, c in Pe.items():
                poly_add(F, P, w, F.mul(x, c))
        kpolys.append({w: c for w, c in P.items() if not F.is_zero(c)})
    base_rank = rank(F, [{widx[w]: c for w, c in P.items()} for P in kpolys], len(words))
    ok = True
    for n in range(1, d):
        T = {}
        U = ((0, n), (1, (-n) % d))
        Uinv = ((1, n), (0, (-n) % d))
        c_one = F.mul(F.sub(F.one, iu), z(n))          # (1-i) z^n
        c_U = F.mul(iu, z(2 * n))                      # i z^{2n}
        for (coef, w) in ((F.one, Uinv), (F.neg(c_U), U), (F.neg(c_one), ())):
            e, cur = 0, w
            for k in range(4):
                poly_add(F, T, reduce_word(cur, d), F.mul(coef, z(e)))
                e2, cur = rho(cur, d)
                e += e2
        T = {w: c for w, c in T.items() if not F.is_zero(c)}
        for w in T:
            if w not in widx:
                ok = False
                print(f"   T_{n}: word {w} not in block -> cannot lie in the support")
                break
        else:
            r2 = rank(F, [{widx[w]: c for w, c in P.items()} for P in kpolys] + [{widx[w]: c for w, c in T.items()}], len(words))
            if r2 != base_rank:
                ok = False
                print(f"   T_{n} NOT in the support span")
    assert ok
    print(f"d={d}: (b) all two-eigenvalue relations T_1..T_{d-1} lie in the SOS support of block (0,0)")
    print(f"d={d}: RIGIDITY support conditions verified")


if __name__ == "__main__":
    for d in [int(t) for t in sys.argv[1:]]:
        check(d)
