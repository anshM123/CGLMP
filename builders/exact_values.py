"""
exact_values.py -- exact algebraic data for the conjectured CGLMP optima.

K_d = [sec(pi (k-l)/(2d))]_{k,l=0}^{d-1}.
  maximally entangled (DKZ):   Pi_ME   = 2 - 1^T K_d 1 / d^2
  all states (Acin et al.):    Pi_ADGL = 2 - lambda_max(K_d)/d
  CGLMP value I = 4 - 2(d Pi - 1)/(d-1).
Minimal polynomials over Q are computed exactly (resultant with the minimal polynomial of
c = cos(pi/(2d)), factorisation over Q, selection of the factor vanishing at the numerical value).
"""
import sys
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
x, c = sp.symbols('x c')


def cheb(j):
    return sp.chebyshevt(j, c)


def data(d):
    mc = sp.minimal_polynomial(sp.cos(sp.pi / (2 * d)), c)
    cnum = mp.cos(mp.pi / (2 * d))
    # symmetric block of K_d under k -> d-1-k : basis e_k + e_{d-1-k}
    half = (d + 1) // 2
    sec = lambda j: 1 / cheb(abs(j))
    Ks = sp.zeros(half, half)
    for i in range(half):
        for j in range(half):
            # (K s_j)_i where s_j = e_j + e_{d-1-j} (or e_j if j = d-1-j)
            v = sec(i - j)
            if j != d - 1 - j:
                v += sec(i - (d - 1 - j))
            Ks[i, j] = sp.together(v)
    # characteristic polynomial numerator in x and c
    p = sp.together((x * sp.eye(half) - Ks).det(method='berkowitz'))
    num, den = sp.fraction(p)
    num = sp.Poly(sp.expand(num), x, c)
    R = sp.resultant(num.as_expr(), mc, c)
    R = sp.Poly(R, x)
    # numeric lambda_max
    Kn = mp.matrix(d, d)
    for i in range(d):
        for j in range(d):
            Kn[i, j] = 1 / mp.cos(mp.pi * (i - j) / (2 * d))
    ev = mp.eigsy(Kn)[0]
    lam = max(ev[i] for i in range(d))
    facs = sp.factor_list(R.as_expr())[1]
    best = None
    for f, m in facs:
        val = abs(sp.Poly(f, x).eval(sp.Float(str(lam), 60)))
        if best is None or val < best[0]:
            best = (val, sp.Poly(f, x))
    mpl = best[1]
    # maximally entangled value s = 1^T K 1 / d^2 in Q(c)
    s = sp.together(sum(sec(i - j) for i in range(d) for j in range(d)) / d ** 2)
    ms = sp.minimal_polynomial(s.subs(c, sp.cos(sp.pi / (2 * d))), x)
    snum = sum(1 / mp.cos(mp.pi * (i - j) / (2 * d)) for i in range(d) for j in range(d)) / d ** 2
    return mpl, lam, ms, snum, mc


if __name__ == "__main__":
    ds = [int(t) for t in sys.argv[1:]] or [2, 3, 4, 5, 6]
    for d in ds:
        mpl, lam, ms, snum, mc = data(d)
        Pi_ad = 2 - lam / d
        I_ad = 4 - 2 * (d * Pi_ad - 1) / (d - 1)
        Pi_me = 2 - snum
        I_me = 4 - 2 * (d * Pi_me - 1) / (d - 1)
        print(f"d={d}: minpoly of cos(pi/2d): {mc}")
        print(f"   lambda_max(K_d) = {mp.nstr(lam, 40)}  (degree {mpl.degree()})")
        print(f"   minpoly: {mpl.as_expr()}")
        print(f"   Pi_ADGL = {mp.nstr(Pi_ad, 40)}   I_ADGL = {mp.nstr(I_ad, 40)}")
        print(f"   s_ME = 1^T K 1/d^2 = {mp.nstr(snum, 40)}  minpoly: {ms}")
        print(f"   Pi_ME = {mp.nstr(Pi_me, 40)}   I_ME = {mp.nstr(I_me, 40)}", flush=True)
