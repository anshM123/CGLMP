"""
check_support.py -- INDEPENDENT exact checker of the rigidity support conditions for the
maximally-entangled CGLMP certificates  ../certificates/maxent/cert_tracial_d{d}.pkl  (read only).

It imports NO module of the lab (no cyclo.py, ncwords.py, verify_*.py, builders).  Re-implemented here:
  * Q(zeta_N), N = 4d, as Q[x]/Phi_N(x); Phi_N from the Moebius product formula
    Phi_N = prod_{k|N} (x^k - 1)^{mu(N/k)}; inverses by the extended Euclidean algorithm;
  * words in the free product Z_d * Z_d * Z_d * Z_d (letters (g, p), g = 0..3 <-> R1..R4, 1 <= p <= d-1);
  * the target polynomials, written out link by link (no generic symmetry routine);
  * positivity of the Gram matrices by interval Cholesky (mpmath.iv) of the real 2k x 2k embedding;
  * exact membership by sparse Gaussian elimination, followed by an exact re-check of the found combination.

Mathematical meaning (see RIGIDITY.md, Lemma 2 and Proposition 3).  The certificate gives
    S_op - lam  =  sum_b sum_{ij} (X_b)_{ij} e_{b,j}^* e_{b,i}   (mod commutators and symmetry relations),
    X_b = N_b Y_b N_b^+ .
If every Y_b is positive definite, the support polynomials are f_{b,u} = sum_i u_i e_{b,i}, u in range(N_b), and
every optimal strategy R satisfies f(R) = 0 for all f in their linear span.  The checker verifies:

 (S0) lam equals the exact DKZ value S(DKZ) (own evaluation with generalised permutation matrices);
 (S1) every kept Y_b is Hermitian and positive definite (rigorous interval Cholesky, 320 bits);
 (S2) for n = 1..d-1, t = 1,2,3:  M_{n,t} = sum_{k=0}^{3} i^{-tk} L_n^(k) lies in the support span, where
        L_n^(0) = R1^n R2^-n, L_n^(1) = R2^n R3^-n, L_n^(2) = R3^n R4^-n, L_n^(3) = w^-n R4^n R1^-n;
      (=> equal links (E_n));
 (S3) for n = 1..d-1:  T_n = sum_k [ (L_n^(k))^-1 - (1-i) z^n - i z^{2n} L_n^(k) ]  lies in the support span;
 (S3') stronger, non-symmetrised: V_n = (L_n^(0))^-1 - (1-i) z^n - i z^{2n} L_n^(0) lies in the support span
      (=> (U_n - z^-n)(U_n - i z^-n) = 0 on the first link directly);
 (S4) consistency: every target polynomial and every sector-0 support polynomial vanishes exactly on DKZ.
Here z = zeta_{4d}, w = z^4, i = z^d.

Usage:  python -B check_support.py 3 4 5 6 7 8 9
"""
import sys
import os
import time
import pickle
from fractions import Fraction
try:
    from gmpy2 import mpq as Q          # fast exact rationals (falls back to fractions)
except ImportError:                      # pragma: no cover
    Q = Fraction
import mpmath

LAB = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates", "maxent"))


# ============================ integer / rational polynomials ============================
def mobius(n):
    res, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    if m > 1:
        res = -res
    return res


def ipoly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out


def ipoly_divexact(a, b):
    """exact division of integer polynomials, b monic (leading coefficient 1)."""
    a = list(a)
    assert b[-1] == 1
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1]
        q[i] = c
        if c:
            for j in range(len(b)):
                a[i + j] -= c * b[j]
    assert all(v == 0 for v in a), "non-exact division"
    return q


def cyclotomic_moebius(N):
    num, den = [1], [1]
    for k in range(1, N + 1):
        if N % k == 0:
            mu = mobius(N // k)
            xk1 = [-1] + [0] * (k - 1) + [1]
            if mu == 1:
                num = ipoly_mul(num, xk1)
            elif mu == -1:
                den = ipoly_mul(den, xk1)
    # normalise signs so that den is monic
    if den[-1] == -1:
        den = [-v for v in den]
        num = [-v for v in num]
    phi = ipoly_divexact(num, den)
    if phi[-1] == -1:
        phi = [-v for v in phi]
    return phi


def qpoly_trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def qpoly_divmod(a, b):
    a = qpoly_trim(a)
    b = qpoly_trim(b)
    if len(a) < len(b):
        return [], a
    q = [Q(0)] * (len(a) - len(b) + 1)
    inv_lead = Q(1) / b[-1]
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] * inv_lead
        q[i] = c
        if c:
            for j in range(len(b)):
                a[i + j] -= c * b[j]
    return q, qpoly_trim(a[:len(b) - 1])


def qpoly_sub(a, b):
    n = max(len(a), len(b))
    return qpoly_trim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)])


def qpoly_mul(a, b):
    if not a or not b:
        return []
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return qpoly_trim(out)


# ============================ the field Q(zeta_N) ============================
class Cyc:
    """Q(zeta_N) = Q[x]/Phi_N(x); elements are tuples of Q of length phi(N) (power basis)."""

    def __init__(self, N):
        self.N = N
        self.phi = cyclotomic_moebius(N)
        self.n = len(self.phi) - 1
        self.zero = tuple(Q(0) for _ in range(self.n))
        self.one = tuple(Q(1) if k == 0 else Q(0) for k in range(self.n))
        self._zp = {}

    def reduce(self, poly):
        p = [Q(v) for v in poly]
        n, phi = self.n, self.phi
        for k in range(len(p) - 1, n - 1, -1):
            c = p[k]
            if c:
                for j in range(n + 1):
                    p[k - n + j] -= c * phi[j]
        p = p[:n] + [Q(0)] * max(0, n - len(p))
        return tuple(p)

    def zp(self, k):
        k %= self.N
        if k not in self._zp:
            self._zp[k] = self.reduce([0] * k + [1])
        return self._zp[k]

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple(x - y for x, y in zip(a, b))

    def neg(self, a):
        return tuple(-x for x in a)

    def mul(self, a, b):
        out = [Q(0)] * (2 * self.n - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        out[i + j] += x * y
        return self.reduce(out)

    def is_zero(self, a):
        return all(x == 0 for x in a)

    def conj(self, a):
        acc = [Q(0)] * self.N
        for k, x in enumerate(a):
            if x:
                acc[(-k) % self.N] += x
        return self.reduce(acc)

    def inv(self, a):
        """extended Euclid: find u with u*a = 1 mod Phi_N."""
        assert not self.is_zero(a), "inverse of zero"
        # invariant: r_j = s_j * a  (mod Phi_N);  r0 = Phi_N = 0*a, r1 = a = 1*a
        r0, r1 = [Q(v) for v in self.phi], qpoly_trim(list(a))
        s0, s1 = [], [Q(1)]
        while r1:
            q, r = qpoly_divmod(r0, r1)
            r0, r1 = r1, r
            s0, s1 = s1, qpoly_sub(s0, qpoly_mul(q, s1))
        # r0 = gcd(Phi_N, a) = s0 * a mod Phi_N; Phi_N irreducible and deg a < deg Phi_N => constant
        assert len(r0) == 1, "not invertible (gcd non-constant)"
        c = Q(1) / r0[0]
        u = self.reduce([x * c for x in s0] + [0])
        assert self.mul(u, a) == self.one, "inverse check failed"
        return u

    def from_cert(self, x):
        vals, den = x
        return self.reduce([Q(v, den) for v in vals])

    def to_iv(self, a):
        """rigorous complex enclosure (re, im) as mpmath.iv intervals."""
        iv = mpmath.iv
        re, im = iv.mpf(0), iv.mpf(0)
        for k, x in enumerate(a):
            if x:
                q = iv.mpf(int(x.numerator)) / int(x.denominator)
                ang = 2 * iv.pi * k / self.N
                re += q * iv.cos(ang)
                im += q * iv.sin(ang)
        return re, im

    def to_complex(self, a):
        import cmath
        return sum(float(x) * cmath.exp(2j * cmath.pi * k / self.N) for k, x in enumerate(a) if x)


# ============================ words and polynomials ============================
def wreduce(w, d):
    out = []
    for (g, p) in w:
        p %= d
        if p == 0:
            continue
        if out and out[-1][0] == g:
            q = (out[-1][1] + p) % d
            out.pop()
            if q:
                out.append((g, q))
        else:
            out.append((g, p))
    return tuple(out)


def padd(K, P, w, c):
    if K.is_zero(c):
        return
    if w in P:
        s = K.add(P[w], c)
        if K.is_zero(s):
            del P[w]
        else:
            P[w] = s
    else:
        P[w] = c


def link(n, k, d):
    """(word, zeta-exponent) of L_n^(k)."""
    if k == 0:
        return ((0, n), (1, -n)), 0
    if k == 1:
        return ((1, n), (2, -n)), 0
    if k == 2:
        return ((2, n), (3, -n)), 0
    return ((3, n), (0, -n)), -4 * n          # w^-n R4^n R1^-n


def link_inv(n, k, d):
    """(word, zeta-exponent) of (L_n^(k))^{-1}."""
    if k == 0:
        return ((1, n), (0, -n)), 0
    if k == 1:
        return ((2, n), (1, -n)), 0
    if k == 2:
        return ((3, n), (2, -n)), 0
    return ((0, n), (3, -n)), 4 * n           # w^n R1^n R4^-n


def targets(K, d):
    """dict name -> polynomial (dict word -> field element)."""
    out = {}
    iu = K.zp(d)                              # i = zeta_{4d}^d
    one_minus_i = K.sub(K.one, iu)
    for n in range(1, d):
        for t in (1, 2, 3):
            P = {}
            for k in range(4):
                w, e = link(n, k, d)
                padd(K, P, wreduce(w, d), K.zp(e - d * t * k))
            out[("M", n, t)] = P
        P = {}
        V = {}
        c1 = K.neg(K.mul(one_minus_i, K.zp(n)))            # -(1-i) z^n
        cU = K.neg(K.mul(iu, K.zp(2 * n)))                 # -i z^{2n}
        for k in range(4):
            wi, ei = link_inv(n, k, d)
            w, e = link(n, k, d)
            padd(K, P, wreduce(wi, d), K.zp(ei))
            padd(K, P, (), c1)
            padd(K, P, wreduce(w, d), K.mul(cU, K.zp(e)))
            if k == 0:
                padd(K, V, wreduce(wi, d), K.zp(ei))
                padd(K, V, (), c1)
                padd(K, V, wreduce(w, d), K.mul(cU, K.zp(e)))
        out[("T", n)] = P
        out[("V", n)] = V
    return out


# ============================ DKZ as generalised permutation matrices ============================
def dkz(d):
    """R_{g+1} = G^{a} X^{-1} G^{-a}, a = 2-g, G = diag(zeta^k): R|k> = zeta^{e[k]} |p[k]>."""
    N = 4 * d
    R = []
    for g in range(4):
        a = 2 - g
        p = [(k - 1) % d for k in range(d)]
        e = [((-a * k) + a * ((k - 1) % d)) % N for k in range(d)]
        R.append((p, e))
    return R


def gp_mul(A, B, N):          # (A B)|k> = A (zeta^{eB[k]} |pB[k]>)
    pA, eA = A
    pB, eB = B
    return [pA[pB[k]] for k in range(len(pB))], [(eB[k] + eA[pB[k]]) % N for k in range(len(pB))]


def gp_pow(A, m, d, N):
    I = (list(range(d)), [0] * d)
    R = I
    for _ in range(m % d):
        R = gp_mul(A, R, N)
    return R


def eval_poly(K, P, R, d):
    """matrix (dict (row,col) -> field) of the polynomial P at the generalised-permutation tuple R."""
    N = 4 * d
    M = {}
    for w, c in P.items():
        G = (list(range(d)), [0] * d)
        for (g, p) in w:
            G = gp_mul(G, gp_pow(R[g], p, d, N), N)
        for k in range(d):
            key = (G[0][k], k)
            M[key] = K.add(M.get(key, K.zero), K.mul(c, K.zp(G[1][k])))
    return {k: v for k, v in M.items() if not K.is_zero(v)}


def S_value(K, R, d):
    """exact S(R) = 2(d-1) + sum_n c_n tr(L_n^(0)+...+L_n^(3)), tr normalised (value in Q(zeta))."""
    val = K.reduce([2 * (d - 1)])
    for n in range(1, d):
        cn = K.neg(K.inv(K.sub(K.one, K.zp(-4 * n))))
        for k in range(4):
            w, e = link(n, k, d)
            M = eval_poly(K, {wreduce(w, d): K.zp(e)}, R, d)
            tr = K.zero
            for (r, c), v in M.items():
                if r == c:
                    tr = K.add(tr, v)
            val = K.add(val, K.mul(cn, K.mul(tr, K.reduce([Q(1) / d]))))
    return val


# ============================ positivity: interval Cholesky ============================
def interval_pd(K, Y, prec=320):
    """rigorous test that the Hermitian matrix Y (list of lists of field elements) is positive definite:
    Cholesky of the real embedding [[A,-B],[B,A]] in interval arithmetic; returns min pivot lower bound or None."""
    iv = mpmath.iv
    old = iv.prec
    iv.prec = prec
    try:
        k = len(Y)
        for a in range(k):
            for b in range(k):
                assert K.is_zero(K.sub(Y[a][b], K.conj(Y[b][a]))), "Y not Hermitian"
        A = [[None] * k for _ in range(k)]
        B = [[None] * k for _ in range(k)]
        for a in range(k):
            for b in range(k):
                A[a][b], B[a][b] = K.to_iv(Y[a][b])
        n = 2 * k
        M = [[None] * n for _ in range(n)]
        for a in range(k):
            for b in range(k):
                M[a][b] = A[a][b]
                M[a][b + k] = -B[a][b]
                M[a + k][b] = B[a][b]
                M[a + k][b + k] = A[a][b]
        L = [[iv.mpf(0)] * n for _ in range(n)]
        minpiv = None
        for j in range(n):
            s = M[j][j]
            for m in range(j):
                s -= L[j][m] * L[j][m]
            if not (s.a > 0):
                return None
            minpiv = s.a if minpiv is None else min(minpiv, s.a)
            Ljj = iv.sqrt(s)
            L[j][j] = Ljj
            for r in range(j + 1, n):
                t = M[r][j]
                for m in range(j):
                    t -= L[r][m] * L[j][m]
                L[r][j] = t / Ljj
        return minpiv
    finally:
        iv.prec = old


# ============================ exact membership ============================
def solve_membership(K, vecs, target):
    """vecs: list of polynomials (dict word -> field); target: polynomial.
    returns coefficient list c with sum_j c_j vecs[j] == target (exact), or None if target not in span."""
    # rows = words; build column-wise sparse representation, eliminate on the augmented system
    words = sorted(set(w for v in vecs for w in v) | set(target))
    widx = {w: i for i, w in enumerate(words)}
    m = len(vecs)
    # system:  sum_j c_j vecs[j][w] = target[w]  for every word w; unknowns c_j
    rows = []
    for w in words:
        row = {}
        for j, v in enumerate(vecs):
            if w in v:
                row[j] = v[w]
        rows.append((row, target.get(w, K.zero)))
    pivots = []           # (col, row dict, rhs) in reduced form
    for (row, rhs) in rows:
        row = dict(row)
        # eliminate existing pivots
        for (c, prow, prhs) in pivots:
            if c in row:
                f = row[c]
                for kk, vv in prow.items():
                    nv = K.sub(row.get(kk, K.zero), K.mul(f, vv))
                    if K.is_zero(nv):
                        row.pop(kk, None)
                    else:
                        row[kk] = nv
                rhs = K.sub(rhs, K.mul(f, prhs))
        if not row:
            if not K.is_zero(rhs):
                return None          # inconsistent: target not in span
            continue
        c = min(row)
        ic = K.inv(row[c])
        row = {kk: K.mul(vv, ic) for kk, vv in row.items()}
        rhs = K.mul(rhs, ic)
        # back-substitute into existing pivots to keep them reduced
        new_piv = []
        for (c2, prow, prhs) in pivots:
            if c in prow:
                f = prow[c]
                prow = dict(prow)
                for kk, vv in row.items():
                    nv = K.sub(prow.get(kk, K.zero), K.mul(f, vv))
                    if K.is_zero(nv):
                        prow.pop(kk, None)
                    else:
                        prow[kk] = nv
                prhs = K.sub(prhs, K.mul(f, rhs))
            new_piv.append((c2, prow, prhs))
        new_piv.append((c, row, rhs))
        pivots = new_piv
    coef = [K.zero] * m
    for (c, prow, prhs) in pivots:
        coef[c] = prhs            # free variables set to 0 (reduced form: other pivot cols eliminated)
    # exact re-check
    acc = {}
    for j, v in enumerate(vecs):
        if K.is_zero(coef[j]):
            continue
        for w, x in v.items():
            padd(K, acc, w, K.mul(coef[j], x))
    keys = set(acc) | set(target)
    assert all(K.is_zero(K.sub(acc.get(w, K.zero), target.get(w, K.zero))) for w in keys), "re-check failed"
    return coef


# ============================ (P1) the SOS identity, re-checked independently ============================
# A symmetry is stored as (anti, imgs), imgs[g] = (e, j, s): alpha(R_g) = zeta^e R_j^s (s = +-1); anti = True for
# anti-automorphisms (they reverse products).  G = group generated by rho, sigma, rev; enumerated explicitly.
def sym_compose(a, b, N):
    """(a o b)(R_g) = a(b(R_g))."""
    anti_a, ia = a
    anti_b, ib = b
    out = []
    for g in range(4):
        e, j, s = ib[g]
        e2, j2, s2 = ia[j]
        out.append(((e + s * e2) % N, j2, s * s2))
    return (anti_a != anti_b, tuple(out))


def sym_group(d):
    N = 4 * d
    rho = (False, ((0, 1, 1), (0, 2, 1), (0, 3, 1), (4, 0, 1)))
    sigma = (False, tuple((0, 3 - g, -1) for g in range(4)))
    rev = (True, tuple((0, g, 1) for g in range(4)))
    ident = (False, tuple((0, g, 1) for g in range(4)))
    G = {ident}
    frontier = [ident]
    while frontier:
        new = []
        for x in frontier:
            for gen in (rho, sigma, rev):
                y = sym_compose(gen, x, N)
                if y not in G:
                    G.add(y)
                    new.append(y)
        frontier = new
    return sorted(G), (rho, sigma, rev)


def sym_apply(a, w, d):
    """alpha(word) = zeta^ph * word'; returns (ph mod N, reduced word')."""
    anti, imgs = a
    N = 4 * d
    ph, out = 0, []
    for (g, p) in w:
        e, j, s = imgs[g]
        ph += e * p
        out.append((j, (s * p) % d))
    if anti:
        out.reverse()
    return ph % N, wreduce(tuple(out), d)


def cyc_canon(w, d):
    w = list(wreduce(w, d))
    while len(w) >= 2 and w[0][0] == w[-1][0]:
        g, q = w[0][0], (w[0][1] + w[-1][1]) % d
        w = w[1:-1]
        if q:
            w = [(g, q)] + w
        w = list(wreduce(tuple(w), d))
    w = tuple(w)
    if len(w) <= 1:
        return w
    return min(w[i:] + w[:i] for i in range(len(w)))


class SymClass:
    """for G-invariant tracial L:  L(w) = zeta^e L(rep)  -> (e, rep), or None if L(w) = 0 is forced."""

    def __init__(self, d, G):
        self.d, self.N, self.G = d, 4 * d, G
        self.memo = {}

    def __call__(self, w):
        d, N = self.d, self.N
        c = cyc_canon(w, d)
        if c in self.memo:
            return self.memo[c]
        seen = {}
        bad = False
        for a in self.G:
            ph, u = sym_apply(a, c, d)          # L(c) = L(a(c)) = zeta^ph L(u)
            u = cyc_canon(u, d)
            if u in seen and seen[u] != ph:
                bad = True
            seen.setdefault(u, ph)
        if bad:
            res = None
        else:
            rep = min(seen)
            res = (seen[rep], rep)
        self.memo[c] = res
        return res


def adj_word(w, d):
    return tuple((g, (-p) % d) for (g, p) in reversed(w))


def S_poly(K, d, lam):
    P = {(): K.sub(K.reduce([2 * (d - 1)]), lam)}
    for n in range(1, d):
        cn = K.neg(K.inv(K.sub(K.one, K.zp(-4 * n))))
        for k in range(4):
            w, e = link(n, k, d)
            padd(K, P, wreduce(w, d), K.mul(cn, K.zp(e)))
    return P


def check_identity(d, verbose=True):
    t0 = time.time()
    with open(os.path.join(LAB, f"cert_tracial_d{d}.pkl"), "rb") as fh:
        cert = pickle.load(fh)
    N = 4 * d
    K = Cyc(N)
    G, gens = sym_group(d)
    lam = K.from_cert(cert["lam"])
    Sop = S_poly(K, d, lam)
    # (i) S_op is G-invariant modulo commutators (so S(alpha.R) = S(R) for every strategy R)
    def cyc_image(P, a):
        out = {}
        for w, c in P.items():
            ph, u = sym_apply(a, w, d)
            padd(K, out, cyc_canon(u, d), K.mul(c, K.zp(ph)))
        return out
    ident = (False, tuple((0, g, 1) for g in range(4)))
    assert ident in G
    base = cyc_image(Sop, ident)
    for a in G:
        img = cyc_image(Sop, a)
        keys = set(img) | set(base)
        assert all(K.is_zero(K.sub(img.get(k, K.zero), base.get(k, K.zero))) for k in keys), "S_op not invariant"
    # (ii) canonical images of both sides
    canon = SymClass(d, G)
    touched = set()
    forced_zero = set()

    def add_canon(acc, w, c):
        r = canon(w)
        if r is None:
            forced_zero.add(cyc_canon(w, d))
            return
        if K.is_zero(c):
            return
        e, rep = r
        touched.add(rep)
        padd(K, acc, rep, K.mul(c, K.zp(e)))
    lhs = {}
    for w, c in Sop.items():
        add_canon(lhs, w, c)
    rhs = {}
    nterms = 0
    for b in cert["keep"]:
        _lab, els, _ent = cert["bdata"][b]
        ker = [[K.from_cert(x) for x in vec] for vec in cert["kernels"][b]]
        Y = [[K.from_cert(x) for x in row] for row in cert["Y"][b]]
        k, n = len(ker), len(els)
        # X = N Y N^+,  N[i][r] = ker[r][i]
        NY = [[K.zero] * k for _ in range(n)]
        for i in range(n):
            for s in range(k):
                acc = K.zero
                for r in range(k):
                    if not K.is_zero(ker[r][i]) and not K.is_zero(Y[r][s]):
                        acc = K.add(acc, K.mul(ker[r][i], Y[r][s]))
                NY[i][s] = acc
        X = [[K.zero] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                acc = K.zero
                for s in range(k):
                    if not K.is_zero(NY[i][s]) and not K.is_zero(ker[s][j]):
                        acc = K.add(acc, K.mul(NY[i][s], K.conj(ker[s][j])))
                X[i][j] = acc
        elsK = [[(K.from_cert(c), wreduce(w, d)) for (c, w) in el] for el in els]
        # SOS = sum_{ij} X_ij e_j^* e_i  ( = sum_{rs} Y_rs f_s^* f_r with f_r = sum_i N_ir e_i )
        for i in range(n):
            for j in range(n):
                if K.is_zero(X[i][j]):
                    continue
                for (cj, wj) in elsK[j]:
                    for (ci, wi) in elsK[i]:
                        add_canon(rhs, wreduce(adj_word(wj, d) + wi, d), K.mul(X[i][j], K.mul(K.conj(cj), ci)))
                        nterms += 1
    keys = set(lhs) | set(rhs)
    bad = [kk for kk in keys if not K.is_zero(K.sub(lhs.get(kk, K.zero), rhs.get(kk, K.zero)))]
    if verbose:
        print(f"d={d}: (P1) |G| = {len(G)}; S_op G-invariant mod commutators; identity S_op - lam = SOS: "
              f"{len(touched)} symmetry classes touched ({len(forced_zero)} further cyclic words forced to 0), "
              f"{len(keys)} with nonzero net coefficient, {nterms} expanded terms: {'OK' if not bad else 'FAILED'} "
              f"[{time.time()-t0:.1f}s]", flush=True)
    assert not bad
    return True


# ============================ main check ============================
def selftest(K):
    """field sanity: zeta has exact multiplicative order N; Phi_N(zeta) = 0 numerically; mul/inv/conj agree
    with complex arithmetic on pseudo-random elements."""
    import cmath
    import random
    N = K.N
    assert K.zp(N) == K.one and K.zp(N // 2) == K.neg(K.one)
    assert all(K.zp(k) != K.one for k in range(1, N))
    zc = cmath.exp(2j * cmath.pi / N)
    assert abs(sum(c * zc ** k for k, c in enumerate(K.phi))) < 1e-9
    rng = random.Random(12345)
    for _ in range(5):
        a = K.reduce([Q(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(K.n)])
        b = K.reduce([Q(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(K.n)])
        if K.is_zero(b):
            continue
        ac, bc = K.to_complex(a), K.to_complex(b)
        assert abs(K.to_complex(K.mul(a, b)) - ac * bc) < 1e-9
        assert abs(K.to_complex(K.inv(b)) - 1 / bc) < 1e-9 * (1 + abs(1 / bc))
        assert abs(K.to_complex(K.conj(a)) - ac.conjugate()) < 1e-9


def check(d, verbose=True):
    t0 = time.time()
    path = os.path.join(LAB, f"cert_tracial_d{d}.pkl")
    with open(path, "rb") as fh:
        cert = pickle.load(fh)
    assert cert["d"] == d
    N = 4 * d
    K = Cyc(N)
    selftest(K)
    report = {"d": d}
    # ---- (S0) lam = S(DKZ)
    R = dkz(d)
    lam = K.from_cert(cert["lam"])
    Sd = S_value(K, R, d)
    assert K.is_zero(K.sub(Sd, lam)), "lam != S(DKZ)"
    # R_g^d = 1 exactly
    for g in range(4):
        P = gp_mul(gp_pow(R[g], d - 1, d, N), R[g], N)
        assert P[0] == list(range(d)) and all(x % N == 0 for x in P[1]), "R^d != 1"
    if verbose:
        print(f"d={d}: (S0) lam = S(DKZ) exactly (= {K.to_complex(lam).real:.15f}); R_g^d = 1", flush=True)
    # ---- (S1) positivity of every kept Gram block
    minpiv = None
    for b in cert["keep"]:
        Y = [[K.from_cert(x) for x in row] for row in cert["Y"][b]]
        mp_ = interval_pd(K, Y)
        assert mp_ is not None, f"block {b}: positive definiteness NOT certified"
        minpiv = mp_ if minpiv is None else min(minpiv, mp_)
    if verbose:
        print(f"d={d}: (S1) all {len(cert['keep'])} kept Gram matrices Y_b Hermitian positive definite "
              f"(interval Cholesky, min pivot lower bound {mpmath.nstr(minpiv, 5)}) [{time.time()-t0:.1f}s]",
              flush=True)
    # ---- support polynomials of the charge-0 blocks
    supp = []
    nblk = 0
    for b in cert["keep"]:
        lab_, els, _ent = cert["bdata"][b]
        ker = cert["kernels"][b]
        if not ker:
            continue
        elsK = [[(K.from_cert(c), wreduce(w, d)) for (c, w) in el] for el in els]
        charges = set(sum(p for (_g, p) in w) % d for el in elsK for (_c, w) in el)
        assert len(charges) == 1, "block not homogeneous in charge"
        if charges != {0}:
            continue
        nblk += 1
        for vec in ker:
            assert len(vec) == len(els)
            P = {}
            for x, el in zip(vec, elsK):
                xv = K.from_cert(x)
                if K.is_zero(xv):
                    continue
                for (c, w) in el:
                    padd(K, P, w, K.mul(xv, c))
            if P:
                supp.append(P)
    if verbose:
        print(f"d={d}: {len(supp)} support polynomials from {nblk} charge-0 blocks", flush=True)
    # ---- (S4a) support polynomials vanish on DKZ
    for P in supp:
        assert not eval_poly(K, P, R, d), "a support polynomial does not vanish on DKZ"
    # ---- targets
    T = targets(K, d)
    for name, P in T.items():
        assert not eval_poly(K, P, R, d), f"target {name} does not vanish on DKZ"
    if verbose:
        print(f"d={d}: (S4) all support polynomials and all {len(T)} targets vanish exactly on DKZ", flush=True)
    # ---- (S5) canonical model: R1c = X (|j> -> |j+1>), Dc = 1 - (1+i)|0><0|, R_{k+1}c = z Dc R_kc, and the
    #      monomial unitary W0 |j> = z^{-2j} |d-1-j> satisfies W0^* R_k^DKZ W0 = R_kc (k = 1..4), W0^* Pi W0 = |0><0|.
    Xc = ([(k + 1) % d for k in range(d)], [0] * d)
    zDc = (list(range(d)), [(1 + 3 * d) % N if k == 0 else 1 for k in range(d)])     # z * diag(-i, 1, ..., 1)
    Rc = [Xc]
    for _ in range(3):
        Rc.append(gp_mul(zDc, Rc[-1], N))
    W0 = ([d - 1 - j for j in range(d)], [(-2 * j) % N for j in range(d)])
    W0inv = ([0] * d, [0] * d)
    for j in range(d):
        W0inv[0][d - 1 - j] = j
        W0inv[1][d - 1 - j] = (2 * j) % N
    assert gp_mul(W0inv, W0, N) == (list(range(d)), [0] * d)
    for k in range(4):
        conj_k = gp_mul(gp_mul(W0inv, R[k], N), W0, N)
        assert conj_k[0] == Rc[k][0] and [x % N for x in conj_k[1]] == [x % N for x in Rc[k][1]], \
            f"W0^* R_{k+1} W0 != canonical"
        Pk = gp_mul(gp_pow(Rc[k], d - 1, d, N), Rc[k], N)
        assert Pk[0] == list(range(d)) and all(x % N == 0 for x in Pk[1]), "canonical R^d != 1"
    if verbose:
        print(f"d={d}: (S5) W0^* R_k^DKZ W0 = canonical R_k (k=1..4) exactly; canonical R_k^d = 1", flush=True)
    res = {}
    for name in sorted(T, key=lambda x: (x[0], x[1:])):
        c = solve_membership(K, supp, T[name])
        res[name] = c is not None
    okM = all(res[("M", n, t)] for n in range(1, d) for t in (1, 2, 3))
    okT = all(res[("T", n)] for n in range(1, d))
    okV = all(res[("V", n)] for n in range(1, d))
    if verbose:
        print(f"d={d}: (S2) chain modes M_(n,t), n=1..{d-1}, t=1,2,3 in support span: {okM}", flush=True)
        print(f"d={d}: (S3) symmetrised two-eigenvalue elements T_n, n=1..{d-1} in support span: {okT}", flush=True)
        print(f"d={d}: (S3') first-link elements V_n, n=1..{d-1} in support span: {okV}", flush=True)
        bad = [k for k, v in res.items() if not v]
        if bad:
            print(f"d={d}: FAILED targets: {bad}")
    assert okM and okT and okV
    # ---- negative controls: slightly wrong polynomials must be REJECTED (the membership test discriminates)
    iu = K.zp(d)
    neg = {}
    for n in range(1, d):
        P = {}                                            # chain mode t=1 without the phase w^-n on L^(3)
        for k in range(4):
            w, e = link(n, k, d)
            padd(K, P, wreduce(w, d), K.zp(-d * k))
        neg[("M-no-phase", n)] = P
        P = {}                                            # T_n with (1-i) replaced by (1+i)
        c1 = K.neg(K.mul(K.add(K.one, iu), K.zp(n)))
        cU = K.neg(K.mul(iu, K.zp(2 * n)))
        for k in range(4):
            wi, ei = link_inv(n, k, d)
            w, e = link(n, k, d)
            padd(K, P, wreduce(wi, d), K.zp(ei))
            padd(K, P, (), c1)
            padd(K, P, wreduce(w, d), K.mul(cU, K.zp(e)))
        neg[("T-wrong-sign", n)] = P
        P = dict(T[("V", n)])                             # V_n + 1
        padd(K, P, (), K.one)
        neg[("V-plus-1", n)] = P
    rejected = all(solve_membership(K, supp, P) is None for P in neg.values())
    if verbose:
        print(f"d={d}: negative controls ({len(neg)} perturbed polynomials) all rejected: {rejected}", flush=True)
    assert rejected
    m = d // 2 + 1
    print(f"d={d}: SUPPORT CONDITIONS VERIFIED (all n <= d-1; rigidity needs n <= floor(d/2)+1 = {m}) "
          f"[{time.time()-t0:.1f}s]", flush=True)
    return True


if __name__ == "__main__":
    args = sys.argv[1:]
    with_identity = "--identity" in args
    ds = [int(t) for t in args if not t.startswith("--")] or [3, 4, 5, 6, 7, 8, 9]
    for d in ds:
        check(d)
        if with_identity:
            check_identity(d)
    print("ALL DONE")
