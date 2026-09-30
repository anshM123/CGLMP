"""
extfield.py -- simple algebraic extension L = F[t]/(p(t)) of a cyclotomic field F = Q(zeta_N)
(cyclo.Field), p monic with coefficients in the real subfield of F (so complex conjugation acts
coefficient-wise).  t is a chosen real root of p (given numerically for evaluation).
Elements: tuples of m F-elements (coefficients of 1, t, ..., t^{m-1}).
"""
import mpmath as mp


class Ext:
    def __init__(self, F, p, troot):
        """p: list of F elements, low->high, monic (p[-1] == 1). troot: mpmath real root."""
        self.F = F
        self.p = p
        self.m = len(p) - 1
        self.troot = troot
        assert F.is_zero(F.sub(p[-1], F.one))
        for c in p:
            assert F.is_real(c), "p must have real coefficients"
        self.zero = tuple([F.zero] * self.m)
        self.one = tuple([F.one] + [F.zero] * (self.m - 1))

    def from_F(self, a):
        return tuple([a] + [self.F.zero] * (self.m - 1))

    def t(self):
        assert self.m >= 2
        return tuple([self.F.zero, self.F.one] + [self.F.zero] * (self.m - 2))

    def is_zero(self, x):
        return all(self.F.is_zero(c) for c in x)

    def add(self, x, y):
        return tuple(self.F.add(a, b) for a, b in zip(x, y))

    def sub(self, x, y):
        return tuple(self.F.sub(a, b) for a, b in zip(x, y))

    def neg(self, x):
        return tuple(self.F.neg(a) for a in x)

    def mulF(self, x, a):
        return tuple(self.F.mul(c, a) for c in x)

    def mul(self, x, y):
        F, m = self.F, self.m
        prod = [F.zero] * (2 * m - 1)
        for i, a in enumerate(x):
            if F.is_zero(a):
                continue
            for j, b in enumerate(y):
                if F.is_zero(b):
                    continue
                prod[i + j] = F.add(prod[i + j], F.mul(a, b))
        # reduce: t^m = -sum_{k<m} p_k t^k
        for k in range(2 * m - 2, m - 1, -1):
            c = prod[k]
            if F.is_zero(c):
                continue
            prod[k] = F.zero
            for j in range(m):
                if not F.is_zero(self.p[j]):
                    prod[k - m + j] = F.sub(prod[k - m + j], F.mul(c, self.p[j]))
        return tuple(prod[:m])

    def conj(self, x):
        return tuple(self.F.conj(a) for a in x)

    def inv(self, x):
        F, m = self.F, self.m
        assert not self.is_zero(x)
        # multiplication matrix columns: x * t^k
        cols = []
        e = self.one
        for k in range(m):
            cols.append(self.mul(x, e))
            e = self.mul(e, self.t()) if m >= 2 else e
        A = [[cols[j][i] for j in range(m)] + [F.one if i == 0 else F.zero] for i in range(m)]
        for c in range(m):
            p = next(r for r in range(c, m) if not F.is_zero(A[r][c]))
            A[c], A[p] = A[p], A[c]
            iv = F.inv(A[c][c])
            A[c] = [F.mul(v, iv) for v in A[c]]
            for r in range(m):
                if r != c and not F.is_zero(A[r][c]):
                    f = A[r][c]
                    A[r] = [F.sub(vr, F.mul(f, vc)) for vr, vc in zip(A[r], A[c])]
        return tuple(A[i][m] for i in range(m))

    def div(self, x, y):
        return self.mul(x, self.inv(y))

    def re(self, x):
        return tuple(self.F.re(a) for a in x)

    def im(self, x):
        return tuple(self.F.im(a) for a in x)

    def is_real(self, x):
        return all(self.F.is_real(a) for a in x)

    def to_complex(self, x, dps=50):
        with mp.workdps(dps):
            tt = mp.mpf(self.troot)
            s = mp.mpc(0)
            pw = mp.mpf(1)
            for a in x:
                s += self.F.to_complex(a, dps) * pw
                pw *= tt
            return s

    def to_cfloat(self, x):
        return complex(self.to_complex(x, 30))

    # --- interface helpers used by the certificate builders
    def from_frac(self, fr):
        return self.from_F(self.F.from_frac(fr))

    def i_unit(self):
        return self.from_F(self.F.i_unit())

    def one_(self):
        return self.one

    def from_int(self, a):
        return self.from_F(self.F.from_int(a))
