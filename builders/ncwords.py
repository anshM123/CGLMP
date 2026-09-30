"""
ncwords.py -- words in free products of cyclic groups Z_d (unitaries U_g with U_g^d = 1),
moment-matrix construction in the Fourier (unitary-power) basis, and a Clarabel SDP wrapper.

Letter = (g, n) with generator g and power n in 1..d-1.  Word = tuple of letters, reduced
(no two adjacent letters with the same generator).
"""
import numpy as np
import scipy.sparse as sp


def reduce_word(w, d):
    out = []
    for (g, n) in w:
        n %= d
        if n == 0:
            continue
        if out and out[-1][0] == g:
            m = (out[-1][1] + n) % d
            out.pop()
            if m != 0:
                out.append((g, m))
        else:
            out.append((g, n))
    return tuple(out)


def dagger(w, d):
    return tuple((g, (-n) % d) for (g, n) in reversed(w))


def charge(w, d, weights=None):
    return sum(n * (1 if weights is None else weights[g]) for g, n in w) % d


def cyclic_canon(w, d):
    """canonical representative of a word under cyclic rotation (for traces)."""
    w = list(reduce_word(w, d))
    # cyclic reduction: merge first and last letter if same generator
    while len(w) >= 2 and w[0][0] == w[-1][0]:
        g = w[0][0]
        m = (w[0][1] + w[-1][1]) % d
        w = w[1:-1]
        if m != 0:
            w = [(g, m)] + w
            w = list(reduce_word(tuple(w), d))
    w = tuple(w)
    if len(w) <= 1:
        return w
    rots = [w[i:] + w[:i] for i in range(len(w))]
    return min(rots)


class SDP:
    """Collects Hermitian PSD blocks whose entries are affine in real variables.
    Each block: dict (i,j) -> list of (coef complex, var index or -1 for constant)."""

    def __init__(self):
        self.nvar = 0
        self.blocks = []
        self.eqs = []   # list of (dict var->coef, rhs) equality constraints on real vars

    def new_var(self):
        self.nvar += 1
        return self.nvar - 1

    def add_block(self, n, entry_fn):
        """entry_fn(i,j) -> list of (complex coef, var or -1); only i<=j requested."""
        ent = {}
        for i in range(n):
            for j in range(i, n):
                ent[(i, j)] = entry_fn(i, j)
        self.blocks.append((n, ent))

    def solve(self, cvec, const=0.0, tol=1e-10, verbose=False, maxiter=400):
        """minimise cvec . y + const subject to blocks >= 0 (Hermitian, embedded as real 2n)."""
        import clarabel
        rows, cols, vals = [], [], []
        b = []
        cones = []
        r0 = 0
        # equality constraints first
        if self.eqs:
            for (coefs, rhs) in self.eqs:
                for v, c in coefs.items():
                    rows.append(r0); cols.append(v); vals.append(c)
                b.append(rhs)
                r0 += 1
            cones.append(clarabel.ZeroConeT(len(self.eqs)))
        for (n, ent) in self.blocks:
            # detect whether block is real
            is_real = all(abs(np.imag(c)) < 1e-15 for lst in ent.values() for (c, v) in lst)
            N = n if is_real else 2 * n

            def real_entry(I, J):
                # returns list of (real coef, var) for the real embedding entry (I,J), I<=J
                if is_real:
                    i, j = I, J
                    lst = ent[(i, j)] if i <= j else ent[(j, i)]
                    return [(np.real(c), v) for (c, v) in lst]
                bi, i = divmod(I, n)
                bj, j = divmod(J, n)
                # H = X + iY ; embedding [[X,-Y],[Y,X]]
                if i <= j:
                    lst = ent[(i, j)]
                    Hre = [(np.real(c), v) for (c, v) in lst]
                    Him = [(np.imag(c), v) for (c, v) in lst]
                else:
                    lst = ent[(j, i)]
                    Hre = [(np.real(c), v) for (c, v) in lst]
                    Him = [(-np.imag(c), v) for (c, v) in lst]
                if bi == bj:
                    return Hre
                if bi == 0 and bj == 1:
                    return [(-c, v) for (c, v) in Him]
                return Him
            # column-major upper triangle with sqrt2 scaling
            for J in range(N):
                for I in range(J + 1):
                    sc = 1.0 if I == J else np.sqrt(2.0)
                    lst = real_entry(I, J)
                    cst = 0.0
                    for (c, v) in lst:
                        if c == 0:
                            continue
                        if v < 0:
                            cst += c
                        else:
                            rows.append(r0); cols.append(v); vals.append(-sc * c)
                    b.append(sc * cst)
                    r0 += 1
            cones.append(clarabel.PSDTriangleConeT(N))
        A = sp.csc_matrix((vals, (rows, cols)), shape=(r0, self.nvar))
        P = sp.csc_matrix((self.nvar, self.nvar))
        q = np.array(cvec, float)
        settings = clarabel.DefaultSettings()
        settings.verbose = verbose
        settings.tol_gap_abs = tol
        settings.tol_gap_rel = tol
        settings.tol_feas = tol
        settings.tol_ktratio = 1e-8
        settings.max_iter = maxiter
        solver = clarabel.DefaultSolver(P, q, A, np.array(b), cones, settings)
        sol = solver.solve()
        # record dual PSD matrices (complex Hermitian form) and primal blocks
        z = np.array(sol.z)
        s = np.array(sol.s)
        off = len(self.eqs)
        self.dual_blocks, self.primal_blocks = [], []
        for (n, ent) in self.blocks:
            is_real = all(abs(np.imag(c)) < 1e-15 for lst in ent.values() for (c, v) in lst)
            N = n if is_real else 2 * n
            m = N * (N + 1) // 2
            Zr = self._smat(z[off:off + m], N)
            Sr = self._smat(s[off:off + m], N)
            off += m
            if is_real:
                self.dual_blocks.append(Zr.astype(complex))
                self.primal_blocks.append(Sr.astype(complex))
            else:
                # real embedding [[X,-Y],[Y,X]] of H = X + iY; the dual of the embedded cone
                # corresponds to Zc = (Z11+Z22) + i(Z21-Z12)  (up to the factor fixed by <Z,M> pairing)
                Z11, Z12, Z21, Z22 = Zr[:n, :n], Zr[:n, n:], Zr[n:, :n], Zr[n:, n:]
                self.dual_blocks.append((Z11 + Z22) + 1j * (Z21 - Z12))
                self.primal_blocks.append(Sr[:n, :n] + 1j * Sr[n:, :n])
        return sol.obj_val + const, np.array(sol.x), sol

    @staticmethod
    def _smat(v, N):
        M = np.zeros((N, N))
        k = 0
        for J in range(N):
            for I in range(J + 1):
                if I == J:
                    M[I, J] = v[k]
                else:
                    M[I, J] = M[J, I] = v[k] / np.sqrt(2.0)
                k += 1
        return M
