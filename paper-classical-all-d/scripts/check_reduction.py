"""Numerical check of the physics side of the paper (Section 2 and Appendix A): the reduction of the CGLMP
expression on maximally entangled states to the Z_4 clock model, and the equal-link (commutative) sector.

All CGLMP values below are computed from outcome probabilities with the expression of Collins, Gisin, Linden,
Massar and Popescu (PRL 88, 040404 (2002)); nothing is taken from the clock-model formulas unless stated.

 (R1) for arbitrary behaviours P(a,b|x,y) (normalisation only):  I_d = 4 - 2 S/(d-1), with the chain labelling
      X1 = A2, Y1 = B2, X2 = A1, Y2 = B1 and S = E m(X1-Y1) + E m(Y1-X2) + E m(X2-Y2) + E m(Y2-X1-1), m(t) = t mod d.
      Here "P(A_a = B_b + k)" is read as the probability that A_a - B_b = k mod d.
 (R2) random maximally entangled PVM strategies (any D, random ranks): S equals the trace formula
      2(d-1) + sum_n c_n [tau(R1^n R2^-n) + tau(R2^n R3^-n) + tau(R3^n R4^-n) + w^-n tau(R4^n R1^-n)],
      and (no-signalling behaviours) S = d Pi - 1, Pi = P(X1<Y1) + P(Y1<X2) + P(X2<Y2) + P(Y2<=X1).
 (R3) covariant strategies built from random NON-commuting order-4 unitaries V_k:  I_d = 4 F(V)/(d(d-1)).
 (R4) twirl + Stone-von Neumann applied to random strategies returns order-4 unitaries V_k with I_d = 4F(V)/(d(d-1)).
 (R5) equal links <=> commuting V: a random strategy has unequal links and non-commuting V; an equal-link strategy
      (commuting V, then a random local unitary) has equal links and its twirl gives commuting V.
 (R6) clock-phase strategies W = diag(z^k i^{-a_k}), all a in Z_4^d (d <= 7): I_d = 4F(a)/(d(d-1)); the maximum is
      I_ME(d), attained exactly at the 4d one-step configurations.
 (R7) the CGLMP/DKZ measurements of Collins et al. (their eqs. (12)-(15)) attain I_ME(d), and in chain order they
      equal G R^(0)_i G^* with G = W_0^{-2}, where R^(0) is the covariant strategy with V = 1.
 (R8) projection form (discussion section): F = (1/2) sum_{x,y in Z_4d} cot(pi(x-y)/4d) tau(A_x B_y) - d/2 with
      A_x = Q_x, B_x = Q_{x+d} for random non-commuting V.
 (R9) Ising form (discussion section): F = (1/2) sum_{n in Z_2d} sum_{u=1}^{d-1} sec(pi u/2d) tau(X_n X_{n+u}) with
      E_k = e^{i pi/4}(X_k - i X_{k+d})/sqrt2, X_{n+2d} = -X_n, for random non-commuting V; and the value d(d-1) of the
      Clifford rotation X_n = cos(pi n/2d) s1 + sin(pi n/2d) s2 (which violates [X_n, X_{n+d}] = 0), compared with F_DKZ.
"""
import itertools
import numpy as np

rng = np.random.default_rng(2026)
mp = np.linalg.matrix_power


def cglmp_standard(P):
    """P[x, y, a, b], x, y in {0, 1} (settings 1, 2).  Collins et al. (2002) CGLMP expression."""
    d = P.shape[2]

    def pd(x, y, t):  # P((A_x - B_y) mod d = t)
        return sum(P[x, y, (j + t) % d, j] for j in range(d))
    I = 0.0
    for k in range(d // 2):
        c = 1 - 2 * k / (d - 1)
        I += c * (pd(0, 0, k) + pd(1, 0, -(k + 1)) + pd(1, 1, k) + pd(0, 1, -k)
                  - pd(0, 0, -k - 1) - pd(1, 0, k) - pd(1, 1, -k - 1) - pd(0, 1, k + 1))
    return I


def chain_S_Pi(P):
    d = P.shape[2]
    a = np.arange(d)
    # E m(A_x - B_y) etc.; chain labelling X1 = A2 (x=1), Y1 = B2 (y=1), X2 = A1 (x=0), Y2 = B1 (y=0)
    m_ab = lambda x, y: sum(P[x, y, i, j] * ((i - j) % d) for i in a for j in a)            # E m(A_x - B_y)
    m_ba = lambda x, y, s: sum(P[x, y, i, j] * ((j - i - s) % d) for i in a for j in a)     # E m(B_y - A_x - s)
    S = m_ab(1, 1) + m_ba(0, 1, 0) + m_ab(0, 0) + m_ba(1, 0, 1)
    lt = lambda x, y: sum(P[x, y, i, j] for i in a for j in a if i < j)    # P(A_x < B_y)
    gt = lambda x, y: sum(P[x, y, i, j] for i in a for j in a if j < i)    # P(B_y < A_x)
    ge = lambda x, y: sum(P[x, y, i, j] for i in a for j in a if j <= i)   # P(B_y <= A_x)
    Pi = lt(1, 1) + gt(0, 1) + lt(0, 0) + ge(1, 0)
    return S, Pi


def rand_unitary(n):
    Z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / abs(np.diag(R)))


def rand_pvm(d, D):
    U = rand_unitary(D)
    lab = rng.integers(0, d, size=D)
    return [U[:, lab == a] @ U[:, lab == a].conj().T for a in range(d)]


def spectral_pvm(R, d):
    """spectral projections of an order-d unitary: P_a = (1/d) sum_n w^{-na} R^n (eigenvalue w^a)."""
    w = np.exp(2j * np.pi / d)
    pw = [mp(R, n) for n in range(d)]
    return [sum(w ** (-n * a) * pw[n] for n in range(d)) / d for a in range(d)]


def probs_from_chain(R, d):
    """R = (R1, R2, R3, R4) = (A2, B2^T, A1, B1^T) as order-d unitaries; returns P[x, y, a, b] on |Phi_D>."""
    D = R[0].shape[0]
    A = {1: spectral_pvm(R[0], d), 0: spectral_pvm(R[2], d)}     # Alice setting index 1 = A2, 0 = A1
    Bt = {1: spectral_pvm(R[1], d), 0: spectral_pvm(R[3], d)}    # transposed Bob projections
    P = np.zeros((2, 2, d, d))
    for x in (0, 1):
        for y in (0, 1):
            for a in range(d):
                for b in range(d):
                    P[x, y, a, b] = (np.trace(A[x][a] @ Bt[y][b]) / D).real
    return P


def S_trace(R, d):
    w = np.exp(2j * np.pi / d)
    D = R[0].shape[0]
    tau = lambda M: np.trace(M) / D
    S = 2 * (d - 1)
    for n in range(1, d):
        c = -1 / (1 - w ** (-n))
        S += c * (tau(mp(R[0], n) @ mp(R[1], -n)) + tau(mp(R[1], n) @ mp(R[2], -n))
                  + tau(mp(R[2], n) @ mp(R[3], -n)) + w ** (-n) * tau(mp(R[3], n) @ mp(R[0], -n)))
    return S


def F_clock(V, d):
    F = 0.0
    for j in range(d):
        for k in range(j + 1, d):
            psi = np.pi * (k - j) / (2 * d)
            h = 1 / np.cos(psi) - 1j / np.sin(psi)
            F += (h * np.trace(V[j] @ V[k].conj().T) / V[j].shape[0]).real
    return F


def F_DKZ(d):
    return sum((d - m) / np.cos(np.pi * m / (2 * d)) for m in range(1, d))


def shift(d):
    X = np.zeros((d, d))
    for k in range(d):
        X[(k + 1) % d, k] = 1          # X|k> = |k+1>
    return X


def covariant(V, d):
    M = V[0].shape[0]
    z = np.exp(2j * np.pi / (4 * d))
    W = np.zeros((d * M, d * M), complex)
    for k in range(d):
        W[k * M:(k + 1) * M, k * M:(k + 1) * M] = z ** k * V[k]
    R1 = np.kron(shift(d), np.eye(M))
    R = [R1]
    for i in range(3):
        R.append(W @ R[-1] @ W.conj().T)
    return R, W


def rand_order4(M):
    """random unitary with V^4 = 1 in a random basis; for M >= 2 at least two distinct eigenvalues (not scalar)."""
    U = rand_unitary(M)
    r = rng.integers(0, 4, size=M)
    if M >= 2 and len(set(r.tolist())) == 1:
        r[1] = (r[0] + 1 + rng.integers(0, 3)) % 4
    return U @ np.diag(1j ** r) @ U.conj().T


def links(R, d, n):
    w = np.exp(2j * np.pi / d)
    return [mp(R[0], n) @ mp(R[1], -n), mp(R[1], n) @ mp(R[2], -n), mp(R[2], n) @ mp(R[3], -n),
            w ** (-n) * mp(R[3], n) @ mp(R[0], -n)]


def twirl_reduce(R, d):
    """rotation twirl + Stone-von Neumann; returns the order-4 unitaries V_k (Appendix A)."""
    w = np.exp(2j * np.pi / d)
    z = np.exp(2j * np.pi / (4 * d))
    D = R[0].shape[0]
    N = 4 * d
    orbit, cur = [], [r.copy() for r in R]
    for r in range(N):
        orbit.append(cur)
        cur = [cur[1], cur[2], cur[3], w * cur[0]]
    Rt = []
    for i in range(4):
        Mx = np.zeros((N * D, N * D), complex)
        for r in range(N):
            Mx[r * D:(r + 1) * D, r * D:(r + 1) * D] = orbit[r][i]
        Rt.append(Mx)
    Sh = np.zeros((N, N))
    for r in range(N):
        Sh[r, (r + 1) % N] = 1                              # W (e_{r+1} (x) v) = e_r (x) v
    W = np.kron(Sh, np.eye(D))
    cov = max(np.abs(W @ Rt[i] @ W.conj().T - Rt[i + 1]).max() for i in range(3))
    cov = max(cov, np.abs(W @ Rt[3] @ W.conj().T - w * Rt[0]).max())
    Z = mp(W, 4)
    P0 = sum(mp(Z, j) for j in range(d)) / d                # projection onto ker(Z - 1)
    u, s, _ = np.linalg.svd(P0)
    M = int(round(s.sum()))
    E = [u[:, :M]]
    for k in range(1, d):
        E.append(Rt[0] @ E[-1])                             # e_{k,mu} = R~_1^k f_mu, in ker(Z - w^k)
    zerr = max(np.abs(Z @ E[k] - w ** k * E[k]).max() for k in range(d))
    V = [z ** (-k) * (E[k].conj().T @ W @ E[k]) for k in range(d)]
    return V, cov, zerr


def comm_norm(V):
    return max(np.abs(V[j] @ V[k] - V[k] @ V[j]).max() for j in range(len(V)) for k in range(len(V)))


if __name__ == "__main__":
    # (R1)
    err = 0.0
    for d in (2, 3, 4, 5, 6, 7):
        for t in range(20):
            P = rng.random((2, 2, d, d))                     # arbitrary (even signalling) behaviours
            P /= P.sum(axis=(2, 3), keepdims=True)
            S, Pi = chain_S_Pi(P)
            err = max(err, abs(cglmp_standard(P) - (4 - 2 * S / (d - 1))))
    print(f"(R1) random behaviours, d = 2..7: max |I_CGLMP - (4 - 2S/(d-1))| = {err:.1e}")
    # (R2)
    err = 0.0
    for d in (3, 4, 5, 6):
        for t in range(6):
            D = int(rng.integers(2, 7))
            pv = [rand_pvm(d, D) for _ in range(4)]           # A1, A2, B1, B2
            w = np.exp(2j * np.pi / d)
            U = [sum(w ** a * p[a] for a in range(d)) for p in pv]
            R = [U[1], U[3].T, U[0], U[2].T]                   # (A2, B2^T, A1, B1^T)
            P = probs_from_chain(R, d)
            I = cglmp_standard(P)
            S = S_trace(R, d)
            Sp, Pi = chain_S_Pi(P)
            err = max(err, abs(I - (4 - 2 * S.real / (d - 1))), abs(S.imag), abs(Sp - (d * Pi - 1)))
    print(f"(R2) random max-ent PVM strategies (d = 3..6, D = 2..6): max |I_CGLMP - (4 - 2 S_trace/(d-1))|, "
          f"|S - (d Pi - 1)| = {err:.1e}")
    # (R3)
    err = 0.0
    cmin = 1e9
    for d in (2, 3, 4, 5, 6):
        for M in (1, 2, 3):
            V = [rand_order4(M) for _ in range(d)]
            cmin = min(cmin, comm_norm(V)) if M > 1 else cmin
            R, W = covariant(V, d)
            P = probs_from_chain(R, d)
            err = max(err, abs(cglmp_standard(P) - 4 * F_clock(V, d) / (d * (d - 1))))
    print(f"(R3) covariant strategies from random order-4 V (M = 1..3, d = 2..6; min commutator norm {cmin:.2f}): "
          f"max |I_CGLMP - 4F(V)/(d(d-1))| = {err:.1e}")
    # (R4), (R5) part 1
    err = e4 = covmax = zmax = 0.0
    lk_min = cm_min = 1e9
    for d in (3, 4, 5):
        for t in range(3):
            D = int(rng.integers(2, 4))
            pv = [rand_pvm(d, D) for _ in range(4)]
            w = np.exp(2j * np.pi / d)
            U = [sum(w ** a * p[a] for a in range(d)) for p in pv]
            R = [U[1], U[3].T, U[0], U[2].T]
            I = cglmp_standard(probs_from_chain(R, d))
            V, cov, zerr = twirl_reduce(R, d)
            covmax, zmax = max(covmax, cov), max(zmax, zerr)
            e4 = max(e4, max(np.abs(mp(v, 4) - np.eye(v.shape[0])).max() for v in V))
            err = max(err, abs(I - 4 * F_clock(V, d) / (d * (d - 1))))
            L = links(R, d, 1)
            lk_min = min(lk_min, max(np.abs(L[i] - L[0]).max() for i in range(1, 4)))
            cm_min = min(cm_min, comm_norm(V))
    print(f"(R4) twirl + Stone-von Neumann on random strategies (d = 3..5): covariance err {covmax:.1e}, "
          f"Z-eigenspace err {zmax:.1e}, max |V^4 - 1| = {e4:.1e}, max |I_CGLMP(R) - 4F(V)/(d(d-1))| = {err:.1e}")
    print(f"(R5a) same random strategies: min over samples of max|L_1^(i) - L_1^(0)| = {lk_min:.2f} (links unequal), "
          f"min commutator norm of the V's = {cm_min:.2f} (V non-commuting)")
    # (R5) part 2: equal-link strategies
    lk = cm = err = 0.0
    for d in (3, 4, 5):
        M = 2
        a = rng.integers(0, 4, size=(M, d))
        U0 = rand_unitary(M)
        V = [U0 @ np.diag(1j ** (-a[:, k])) @ U0.conj().T for k in range(d)]   # commuting order-4 unitaries
        R, W = covariant(V, d)
        G = rand_unitary(d * M)
        R = [G @ r @ G.conj().T for r in R]                  # random local unitary G (x) conj(G)
        for n in range(1, d):
            L = links(R, d, n)
            lk = max(lk, max(np.abs(L[i] - L[0]).max() for i in range(1, 4)))
        V2, cov, zerr = twirl_reduce(R, d)
        cm = max(cm, comm_norm(V2))
        err = max(err, abs(cglmp_standard(probs_from_chain(R, d)) - 4 * F_clock(V2, d) / (d * (d - 1))))
    print(f"(R5b) equal-link strategies (commuting V, random local unitary, d = 3..5): max link difference {lk:.1e}, "
          f"commutator norm of the recovered V's {cm:.1e}, value err {err:.1e}")
    # (R6)
    for d in range(2, 8):
        IME = 4 * F_DKZ(d) / (d * (d - 1))
        best, arg, err = -1e9, [], 0.0
        onestep = {tuple((c + (1 if k >= r else 0)) % 4 for k in range(d)) for c in range(4) for r in range(d)}
        z = np.exp(2j * np.pi / (4 * d))
        for a in itertools.product(range(4), repeat=d):
            V = [np.array([[1j ** (-a[k])]]) for k in range(d)]
            if d <= 5 or a[0] == 0:                          # probabilities for all a when d <= 5
                R, _ = covariant(V, d)
                I = cglmp_standard(probs_from_chain(R, d))
                err = max(err, abs(I - 4 * F_clock(V, d) / (d * (d - 1))))
            else:
                I = 4 * F_clock(V, d) / (d * (d - 1))
            if I > best + 1e-9:
                best, arg = I, [a]
            elif I > best - 1e-9:
                arg.append(a)
        if d > 5:   # global phase a -> a + c maps I to itself; complete the argmax list
            arg = {tuple((x + c) % 4 for x in t) for t in arg for c in range(4)}
        print(f"(R6) d={d}: max over a in Z_4^d of I_d = {best:.12f}, I_ME(d) = {IME:.12f}, #maximisers = {len(set(arg))}"
              f" (4d = {4*d}), all one-step: {set(arg) <= onestep}; max |I_CGLMP - 4F(a)/(d(d-1))| = {err:.1e}")
    # (R7)  CGLMP measurements exactly as described in Collins et al. (2002), eqs. (12)-(15):
    #       phases exp(2 pi i alpha_a j/d), exp(2 pi i beta_b j/d), then a DFT on each side, then computational basis.
    for d in range(2, 9):
        w = np.exp(2j * np.pi / d)
        j = np.arange(d)
        al, be = (0.0, 0.5), (0.25, -0.25)
        # measurement vectors |k>_{A,a} = d^{-1/2} sum_j exp(-2 pi i j (k + alpha_a)/d)|j>,
        #                     |l>_{B,b} = d^{-1/2} sum_j exp(+2 pi i j (l - beta_b)/d)|j>
        psiA = [np.array([np.exp(-2j * np.pi * j * (k + al[x]) / d) / np.sqrt(d) for k in range(d)]).T for x in (0, 1)]
        psiB = [np.array([np.exp(2j * np.pi * j * (l - be[y]) / d) / np.sqrt(d) for l in range(d)]).T for y in (0, 1)]
        P = np.zeros((2, 2, d, d))
        ferr = 0.0
        for x in (0, 1):
            for y in (0, 1):
                for k in range(d):
                    for l in range(d):
                        amp = np.sum(psiA[x][:, k].conj() * psiB[y][:, l].conj()) / np.sqrt(d)
                        P[x, y, k, l] = abs(amp) ** 2
                        sarg = np.sin(np.pi * (k - l + al[x] + be[y]) / d)
                        ferr = max(ferr, abs(P[x, y, k, l] - 1 / (2 * d ** 3 * sarg ** 2)))   # eq. (15)
        I = cglmp_standard(P)
        UA = [psiA[x] @ np.diag(w ** j) @ psiA[x].conj().T for x in (0, 1)]
        UB = [psiB[y] @ np.diag(w ** j) @ psiB[y].conj().T for y in (0, 1)]
        Rstd = [UA[1], UB[1].T, UA[0], UB[0].T]
        R0, _ = covariant([np.eye(1)] * d, d)
        G = np.diag(np.exp(2j * np.pi * j / (4 * d)) ** (-2))          # G = W_0^{-2}
        diff = max(np.abs(Rstd[i] - G @ R0[i] @ G.conj().T).max() for i in range(4))
        print(f"(R7) d={d}: CGLMP value of the CGLMP/DKZ measurements = {I:.12f}, I_ME(d) = {4*F_DKZ(d)/(d*(d-1)):.12f} "
              f"(probabilities vs eq. (15): err {ferr:.1e}); max |R_i(DKZ) - W0^-2 R^(0)_i W0^2| = {diff:.1e}")
    # (R8)
    err = 0.0
    for d in (2, 3, 4, 5):
        N = 4 * d
        for M in (2, 3):
            E = [rand_order4(M) for _ in range(d)]      # E_k = V_k^*
            V = [e.conj().T for e in E]
            Q = {}
            for k in range(d):
                for a in range(4):
                    Q[(k - d * a) % N] = sum((1j ** (-a * n)) * mp(E[k], n) for n in range(4)) / 4  # eigenvalue i^a
            k_ = lambda u: 0.0 if u % N == 0 else 1 / np.tan(np.pi * (u % N) / N)
            val = sum(k_(x - y) * np.trace(Q[x] @ Q[(y + d) % N]).real / M for x in range(N) for y in range(N))
            err = max(err, abs(F_clock(V, d) - (val / 2 - d / 2)))
    print(f"(R8) projection form F = (1/2) sum cot(pi(x-y)/4d) tau(A_x B_y) - d/2, random non-commuting E (d = 2..5): "
          f"max err {err:.1e}")

    # (R9)
    err = 0.0
    for d in (2, 3, 4, 5, 6):
        for M in (2, 3):
            E = [rand_order4(M) for _ in range(d)]
            V = [e.conj().T for e in E]
            P = []
            for k in range(d):
                P.append([sum((1j ** (-a * n)) * mp(E[k], n) for n in range(4)) / 4 for a in range(4)])  # eigenvalue i^a
            X = {}
            for k in range(d):
                X[k] = P[k][0] + P[k][1] - P[k][2] - P[k][3]
                X[k + d] = P[k][0] - P[k][1] - P[k][2] + P[k][3]
            recon = max(np.abs(np.exp(1j * np.pi / 4) * (X[k] - 1j * X[k + d]) / np.sqrt(2) - E[k]).max() for k in range(d))
            Xn = lambda n: X[n % (2 * d)] * (1 if (n // (2 * d)) % 2 == 0 else -1)
            val = 0.5 * sum(1 / np.cos(np.pi * u / (2 * d)) * np.trace(Xn(n) @ Xn(n + u)).real / M
                            for n in range(2 * d) for u in range(1, d))
            err = max(err, abs(F_clock(V, d) - val), recon)
    s1 = np.array([[0, 1], [1, 0]]); s2 = np.array([[0, -1j], [1j, 0]])
    cliff = []
    for d in (2, 3, 5, 8, 13):
        Xc = lambda n: np.cos(np.pi * n / (2 * d)) * s1 + np.sin(np.pi * n / (2 * d)) * s2
        val = 0.5 * sum(1 / np.cos(np.pi * u / (2 * d)) * np.trace(Xc(n) @ Xc(n + u)).real / 2
                        for n in range(2 * d) for u in range(1, d))
        cliff.append((d, float(round(val, 10)), d * (d - 1), float(round(F_DKZ(d), 6))))
    print(f"(R9) Ising form, random non-commuting E (d = 2..6): max err {err:.1e};  Clifford rotation "
          f"(d, value, d(d-1), F_DKZ): {cliff}")
