"""
twirl_check.py -- numerical check of the rotation-twirl reduction (workstream B).

Tracial (maximally entangled) CGLMP, Gill form (to be MINIMISED):
  S(R) = 2(d-1) + sum_{n=1}^{d-1} c_n [tau(R1^n R2^-n) + tau(R2^n R3^-n) + tau(R3^n R4^-n) + w^-n tau(R4^n R1^-n)],
  c_n = -1/(1 - w^-n),  w = exp(2 pi i/d),  R_i^d = 1.
Rotation rho: (R1,R2,R3,R4) -> (R2,R3,R4,w R1) leaves S invariant; rho^4 = global phase w, so rho generates Z_{4d}.
Twirled strategy R' = direct sum over Z_{4d}; W = cyclic shift of the summands:
  W R'_i W^-1 = R'_{i+1} (i=1,2,3),  W R'_4 W^-1 = w R'_1,  W^{4d} = 1.
Then Z := W^4 satisfies Z R'_1 Z^-1 = w R'_1 (Weyl pair), so C^{D'} = C^d (x) C^M, R'_1 = X (x) 1, Z = Z0 (x) 1 and
W = sum_k |k><k| (x) W_k with W_k^4 = w^k.  Writing W_k = zeta^k V_k (zeta = exp(2 pi i/(4d)), k = 0..d-1),
V_k^4 = 1 and
  S = 2(d-1) - (2/d) F,   F = sum_{0<=j<k<=d-1} Re[ h_{k-j} tr(V_j V_k^*) ],   h_m = sec(psi_m) - i csc(psi_m),
  psi_m = pi m/(2d).
This script checks all of this numerically on random strategies.
"""
import numpy as np

rng = np.random.default_rng(7)


def rand_unitary(D):
    Z = rng.normal(size=(D, D)) + 1j * rng.normal(size=(D, D))
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / abs(np.diag(R)))


def rand_order_d(d, D, ranks=None):
    """random unitary with spectrum in d-th roots of unity (random ranks)."""
    U = rand_unitary(D)
    lab = rng.integers(0, d, size=D) if ranks is None else np.repeat(np.arange(d), ranks)
    w = np.exp(2j * np.pi / d)
    return U @ np.diag(w ** lab) @ U.conj().T


def S_value(R, d):
    w = np.exp(2j * np.pi / d)
    D = R[0].shape[0]
    tau = lambda M: np.trace(M) / D
    mp = np.linalg.matrix_power
    S = 2 * (d - 1)
    for n in range(1, d):
        c = -1 / (1 - w ** (-n))
        S += c * (tau(mp(R[0], n) @ mp(R[1], -n)) + tau(mp(R[1], n) @ mp(R[2], -n))
                  + tau(mp(R[2], n) @ mp(R[3], -n)) + w ** (-n) * tau(mp(R[3], n) @ mp(R[0], -n)))
    return S


def twirl(R, d):
    """direct sum over the Z_{4d} orbit of rho; returns R' and the shift W."""
    w = np.exp(2j * np.pi / d)
    D = R[0].shape[0]
    orbit = []
    cur = [r.copy() for r in R]
    for k in range(4 * d):
        orbit.append(cur)
        cur = [cur[1], cur[2], cur[3], w * cur[0]]
    N = 4 * d
    Rp = []
    for i in range(4):
        M = np.zeros((N * D, N * D), complex)
        for k in range(N):
            M[k * D:(k + 1) * D, k * D:(k + 1) * D] = orbit[k][i]
        Rp.append(M)
    # W maps summand k+1 -> summand k  (so that W R'_i W^-1 has block k equal to (rho^{k+1} R)_i = (rho^k R)_{i+1})
    Sh = np.zeros((N, N))
    for k in range(N):
        Sh[k, (k + 1) % N] = 1
    W = np.kron(Sh, np.eye(D))
    return Rp, W


def weyl_decompose(R1, Z, d):
    """find unitary basis change so that R1 = X (x) 1, Z = Z0 (x) 1.  Returns list of isometries E_k: C^M -> eigenspace
    of Z with eigenvalue w^k, with R1 E_k = E_{k+1} (X|k> = |k+1>)."""
    w = np.exp(2j * np.pi / d)
    ev, U = np.linalg.eig(Z)
    # Z unitary: use Schur/eigh of Hermitian parts to get orthonormal eigenspaces
    Dp = Z.shape[0]
    E = []
    # eigenspace of eigenvalue w^0 via projector (1/d) sum_j w^{0} Z^j
    P0 = sum(np.linalg.matrix_power(Z, j) for j in range(d)) / d
    u, s, vh = np.linalg.svd(P0)
    M = int(round(s.sum().real))
    E0 = u[:, :M]
    # Z R1 Z^-1 = w R1 -> R1 maps eigenvalue w^k of Z to eigenvalue w^{k-1}? check both
    E = [E0]
    for k in range(1, d):
        E.append(R1 @ E[-1])
    # determine the Z eigenvalue on each E_k
    lam = [np.vdot(E[k][:, 0], Z @ E[k][:, 0]) for k in range(d)]
    return E, lam


def reduced_F(V, d):
    F = 0.0
    for j in range(d):
        for k in range(j + 1, d):
            m = k - j
            psi = np.pi * m / (2 * d)
            h = 1 / np.cos(psi) - 1j / np.sin(psi)
            M = V[j].shape[0]
            F += (h * np.trace(V[j] @ V[k].conj().T) / M).real
    return F


if __name__ == "__main__":
    for d in (3, 4, 5):
        for trial in range(2):
            D = int(rng.integers(2, 5))
            R = [rand_order_d(d, D) for _ in range(4)]
            S0 = S_value(R, d)
            Rp, W = twirl(R, d)
            S1 = S_value(Rp, d)
            w = np.exp(2j * np.pi / d)
            errs = [np.abs(W @ Rp[i] @ W.conj().T - Rp[i + 1]).max() for i in range(3)]
            errs.append(np.abs(W @ Rp[3] @ W.conj().T - w * Rp[0]).max())
            Z = np.linalg.matrix_power(W, 4)
            eW = np.abs(np.linalg.matrix_power(W, 4 * d) - np.eye(W.shape[0])).max()
            eZ = np.abs(Z @ Rp[0] @ Z.conj().T - w * Rp[0]).max()
            E, lam = weyl_decompose(Rp[0], Z, d)
            # lam[k] should be w^{a_k}; find labelling
            ks = [int(round(np.angle(l) / (2 * np.pi / d))) % d for l in lam]
            # W_k blocks: W restricted to eigenspace of Z with eigenvalue w^k
            zeta = np.exp(2j * np.pi / (4 * d))
            Vs = [None] * d
            for idx in range(d):
                k = ks[idx]
                Wk = E[idx].conj().T @ W @ E[idx]
                Vs[k] = Wk / zeta ** k
            e4 = max(np.abs(np.linalg.matrix_power(v, 4) - np.eye(v.shape[0])).max() for v in Vs)
            F = reduced_F(Vs, d)
            S2 = 2 * (d - 1) - 2 * F / d
            print(f"d={d} D={D}: S={S0.real:.12f} (im {abs(S0.imag):.1e})  S_twirl={S1.real:.12f}  "
                  f"covariance err={max(errs):.1e} W^4d err={eW:.1e} Weyl err={eZ:.1e}  "
                  f"Z labels={ks}  V^4 err={e4:.1e}  S_reduced={S2:.12f}")
    # DKZ
    import sys
    sys.path.insert(0, r"<project>/iqoqi\programs\oqp27B")
    from dkz_tracial import unitaries
    for d in (3, 4, 5, 6):
        R = unitaries(d)
        k = np.arange(d)
        Kd = 1 / np.cos(np.pi * (k[:, None] - k[None, :]) / (2 * d))
        Sme = 2 * d - 1 - Kd.sum() / d
        Fdkz = reduced_F([np.eye(1)] * d, d)
        print(f"DKZ d={d}: S={S_value(R, d).real:.12f}  S_ME={Sme:.12f}  from V=1: {2*(d-1)-2*Fdkz/d:.12f}")
