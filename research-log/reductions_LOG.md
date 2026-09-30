# Workstream B (symmetry / operator / representation theory) -- log

Folder: iqoqi/programs/oqp27B_all/B_operator/.  Python: python.
Conventions as in ../../oqp27B/LOG.md sections 2-4 and the shared board: tau = normalised trace, R1=A1, R2=B1^T, R3=A2,
R4=B2^T (R_i^d = 1), Gill form S = d*Pi - 1 (minimise), I_d = 4 - 2S/(d-1).  w = exp(2 pi i/d), z = exp(2 pi i/(4d)),
psi_m = pi m/(2d).

## 1. Rotation twirl + Stone-von Neumann: EXACT reduction to a Z_4 clock model (PROVED, all d)  [2026-09-29 ~16:00]
* rho(R1,R2,R3,R4) = (R2,R3,R4,wR1) and the global phase R -> wR leave S invariant; rho^4 = global phase w; so rho
  generates Z_{4d}.  Twirl: R~ = direct sum over the Z_{4d} orbit; W = cyclic shift of the summands:
  W R~_i W^* = R~_{i+1}, W R~_4 W^* = w R~_1, W^{4d} = 1, S(R~) = S(R).
* Z := W^4 satisfies Z R~_1 Z^* = w R~_1, Z^d = 1 -> Stone-von Neumann for Z_d: R~_1 = X (x) 1_M, Z = Z_0 (x) 1_M.
  W commutes with Z, so W = sum_k |k><k| (x) W_k, W_k^4 = w^k.  V_k := z^{-k} W_k has V_k^4 = 1.
* Link traces: tau(R~_1^n R~_2^-n) = (1/d) z^-n [ sum_{k>=n} tr(V_{k-n}V_k^*) + i sum_{k<n} tr(V_{k-n+d} V_k^*) ]; all four
  chain terms coincide (covariance).  Result:
        S = 2(d-1) - (2/d) F(V),   I_d = 4 F(V)/(d(d-1)),
        F(V) = sum_{0<=j<k<=d-1} Re[ h_{k-j} tr(V_j V_k^*) ],  h_m = sec psi_m - i csc psi_m.
  Conversely every family of order-4 unitaries V_k is a max-ent strategy (R_i = W^{i-1}(X (x) 1)W^{-(i-1)}).
  DKZ = (V_k = 1): F_DKZ = sum_m (d-m) sec psi_m, i.e. I_ME(d).
* Checked numerically on random strategies (twirl_check.py): S, S(twirl), covariance errors ~1e-16, V^4 = 1, S from F.
* Equivalent forms (check_forms.py): E_k = V_k^*: F = sum_{j!=k} kappa(k-j) tau(E_j^* E_k), kappa(m) = h_m/2 (m>0),
  kappa(-m) = conj kappa(m); Et_k = z^-k E_k (Et_k^4 = w^-k): F = d(d-1) - 2 sum_{q=0}^{d-1} q ||Ehat_q||^2 (mean aliased
  frequency).  Y_k = z^k E_k (anti-periodic Y_{k+d} = -Y_k): F = sum_{j!=k} (-i/sin(pi(k-j)/d)) tau(Y_j^* Y_k)
  (anti-periodic discrete Hilbert kernel, eigenvalues d-1-2p, p = 0..d-1).
* Symmetries of the reduced problem: shift E_k -> E_{k+1}, E_{d-1} -> i E_0 (order 4d, s^d = global phase i); antiunitary
  E_k -> conj(E_{-k}).
* Cyclic (arc) form of the CLASSICAL problem: with increments t_r = s_{r+1}-s_r (sum t = 1 mod 4) and u_r = z i^{t_r}:
  F = sum over all d(d-1) arcs A of Im(prod_{r in A} u_r)/sin(pi |A|/d).

## 2. Dead ends (fast kills)
* Bochner relaxation (drop V^4 = 1, keep only stationarity after a second twirl): max = d(d-1) (i.e. S >= 0): useless.
  (bochner_look.py; Phi(theta_l) = d(d-1) - 2 d l exactly.)
* Level-1 (projector Gram) relaxation of the reduced problem: loose (d=3: 4.732 vs 4.309; d=9: 61.96 vs 52.86).
  (reduced_l1.py)
* Per-pair bounds: sum_{j<k} max(sec, csc) > F_DKZ; the pairs cannot be optimised separately.
* 'Simultaneous majorisation' (DKZ maximising every low-pass energy C_r = sum_{p<r}|zhat_p|^2): false for r = d-1.
* Second twirl of the reduced problem (over its own Z_{4d}) gives back the original 4-cycle problem (Fourier duality
  Z_d <-> Z_4); no new information.

## 3. Classical reduced problem (classical_z4.py, brute force d = 2..10)
Max over s in Z_4^d = F_DKZ exactly; maximisers = the d one-step configurations s = (0^a,1^(d-a)) (+ global phase);
gap to the next value: 2.31, 2.16, 2.10, 2.07, 2.05, 2.04, 2.03, 2.02 (d = 3..10).

## 4. Level-2 relaxation of the reduced problem (symmetry-reduced; reduced_sym.py / reduced_sym2.py / run_ipm.py)
Words 1, E_0^a, E_0^a E_k^b (normalised: first letter at site 0); functional twirled over Z_{4d}; 4d momentum blocks
B_p[u,v] = sum_{r<d} e^{2 pi i p r/(4d)} L(u^* s^r(v)), p = -charge mod 4; sizes 3d-2 (p=0), 3d-3 (p=0 mod 4), 2d-1 else.
Formulation checked on DKZ and random strategies (check_feas.py: PSD, objective = F).  Solved with ipm_fast.py (own
re-implementation of the lab IPM with sparse vectorised Schur complement):
   d:            3        4        5        6        7        8        9
   bound-F_DKZ: +2.8e-9  +1.9e-10 +5.2e-8  +9.2e-8  +2.2e-8  +7.7e-8  +5.0e-7    (tight to solver accuracy)
DKZ moment blocks all have rank exactly 1 (dkz_kernels.py): the optimum is a classical point of the clock model.

## 5. More numerics on the reduced relaxation (2026-09-29 evening)
* Level 2 (full, all words E_0^a E_k^b) tight to solver accuracy for d = 3..11:
  bound - F_DKZ = +2.6e-9, +2.7e-9, +5.7e-8, +9.3e-8, -1.2e-9, +6.8e-8, +5.0e-7, -1.1e-7, +2.5e-7 (d=3..11).
  (d = 12 exceeds the 1 GB cap with ipm_fast.)  Files: run_ipm.py, run_ipm_*.log, ipm_fast.py (memory-lean).
* Which blocks are needed (block_subset.py): charge-0 blocks (p = 0 mod 4) alone: tight d <= 5, NOT d = 6 (+4.5e-4);
  charges {0,2}: tight d <= 8; charges {0,1,3}: tight d <= 8.
* Correlative sparsity (sparse_test.py): PSD only on 3-site cliques {0,k,l}: tight d = 3,4,5; fails d = 6 (+0.68),
  7 (+2.64), 8 (+5.51).  2-site cliques: loose already at d = 3.  => no 3-body decomposition for d >= 6.
* Short-range words (range_test.py; words E_0^a E_k^b with cyclic distance <= r): tight only when r >= floor(d/2)
  (i.e. the full word set).  Long-range words are essential (Hilbert-type kernel).
* The SOS Gram blocks from the IPM (get_sos.py, analyze_sos.py, dkz_vec.py): rank n_p - 1 exactly, kernel = the
  DKZ vector (strict complementarity).  Closed form of the DKZ block vector (p = 0): v[E_0^a E_k^{-a}] ~ (d-k) + i^a k.

## 6. Dead end: sequential-measurement rounding (seq_round.py)
F(E) <= E_seq[F(s)] (s = outcomes of Lueders measurements of E_0,...,E_{d-1} in this order on the maximally mixed
state) is FALSE (random E: F - E_seq F up to +3.0 at d=4; also near the optimum).  Non-adjacent pairs interfere.

## 7. THEOREM L (local optimality of DKZ for ALL d) -- PROVED  [local_hessian.py = numerical confirmation]
Setting: reduced problem, point E^(0) = direct sum over one-step configurations s (s_k = c + [k >= a]) with arbitrary
multiplicities (this is the image of DKZ (x) 1_M under the reduction).  Perturb E_k(eps) = e^{i eps H_k} E_k e^{-i eps H_k}
(every nearby order-4 tuple is of this form).  Then F(eps) = F_DKZ + eps^2 F2 + O(eps^3) with
   F2 = -(1/N) sum_{config pairs s != s'} sum_{j<k} w_jk(s,s') || h_j - h_k ||_HS^2,   h_k := (H_k)_{s s'} block,
   w_jk(s,s') = 2 Re[kappa(k-j) conj(D_j) D_k],  D_k = i^{s_k} - i^{s'_k}
             = J(s_k-s_j) + J(s'_k-s'_j) - J(s'_k-s_j) - J(s_k-s'_j),  J(t) = J_{k-j}(t) = 2 sin(pi((k-j)/d + t)/2)/sin(pi(k-j)/d).
(First order vanishes since E^(0) is diagonal.)  Enumerating x = s_j, eps = s_k - s_j in {0,1}, eps' = s'_k - s'_j in {0,1},
delta = s'_j - s_j != 0, and requiring D_j, D_k != 0 gives
   w in { 2 sec, 4 sec, 2 csc, 4 csc, 2(sec + csc) }(psi_{k-j})  >=  2 sec(pi/(2d)) > 0.
Hence F2 <= 0, with F2 = 0 iff h is constant on supp(D) for every pair (the graph on supp(D) is complete with positive
weights), i.e. iff (H_k) is, modulo the commutants of the E_k, a common H: the tangent space of the unitary-conjugation
orbit.  Morse-Bott => E^(0) is a STRICT local maximum of F modulo unitary conjugation, for every d; via the (continuous)
reduction map, DKZ (x) 1_M is a strict local maximum of I_d on maximally entangled states modulo local unitaries, for
every d and M.  Minimal pair weight 2 sec(pi/(2d)) (numerically confirmed d = 3..8), matching the numerically observed
Hessian gap 2 sec(pi/(2d))/d of the lab (normalisation not re-derived).

## 8. THEOREM B-T3: the reduction for ALL pure states (weighted clock model) -- PROVED; explains the sec kernel
State psi = sum_x sqrt(lam_x)|xx> (Schmidt form WLOG), <A (x) B> = beta(A, B^T), beta(X,Y) = Tr(L^{1/2} X L^{1/2} Y),
L = diag(lam); beta is symmetric, so rho is still a symmetry.  Twirl with the direct-summed state L~ = (+)_orbit L/(4d):
W commutes with L~, hence with Z = W^4; Stone-von Neumann as before; L~ = sum_k |k><k| (x) L_k, [L_k, V_k] = 0.  Result:
     S = 2(d-1) - 2 sum_{j<k} Re[ h_{k-j} Tr(Y_j^* Y_k) ],   I_d = 4/(d-1) * sum_{j<k} Re[ h_{k-j} Tr(Y_j^* Y_k) ],
     Y_k = E_k L_k^{1/2}: normal operators whose spectrum lies on the four rays R_+ {1,i,-1,-i}, sum_k Tr(Y_k^* Y_k) = 1,
and conversely (any such Y's come from a pure state + PVMs).  Max-ent <-> L_k = 1/(dM).  Checked numerically
(twirl_allstates.py: S from the bipartite formula, after the twirl, and from the Y's agree to 1e-12; off-diagonal blocks of
L~ and [L_k, V_k] vanish to 1e-16).
* Classical aligned configuration Y_k = y_k (y >= 0, |y| = 1): value sum_{j<k} sec(psi_{k-j}) y_j y_k = (y^T K_d y - 1)/2 with
  K_d = [sec(pi(k-l)/(2d))] (diagonal 1) = the REAL PART of the clock-model couplings h on the aligned configuration.
  Uniform weights (forced by max-ent) -> (1^T K_d 1/d - 1)/2 -> I_ME(d); free weights -> (lambda_max(K_d) - 1)/2 ->
  I*(d) = 2(lambda_max(K_d)-1)/(d-1) (ADGL).  This is WHY sec kernels appear (task item 5): K_d is the clock-model coupling
  matrix restricted to the aligned (DKZ) configuration; the csc (imaginary) parts drop out there.
* Consequently the all-states Tsirelson conjecture <=> the weighted clock model's quantum value equals its classical
  value (lambda_max(K_d)-1)/2, attained at the aligned configuration with Perron weights.

## 9. Classical Sherali-Adams and the Ising split (2026-09-29 19:40)
* classical_sa.py: LP over triple-consistent laws of differences; exact for d <= 8, fails d >= 9 (+0.54, +2.52, +5.23, +8.53
  for d = 9..12).  Pairs-only LP ~ 2 F_DKZ.  => even classically no fixed-size local certificate.
* Ising split (exact identity, B-L3 on the board): F_DKZ - F = sum_{j<k} delta'(1 - tau(A_j A_k)) + 2 sum (sec pi(2) + csc pi(3)),
  A = E^2 involutions, delta' = (sec - csc)/2.  Verified by expanding in pi(t).

## 10. PROOFS (complete arguments for the theorems claimed by B)

### P1. Theorem B-T1 (exact reduction, max-ent).
Setting: R_i in U(D), R_i^d = 1, tau = Tr/D.  S(R) := 2(d-1) + sum_{n=1}^{d-1} c_n Lambda_n(R),
Lambda_n(R) := tau(R1^n R2^-n) + tau(R2^n R3^-n) + tau(R3^n R4^-n) + w^-n tau(R4^n R1^-n), c_n = -1/(1-w^-n); I_d = 4 - 2S/(d-1).
(a) Invariances.  For rho R := (R2,R3,R4,wR1) the four summands of Lambda_n are permuted cyclically (the last one:
    w^-n tau((wR1)^n R2^-n) = tau(R1^n R2^-n)), so S(rho R) = S(R); S(wR) = S(R) (links invariant); S is affine under direct
    sums with weights D_i/sum D_i.  rho^4 R = wR, hence rho^{4d} = id.
(b) Twirl.  R~ := (+)_{r in Z_{4d}} rho^r R on C^{4d} (x) C^D; S(R~) = S(R).  With the summand shift W (W e_{r+1} (x) v = e_r (x) v)
    one has W R~_i W^* = R~_{i+1} (i <= 3) and W R~_4 W^* = w R~_1 because (rho^{r+1}R)_i = (rho^r R)_{i+1}, (rho^{r+1}R)_4 =
    w (rho^r R)_1.  W^{4d} = 1.  Z := W^4 satisfies Z R~_1 Z^* = w R~_1 and Z^d = 1.
(c) Stone-von Neumann (elementary): H_k := ker(Z - w^k); R~_1 H_k = H_{k+1}; choose an ONB f_mu of H_0 and put
    e_{k,mu} = R~_1^k f_mu (k = 0..d-1): an ONB of C^{4dD} (so 4dD = dM) in which R~_1 = X (x) 1_M, Z = Z_0 (x) 1_M.
    W commutes with Z = W^4, so W = sum_k |k><k| (x) W_k, W_k^4 = w^k; put V_k := z^-k W_k (k = 0..d-1), V_k^4 = 1.
(d) Evaluation.  R~_{i} = W^{i-1} R~_1 W^{-(i-1)} and W^4 R~_1^{-n} W^-4 = w^-n R~_1^-n give that all four summands of Lambda_n
    equal tau(L_n), L_n := R~_1^n W R~_1^-n W^* = sum_k |k><k| (x) W_{k-n} W_k^*  (indices mod d).  With W_k = z^k V_k and
    z^d = i:  tau(L_n) = (1/d) z^-n [ sum_{k>=n} tr(V_{k-n} V_k^*) + i sum_{k<n} tr(V_{k-n+d} V_k^*) ].
    Using 4 c_n z^-n = i e^{i psi_n}/(sin psi_n cos psi_n) (psi_n = pi n/2d) and, for the second sum, n = d - (j-k),
    psi_{d-m} = pi/2 - psi_m, tr(V_j V_k^*) = conj tr(V_k V_j^*), both sums combine into
       4 sum_n c_n tau(L_n) = -(2/d) sum_{j<k} Re[ (sec psi_{k-j} - i csc psi_{k-j}) tr(V_j V_k^*) ],
    i.e. S = 2(d-1) - (2/d) F(V) and I_d = 4F/(d(d-1)).
(e) Converse: for any order-4 V_k on C^M put W := sum_k |k><k| (x) z^k V_k, R_i := W^{i-1}(X (x) 1)W^{-(i-1)} (i = 1..4); these are
    order-d unitaries (conjugates of X (x) 1), W R_4 W^* = W^4 (X(x)1) W^-4 = w (X (x) 1) because W^4 = Z_0 (x) 1, and step (d)
    applies verbatim (it used only this covariant form).  Hence sup_R I_d = 4 sup_V F/(d(d-1)).  []

### P2. Theorem B-T2 (local optimality, all d).
Let E^0_k = sum_s i^{s_k} Pi_s (Pi_s = projection onto the multiplicity space of the one-step configuration s; any finite
multiset of one-step configurations s_k = c + [k >= a]).  (i) Nearby tuples: if E_k^4 = 1 and ||E_k - E^0_k|| < sqrt2 then the
spectral projections of E_k are close to those of E^0_k and have the same ranks, so E_k = U_k E^0_k U_k^* with U_k near 1
(e.g. U_k = e^{iH_k}, H_k Hermitian, unique modulo the commutant of E^0_k).  (ii) Expansion: E_k(eps) = E_k + eps A_k + eps^2 B_k +
O(eps^3), A_k = i[H_k,E_k], B_k = -(1/2)[H_k,[H_k,E_k]].  tau(A_j^* E_k) = 0 = tau(E_j^* A_k) (A has zero diagonal blocks, E is
block-diagonal), so the first variation of every tau(E_j^*E_k) vanishes.  Second order: with h_k := (H_k)_{s s'} and
D_k := i^{s_k} - i^{s'_k}, the (s,s')+(s',s) contributions give
     d^2/(2 d eps^2) tau(E_j^* E_k) = (1/N) sum_{s<s'} conj(D_j) D_k ( <h_j,h_k> + <h_k,h_j> - ||h_j||^2 - ||h_k||^2 )
                                  = -(1/N) sum_{s<s'} conj(D_j) D_k ||h_j - h_k||_HS^2,
so F2 = -(1/N) sum_{s<s'} sum_{j<k} w_jk(s,s') ||h_j - h_k||^2 with w_jk = 2 Re[kappa(k-j) conj(D_j) D_k]
= J(s_k-s_j) + J(s'_k-s'_j) - J(s'_k-s_j) - J(s_k-s'_j).  (iii) Weights: for one-step s, s' and j<k, eps := s_k - s_j,
eps' := s'_k - s'_j in {0,1}; delta := s'_j - s_j.  D_j != 0 iff delta != 0, D_k != 0 iff delta + eps' - eps != 0.  With
J(0) = sec, J(1) = csc, J(2) = -sec, J(3) = -csc (at psi = psi_{k-j}):
   (0,0): w = 2sec - J(delta) - J(-delta) = 2sec, 4sec, 2sec (delta = 1,2,3);  (1,1): 2csc - J(delta+1) - J(1-delta) = 2csc, 4csc, 2csc;
   (0,1): delta in {1,2}: w = sec + csc - J(delta+1) - J(-delta) = 2(sec+csc);  (1,0): delta in {2,3}: 2(sec+csc).
All weights are >= 2 min(sec psi, csc psi) >= 2 sec(pi/2d) > 0.  (iv) Hence F2 <= 0; F2 = 0 iff for each pair (s,s') the vector
h is constant on supp D (the weighted graph on supp D is complete with positive weights), iff (H_k) agrees with a single H
modulo the commutants, i.e. iff the perturbation is tangent to the orbit {U E^0 U^*}.  Morse-Bott (critical manifold = the
unitary orbit, Hessian negative definite on the normal space) => strict local maximum modulo unitary equivalence.  Since the
reduction map R -> V of P1 is continuous (the basis e_{k,mu} = R~_1^k f_mu depends continuously on R, with Z and f_mu fixed),
DKZ (x) 1_M is a strict local maximum of I_d over max-ent strategies of the same dimension, modulo local unitaries.  []

### P3. Theorem B-T3 (all pure states).
Identical to P1 with tau replaced by the symmetric bilinear form beta(X,Y) = Tr(L^{1/2} X L^{1/2} Y) (Schmidt form of psi;
<A (x) B> = beta(A, B^T)); the twirled state L~ = (+)_r L/(4d) commutes with W, hence with Z; so L~ = sum_k |k><k| (x) L_k and,
because L~ commutes with W, [L_k, W_k] = 0.  In step (d) beta~(R~_1^n, R~_2^-n) = sum_k Tr(L_k^{1/2} L_{k-n}^{1/2} W_{k-n} W_k^*),
and Tr(L_k^{1/2} L_j^{1/2} V_j V_k^*) = Tr(Y_j^* Y_k) with Y_k := V_k^* L_k^{1/2} = E_k L_k^{1/2} (L_k commutes with V_k).  The same
trigonometric algebra gives S = 2(d-1) - 2 sum_{j<k} Re[h_{k-j} Tr(Y_j^* Y_k)] with sum_k Tr(L_k) = 1.  Converse as in P1(e)
with the state (+)_k L_k.  Aligned classical configuration Y_k = y_k 1 (y_k >= 0): sum_{j<k} sec(psi_{k-j}) y_j y_k =
(y^T K_d y - 1)/2 (diag K_d = 1), maximised by the Perron vector: I = 2(lambda_max(K_d)-1)/(d-1) (ADGL).  []

## 11. Landscape (dead end for a 'benign landscape' proof)  [landscape.py, self_laplacian.py]
* WLOG balanced spectra (global-phase twirl).  Riemannian ascent over E_k = U_k Z_4 U_k^* (M = 4), random starts:
  d = 4,5,6: best = F_DKZ; d = 7: 13/16 runs reach F_DKZ (commuting), 3/16 stop at F = 28.8333 with ||[E_j,E_k]|| = 1.41:
  a NON-classical local maximum below F_DKZ.  Random (unbalanced) spectra: many poor local maxima, some non-commuting.
* Classical mixtures: at a balanced sum of the phase copies of one configuration s, the Hessian blocks are weighted Laplacians
  with weights |1-i^c|^2 J_{k-j}(s_k - s_j); configurations with PSD 'self-Laplacian' = one-step ones for d = 3,4,6 but also
  non-optimal ones for d = 5 (F = 10, borderline), 7 (e.g. s = (0,3,0,0,0,0,0), F = 28.68, lambda_2 = 0.70), 8 (8 extra).
  => spurious local maxima exist (both classical mixtures and non-commuting ones); no 'every local max is classical/optimal'
  argument can work.  The global statement genuinely needs a certificate (or a global inequality).

## 12. SATWAP in the reduced picture (negative check; satwap_reduced.py)
SATWAP (m = 2, weights alpha_k, beta_k of Salavrakos et al. as reconstructed) is also a rho-symmetric 4-cycle expression, so it
reduces to a clock model with coupling H_jk = (4/d) what_{k-j} z^{-(k-j)}.  Its top eigenvector is again the non-Z_4 plane wave
z^k (phases pi k/(2d)), and d*lambda_max (1.82, 2.01, 2.13, ... for d = 3,4,5) exceeds the value at V == 1 (4/3, 3/2, 8/5 =
2(d-1)/d).  So the reduced Bochner relaxation is loose for SATWAP too: SATWAP's per-n SOS lives in the ORIGINAL power
structure, it does not come from a level-1 argument in the reduced picture.  (No claim about SATWAP made on the board.)

## 13. Classical problem as a pinned particle system (Aubry-Mather analogy; partial)
Lift x_k := k + d s_k (s_k in Z, s_{k+d} = s_k + 1, so x_{k+d} = x_k + 2d, x_k = k mod d).  Then
     F(s) = sum_{j in Z_d} sum_{m=1}^{d-1} sin(pi (x_{j+m} - x_j)/(2d)) / sin(pi m/d)
          = (d/pi) sum_{j in Z_d} sum_{n != 0} sin(pi (x_{j+n} - x_j)/(2d)) / n        (lattice Hilbert-transform form).
Ideal (Bochner) x_k = 2k is forbidden by the pinning x_k = k mod d.  ORDERED (monotone) lifted configurations have increments
x_{k+1}-x_k in 1 + dZ, positive, summing to 2d per period => exactly one increment 1+d: ORDERED <=> ONE-STEP (DKZ orbit).
Pair terms are supermodular (sin concave) exactly when all differences lie in [0,2d], i.e. again only on the one-step set;
the interaction is periodic in the difference (Z_4), so Aubry's crossing lemma / twist-map theory do not apply globally.
Dead end as a proof, but it explains the optimal set: Beatty/Sturmian rounding x_k = k + d floor((k+a)/d) of the slope-2 line.
* Diagonal SOS ansatz G_p = Pi^perp Diag(x) Pi^perp: infeasible (d = 3,4,5; diag_ansatz.py).
* Fourier look at the analytic-centre SOS (fourier_look.py): eigenvectors only partially concentrated in the distance
  frequency; no visible structure.

## 14. Classical telescoping lemma (PROVED, all d) -- covers all near-optimal classical configurations  [lift_lemma.py]
Arc/csc form (exact, any classical s in Z_4^d, cyclic with s_{k+d} = s_k + 1):
    F(s) = sum_{j in Z_d} sum_{m=1}^{d-1} csc(pi m/2d) chi(n_{j+m} - n_j),   chi(t) = [t = 1 mod 4] - [t = 3 mod 4],
for ANY integer lift n of s with n_{k+d} = n_k + 1.  Telescoping: sum_{j in Z_d} (n_{j+m} - n_j) = m for every lift and m.
Hence  F(s) - F_DKZ = sum_m csc(pi m/2d) sum_j [chi(D_m(j)) - D_m(j)],  D_m(j) := n_{j+m} - n_j,  and chi(D) - D <= 0 for D >= -1.
LEMMA: if some lift satisfies n_k - n_j in [-1, 2] for all 0 <= j < k <= d-1 (equivalently all D_m(j) >= -1), then
F(s) <= F_DKZ, with equality iff all D_m(j) in {0,1}, i.e. iff s is a one-step (DKZ) configuration.
Numerics (greedy lift): the liftable class has exactly d 2^{d-1} elements (d = 3..9) and contains every configuration with
F > F_DKZ - 4.5 (d=4) ... (all configs within the top ~15% of the range); the best NON-liftable configuration has
F = -0.31, 3.93, 10.0, 15.49, 23.69, 32.43, 43.69 for d = 3..9 (vs F_DKZ = 4.31, 8.69, 14.55, 21.90, 30.74, 41.05, 52.86).
Per-m versions are false (all-(+1) configuration, m = 1), so any complete proof must couple different m through the csc
weights.  The telescoping identity is the classical shadow of A's grade-sharp additivity (A-L1).

## 15. Quantum telescoping attempt (dead end)  [trace_ineq.py]
Operator form of the telescoping: for Hermitian lifts N_k (integer spectra, e^{i pi N_k/2} = E_k, N_{k+d} = N_k + 1) one has the
exact OPERATOR identity sum_{j in Z_d} (N_{j+m} - N_j) = m 1, hence
    F_DKZ - F = sum_m csc(pi m/2d) sum_j tau[ (N_{j+m} - N_j) - Im(E_j^* E_{j+m}) ]   (any lifts).
The termwise trace inequality  Im tr(e^{-i pi A/2} e^{i pi B/2}) <= tr(B - A)  for integer-spectrum A, B with B - A >= -1 is FALSE
(random search: violation 0.18 with spectra {-1,0,1}).  With the twirl the per-m sums are pi_m(1)-pi_m(3) <= m/d, false in
general; any proof must couple different m through the csc weights.

## 16. SUMMARY (B, 2026-09-29 ~21:30)
PROVED, all d:
  B-T1  exact reduction: max-ent CGLMP_d  ==  tracial Z_4 clock model  F(V) = sum_{j<k} Re[h_{k-j} tr(V_j V_k^*)],
        h_m = sec(pi m/2d) - i csc(pi m/2d), V_k^4 = 1;  I_d = 4F/(d(d-1));  DKZ <-> V == 1 (classical ground state).
  B-T3  same for all pure states: weighted clock model (Y_k = E_k L_k^{1/2}); K_d = coupling on the aligned configuration;
        ADGL = aligned classical optimum with Perron weights.
  B-T2  DKZ (x) 1_M is a STRICT local maximiser of I_d on max-ent states (modulo local unitaries), every d, M; explicit Hessian
        = sum over pairs of one-step configurations of complete-graph Laplacians with weights in {2sec,4sec,2csc,4csc,2(sec+csc)}.
  B-L3  Ising split identity;  B-L4 classical telescoping lemma (all near-optimal classical configurations);  B-R1 lifted
        (Aubry-Mather) form: ordered <=> one-step.
NUMERICAL: reduced level-2 relaxation tight d = 3..11; no k-local certificate (quantum 3-cliques fail d >= 6, classical SA3
  fails d >= 9); spurious local maxima (d = 7); SOS Gram blocks dense.
OPEN: the global all-d inequality F(V) <= F_DKZ for non-commuting order-4 V (equivalently OPT_d); even the classical
  (commuting) case is open for d >= 11 (brute force d <= 10; SA3 LP d <= 8).
