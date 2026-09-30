# OQP 27B — The power of CGLMP inequalities: working log

Folder: `iqoqi/programs/oqp27B/`. Python: `python` (numpy, scipy, clarabel, scs, mpmath, sympy).

## 0. Problem (IQOQI OQP 27B)
Show that the Durt–Kaszlikowski–Zukowski (DKZ) / CGLMP Fourier-type observables are necessarily
optimal for the CGLMP inequality on the maximally entangled state; show they give the highest
noise resistance and the best Kullback–Leibler discrimination. Known: numerical (NPA) optimality;
analytic SOS certificates of the all-state Tsirelson bound for d = 3, 4 (Ioannou–Rosset,
arXiv:2112.10803).

## 1. Novelty check (2026-09-29)
- arXiv full-text search "CGLMP" (newest first, 43 hits back to 2015): no 2022–2026 paper gives an
  analytic/exact CGLMP maximum for d >= 5 or for all d, nor a proof of DKZ optimality on the
  maximally entangled state. Only 2112.10803 (d = 3, 4, all states).
- arXiv:2606.21626 (Coccia–Padovan–Vallone, Jun 2026, SOS Tsirelson bounds in arbitrary dimension):
  treats SATWAP / MUB-type expressions, no CGLMP.
- arXiv:2601.02893 (Hsu et al. 2026): symmetry vs dimension trade-off; no exact CGLMP maxima.
- arXiv:2604.03700: CHSH mod 3 exact value; no CGLMP.
- Zohren–Reska–Gill–Westra EPL 90, 10002 (2010), arXiv:1003.0616: d -> infinity only (bound 0 in
  the probability form); not a finite-d result.
- Cha & Lee (SNU) 2026 OQP work (per lit/S7): OQP 37, 40, trine, spin-1 CHSH — nothing on OQP 27.
- Semantic Scholar citations of 2112.10803 (re-queried later today, 11 papers): 2609.05555 (I3322),
  2608.16512 (MoMPy, no CGLMP), 2606.19033, 2601.02893, 2504.00162, 2310.00612, 2307.15688, 2303.02127,
  2210.07263, 2203.05372, OJMO entry -- none gives exact CGLMP values for d >= 5 or max-ent optimality.
Verdict: open as stated (all-d, and d >= 5 exact).

## 2. Notation / reformulation
Probability form (all no-signalling boxes):
  Pi := P(X1<Y1) + P(Y1<X2) + P(X2<Y2) + P(Y2<=X1) >= 1 (local),
  Gill form S = E m(X1-Y1)+E m(Y1-X2)+E m(X2-Y2)+E m(Y2-X1-1) = d*Pi - 1,
  CGLMP I_d = 4 - 2(d*Pi - 1)/(d-1)   (local bound I_d <= 2).
(verified numerically: conventions2.py, see section 3)

## 3. Closed forms, conventions (verified: basics.py, conventions2.py, dkz_scan.py, exact_values.py)
- Standard CGLMP (literal CG-form, Collins et al.) satisfies I_d = 4 - 2(d*Pi_std - 1)/(d-1) with
  Pi_std = P(X1<Y1)+P(Y1<=X2)+P(X2<Y2)+P(Y2<X1). We use the equivalent chain form (settings swapped)
  Pi = P(X1<Y1)+P(Y1<X2)+P(X2<Y2)+P(Y2<=X1) (local min 1), and the Gill form S = d*Pi - 1.
- DKZ measurements in this convention: Alice |a>_x = d^-1/2 sum_k w^{k(a+alpha_x)}|k>, alpha=(1/2,0);
  Bob |b>_y = d^-1/2 sum_k w^{-k(b+beta_y)}|k>, beta=(1/4,-1/4).
- K_d := [sec(pi(k-l)/(2d))]_{k,l=0..d-1}. Then (Zohren et al. formula, re-derived numerically):
    DKZ on max. entangled:  Pi_ME   = 2 - 1^T K_d 1/d^2 ;   I_ME   = 2(1^T K_d 1/d - 1)/(d-1)
    Acin et al. optimum:    Pi_ADGL = 2 - lambda_max(K_d)/d ;  I_ADGL = 2(lambda_max(K_d) - 1)/(d-1)
- Minimal polynomials (exact_values_out.txt):
    d=3: lambda_max(K_3) root of 3x^2-12x+1 (I = 1+sqrt(11/3)); s_ME=1^TK1/9 root of 243x^2-378x+83
    d=4: lambda_max degree 8: x^8-8x^7-12x^6+184x^5-90x^4-760x^3-140x^2+72x+1 ; s_ME: 256x^4-256x^3-288x^2+208x+31
    d=5: lambda_max degree 6: 5x^6-30x^5-185x^4+940x^3+819x^2-526x+1 ; s_ME: 390625x^4-62500x^3-596250x^2-80100x+24641
    d=6: lambda_max degree 12 (see file) ; s_ME: 944784x^4-1469664x^3-538488x^2+987768x+94681
  Values: I_ADGL = 2.9148542155, 2.9726982671, 3.0157104755, 3.0497004192 (d=3..6);
          I_ME   = 2.8729340512, 2.8962432185, 2.9105448081, 2.9202036064 (d=3..6).

## 4. Numerical relaxations (own Fourier/charge-reduced moment code + own IPM, ipm.py)
- Bipartite NPA level 1+AB (npa_bip.py): equals ADGL to <=5e-9 for d = 2..7 (run_ipm_test.py bip).
- Maximally entangled of ANY dimension D <-> tracial problem (tau(A B^T)); R1=A1,R2=B1^T,R3=A2,R4=B2^T.
  * tracial level 1: very loose (Pi even negative at d=5).
  * tracial level 2 (all words R_i^n R_j^m), with Z_d, reversal, chain rotation rho, reflection sigma
    symmetry (tracial_sym2.py): equals DKZ value to <=3e-9 for d = 3..8 (ipm_full_78.txt).
  * charge-0 word subset {1, R_i^q, R_i^n R_j^-n}: tight for d<=5, NOT for d=6 (gap 2.5e-5).
- DKZ max-ent optimum is the discrete geodesic R_i = G^{3-i} X^{-1} G^{i-3}, G = diag(zeta_{4d}^k):
  each step conjugates by the fractional clock Z^{-1/4}; R_5 = w R_1.
- Strict complementarity observed: per block rank(X)+rank(Z)=n; DKZ primal rank d per charge sector.

## 5. EXACT CERTIFICATES (2026-09-29)
Pipeline (both problems): exact symmetry-reduced moment relaxation over Q(zeta_{4d}) (cyclo.py);
exact optimal point (DKZ as generalised permutation matrices; Perron vector over L = Q(zeta_{4d})(lambda*));
exact kernels of the symmetrised optimal moment blocks -> face X_b = N_b Y_b N_b^+; numerical reduced SDP
(ipm.py); rounding of Y to Gaussian rationals (2^-40); exact projection onto the SOS identities; exact LDL.
Independent checkers (Fraction arithmetic, own word/symmetry code): verify_tracial.py, verify_bip.py.
- Maximally entangled, ANY local dimension (tracial level 2): cert_tracial_d{3,4,5}.pkl VERIFIED
  (verify_tracial.py): I_d <= I_ME(d) for all max-ent states and PVMs; DKZ attains it.
- All states (NPA 1+AB): cert_bip_d{3,4,5}.pkl VERIFIED (verify_bip.py): I_d <= 2(lambda_max(K_d)-1)/(d-1).
  d=5 is new (Ioannou-Rosset stopped at d=4): I_5* = (lambda-1)/2, lambda root of
  5x^6-30x^5-185x^4+940x^3+819x^2-526x+1 (largest), I_5* = 3.01571047552267328552398380196.
- 2026-09-29: verify_tracial.py 5 6 -> VERIFIED (d=5: 122 classes/33889 terms, min pivot 9.66e-4;
  d=6: 220 classes/67001 terms, min pivot 1.81e-4). Max-ent optimality of DKZ now proven for d=3..6.
- Exact linear-algebra swell for the Perron field (d=6 bipartite, degree-24 field) solved with a
  multi-modular solver (modsolve.py: primes splitting completely in L, vectorised elimination over all
  embeddings, CRT + rational reconstruction, exact re-check). d=4,5 re-certified with it in 2s/31s.
- Word-set reduction: the tracial level-2 relaxation restricted to ADJACENT pairs R_i^n R_{i+-1}^m
  (plus 1, R_i^n) is still tight numerically for d=5..8 (subset_test.py); used for d=7,8 data.
- bipartite d=6: cert_bip_d6.pkl (pivot-subsystem projection + modsolve, 63 s) VERIFIED by verify_bip.py:
  I_6* = 2(lambda_max(K_6)-1)/5 = 3.04970041924084822234629611391, lambda_max(K_6) of degree 12.
  (Tsirelson bound for d=6 now proven.)

## 6. Rigidity (self-testing) for the maximally entangled case -- argument
Support of the tracial SOS in sector 0 (all Y_b > 0 on the kernels) gives, for ANY optimal max-ent
strategy (faithful trace => f(R)=0 for every f in the SOS support):
 (E_n) kappa != 1 blocks have full kernels => all chain modes of the link unitaries vanish =>
       R1^n R2^-n = R2^n R3^-n = R3^n R4^-n = w^-n R4^n R1^-n =: U_n   (equal links), n=1..d-1;
 (T_n) kappa = 1 block kernel = span of the two-eigenvalue relations => (U_n - z^-n)(U_n - i z^-n) = 0,
       i.e. z^n U_n = 1 + (i-1) Pi^(n), Pi^(n) a projection  (z = zeta_{4d}).
Derivation: R2 = z D R1 with D = 1-(1+i)Pi, Pi := Pi^(1); z^n U_n = prod_{j=n-1..0} (1+(i-1) R1^j Pi R1^-j).
Two-valuedness of these products for all n (Jordan lemma for two projections: a product of two
"quarter-phase" unitaries 1+(i-1)P, 1+(i-1)Q has spectrum in {1,i} iff PQ = 0) gives inductively that
P_j := R1^j Pi R1^-j (j=0..d-1) are mutually orthogonal, and U_d = 1 gives sum_j P_j = 1. Hence
(R1, Pi) is a system of imprimitivity: R1 = X (x) 1_M (cyclic shift), Pi = |0><0| (x) 1_M, and
R_{i+1} = z D R_i: every optimal max-ent strategy is unitarily equivalent to (DKZ) (x) 1_M, i.e. all
four measurements are Fourier bases dressed by diagonal unitaries -- the DKZ form is NECESSARY.
(To be machine-checked per d: kernel structure of the sector-0 blocks, verify_rigidity.py.)
- tracial d=7 (adjacent word set): first certificate built with the old exact solver (1363 s); later replaced (see below).
- KL strength (vDGG, uniform settings), kl_strength.py: d=2 CHSH 0.046274 bits (matches vDGG 0.0462);
  d=3 DKZ+max-ent 0.057783 bits; DKZ + optimised symmetric Schmidt vector 0.076850 bits
  (lambda ~ (0.6475, 0.4019, 0.6475)). Global random-restart search over all PVMs+states running
  (numerical only).
- Exact values d=7,8 (exact_values_78.txt): lambda_max(K_7) degree 12, lambda_max(K_8) degree 32;
  s_ME degree 6 (d=7), 8 (d=8). I_ADGL(7)=3.0776483114342576, I_ADGL(8)=3.1012805879054686;
  I_ME(7)=2.9271609424245431, I_ME(8)=2.9324096087044598.
- Closed form (max-ent): I_ME(d) = 4/(d(d-1)) * sum_{j=1}^{d-1} (d-j) sec(pi j/(2d));  v_crit = 2/I_ME(d).
- tracial d=8 (adjacent words, modular projection): cert_tracial_d8.pkl built (550 s); verification running.
- verify_rigidity.py 3..8: support conditions (a),(b) hold for all d=3..8.
- Regenerated certificates re-verified: verify_tracial 3 4 OK; verify_bip 3 4 5 6 OK (final files).
- Tightness pattern (numerical, IPM): adjacent-pair level-2 tracial relaxation tight for d<=9
  (d=9: -3.6e-9) but NOT for d=10 (-1.34e-5, gap 1e-10). Primal search (maxent_primal.py, Riemannian
  gradient descent, random restarts): d=10 D=10 (12 restarts), d=10 D=20 (6), d=9 D=9 (8) all converge
  to the DKZ value (|diff| < 2e-13) -> DKZ still optimal numerically at d=10; the reduced word set is
  simply too small there. Full level 2 at d=10: running.
- KL local test d=3 (kl_local_d3.log): 8 perturbed starts (sigma 0.05, 0.12) with local optimisation over
  ALL projective measurements and pure states never exceed 0.07685015 bits (DKZ + optimised state);
  far starts fall into the embedded-CHSH optimum 0.04627 bits. Numerical evidence only.
- Full level-2 tracial relaxation at d=10 (nvar 1120): tight to -1.6e-9 (IPM gap 9e-10), whereas the
  adjacent-pair subset fails at d=10. Supports "tracial level 2 is tight for all d" (numerically d<=10).
- Hessian at DKZ (D=d, hessian_dkz.py): stationary (|grad| ~ 1e-12); Hessian PSD with exactly
  d^2+4d-1 zero modes (= symmetry orbit: global unitary + per-vector phases) for d=3..10, i.e. DKZ is a
  NON-DEGENERATE local optimum for d<=10 (numerical). Spectral gap = 2 sec(pi/(2d))/d to 5 digits
  (d=3..10); d=3 spectrum (x d): 4/sqrt3, 2sqrt3, 4+2/sqrt3, 2+2sqrt3, 4sqrt3, 2+10/sqrt3, 6+2sqrt3, 8+4/sqrt3.
- verify_tracial.py 8 -> VERIFIED (344 classes, 88369 expanded terms, min certified pivot 2.13e-4):
  max-ent optimality of DKZ proven for d=8 (I_ME(8) = 2.93240960870445984824840081342).
- bipartite d=7: cert_bip_d7.pkl built (1223 s, min LDL pivot 1.1e-8 > 0 exact); verification running.
- d=7 tracial re-certified with the modular pivot method (smaller numbers) for faster checking;
  (the old-method file was discarded; its verification had not finished).
- verify_tracial.py 7 (modular-method certificate, 108 s to build) -> VERIFIED (229 classes, 55777 terms,
  min pivot 3.48e-4); verify_rigidity.py 7 OK.  => max-ent optimality + rigidity proven for d = 3..8.
- Cross-n coupling is essential: splitting the kappa!=1 sector-0 blocks by the power n (or by {n,d-n})
  destroys tightness (d=4: -1.2e-2; d=5: -2.2e-3), split_test.py.  No per-n certificate structure.
## 8. Results table (values from exact algebraic numbers)

CGLMP value I_d (local bound 2).  K_d = [sec(pi(k-l)/(2d))]_{k,l=0..d-1}.
  max. entangled (any local dimension D, PVMs):  I_ME(d) = 4/(d(d-1)) * sum_{j=1}^{d-1} (d-j) sec(pi j/(2d))
  all states (Tsirelson bound):                    I*(d)   = 2(lambda_max(K_d) - 1)/(d-1)
  limit d -> infinity:  I_ME -> 32 G/pi^2 = 2.96981498168617730363 (G = Catalan's constant);  I* -> 4.

 d | I_ME(d)                          status (max-ent)            | I*(d)                            status (all states)
 3 | 2.872934051172335372024396748   EXACT CERT, verified        | 2.914854215512676219950203823   EXACT CERT, verified (known: Ioannou-Rosset)
 4 | 2.896243218458708353238334342   EXACT CERT, verified        | 2.972698267102243830469832428   EXACT CERT, verified (known: Ioannou-Rosset)
 5 | 2.910544808072077472513641868   EXACT CERT, verified        | 3.015710475522673285523983802   EXACT CERT, verified (new)
 6 | 2.920203606379063954194027462   EXACT CERT, verified        | 3.049700419240848222346296114   EXACT CERT, verified (new)
 7 | 2.927160942424543083960828253   EXACT CERT, verified        | 3.077648311434257628007210561   EXACT CERT, verified (new)
 8 | 2.932409608704459848248400813   EXACT CERT, verified        | 3.101280587905468566855035285   EXACT CERT, verified (new)
 9 | 2.936509531681428021417975460   EXACT CERT, verified        | 3.121684417680055775394584842   numerical (NPA)
10 | 2.939800297927745040587961209   numerical (level 2)         | 3.139587407734804908553579701   --

Minimal polynomials (exact_values_out.txt, exact_values_78.txt):
  s_ME = 1^T K_d 1 / d^2  (I_ME = 2(d s_ME - 1)/(d-1)):
    d=3: 243x^2 - 378x + 83
    d=4: 256x^4 - 256x^3 - 288x^2 + 208x + 31
    d=5: 390625x^4 - 62500x^3 - 596250x^2 - 80100x + 24641
    d=6: 944784x^4 - 1469664x^3 - 538488x^2 + 987768x + 94681
    d=7: 96889010407x^6 - 177959406870x^5 - 75582305911x^4 + 261073013516x^3 - 91441660023x^2 - 15030609174x + 2192217831
    d=8: 134217728x^8 - 134217728x^7 - 369098752x^6 + 333447168x^5 + 276267008x^4 - 213024768x^3 - 43293696x^2 + 16512896x + 1324033
    d=9: 847288609443x^6 - 1066956026706x^5 - 1374180474483x^4 + 1401645876804x^3 + 614410966413x^2 - 327872234322x + 20985409187
  lambda = lambda_max(K_d)  (I* = 2(lambda - 1)/(d-1)), largest real root of:
    d=3: 3x^2 - 12x + 1                      (I* = 1 + sqrt(11/3))
    d=4: x^8 - 8x^7 - 12x^6 + 184x^5 - 90x^4 - 760x^3 - 140x^2 + 72x + 1
    d=5: 5x^6 - 30x^5 - 185x^4 + 940x^3 + 819x^2 - 526x + 1
    d=6: 81x^12 - 972x^11 - 5562x^10 + 89532x^9 + 65583x^8 - 2568024x^7 + 1164756x^6 + 21967128x^5
         - 3527009x^4 - 6876380x^3 - 528634x^2 + 3884x + 1
    d=7: 49x^12 - 784x^11 - 4214x^10 + 110544x^9 - 98833x^8 - 4222624x^7 + 10647308x^6 + 25631648x^5
         - 11724545x^4 - 20568912x^3 + 2359882x^2 - 32368x + 1
    d=8: degree 32 (exact_values_78.txt)
    d=9: degree 15 (exact_values_9.txt)

White-noise critical visibility for the CGLMP witness: v = 2/I (max-ent: 0.696152 (d=3), 0.690550, 0.687157,
0.684884, 0.683256, 0.682033 (d=8)); all states: 0.686141 (d=3), 0.672789, 0.663194, 0.655802, 0.649847, 0.644895.
- hessian_fast.py (Hessian = finite differences of the analytic Riemannian gradient) for
  d = 3,4,5,8,10,12,16,20,25,30: DKZ is stationary, no negative directions, the zero modes are exactly the
  d^2+4d-1 symmetry directions, and the smallest nonzero Hessian eigenvalue equals 2 sec(pi/(2d))/d to
  6 digits in every case -> DKZ is a NON-DEGENERATE LOCAL optimum on |Phi_d> for all tested d<=30
  (numerical; conjectured spectral gap 2 sec(pi/(2d))/d for all d).
- verify_bip.py 7 -> VERIFIED (91 classes, 16177 terms, min certified pivot 1.12e-8):
  Tsirelson bound for d=7 proven: I_7* = 2(lambda_max(K_7)-1)/6 = 3.07764831143425762800721056066.
- tracial d=9 (adjacent words): cert_tracial_d9.pkl built (1091 s, min LDL pivot 5.87e-5);
  verify_rigidity 9 OK; sanity check OK; exact verification running.
- End-to-end numerical sanity checks (sanity_tracial.py d=3..9, sanity_bip.py d=3..7): on random
  strategies (symmetrised over the group images) S - lam equals the SOS expression to <=1e-13 and is >0.
- Full level-2 tracial relaxation d=9: tight (+3e-10).

## 9. FINAL STATEMENT OF RESULTS (rigour levels)
Notation: CGLMP value I_d (local bound 2); K_d = [sec(pi(k-l)/(2d))]_{k,l=0..d-1}.

(A) Maximally entangled states -- COMPUTER-ASSISTED EXACT PROOF, d = 3,4,5,6,7,8,9.
    For every local dimension D, every |Phi_D> and all projective d-outcome measurements:
        I_d <= I_ME(d) = 4/(d(d-1)) sum_{j=1}^{d-1} (d-j) sec(pi j/(2d)),
    attained by the DKZ/CGLMP Fourier measurements on |Phi_d>. Certificates: tracial (Tr(A B^T)/D) moment
    relaxation of level 2 (full words d<=6, adjacent-pair words d=7,8,9), symmetry-reduced
    (Z_d charge, chain rotation rho with twist, reflection sigma, transpose), exact SOS over Q(zeta_{4d})
    (Gram matrices on the exact kernel of the optimal moment blocks), checked by verify_tracial.py
    (independent Fraction arithmetic; identity modulo symmetry+cyclic relations; exact LDL with
    interval-certified positive pivots; exact attainment by DKZ).  Previously: numerical only (even d=3).
(B) Rigidity ("necessarily of DKZ form") -- PROOF for the same d: any maximally entangled strategy attaining
    I_ME(d) is, up to local unitaries (W (x) conj W, which fix |Phi_D>), DKZ (x) identity on an auxiliary
    maximally entangled state; in particular D is a multiple of d and all four measurement bases are
    Fourier (DFT) bases dressed by diagonal unitaries.  Proof = SOS support conditions (verify_rigidity.py)
    + the analytic argument of section 6 (equal links, two-eigenvalue relations, Jordan lemma,
    system of imprimitivity).
(C) All states (Tsirelson bound) -- COMPUTER-ASSISTED EXACT PROOF, d = 3,4,5,6,7,8:
        I_d <= I*(d) = 2(lambda_max(K_d) - 1)/(d-1)   (the Acin-Durt-Gisin-Latorre value),
    attained by DKZ measurements and psi = sum_k v_k|kk>, v = Perron vector of K_d. NPA level 1+AB,
    symmetry-reduced, exact SOS over L = Q(zeta_{4d})(lambda*), checked by verify_bip.py (independent
    arithmetic, interval-certified root and pivots, exact attainment via the restricted Bell operator
    (2d-1)1 - K_d). d=3,4 reproduce Ioannou-Rosset; d=5,6,7 are new.
(D) Noise resistance (CGLMP witness, white noise): critical visibility 2/I_ME(d) on maximally entangled
    states and 2/I*(d) over all states are optimal (corollary of (A),(C)) for the same d.
    d = 3, two maximally entangled qutrits (DKZ's setting, D = 3), ALL Bell inequalities: the facets of
    the (2,2,3) local polytope are (up to relabelling) positivity, lifted CHSH and CGLMP_3 (Collins-Gisin
    2004 classification, taken from the literature). CGLMP_3 facets: threshold >= 2/I_ME(3) = 0.696152
    by (A). Lifted CHSH: for +-1 observables on C^3 and |Phi_3>, CHSH <= (||B1+B2||_1+||B1-B2||_1)/3
    <= (4 sqrt2 + 2)/3 (Hoelder + Jordan lemma: one 2-dim block (<= 4 sqrt2) + one 1-dim block (2)), and
    the white-noise CHSH value of any lifting is in [-2/9, 2/9]; hence threshold >= 4/(3 sqrt2 + 1)
    = 0.762974 > 0.696152. So DKZ measurements realise the highest noise resistance for d = 3 (given the
    cited facet classification). General d: needs the facet list (Problem 27A) -- not shown.
(E) KL / statistical strength (vDGG): NUMERICAL only. d=3: DKZ + optimised symmetric state 0.076850 bits
    (DKZ + max-ent 0.057783); d=4: 0.097652 bits (max-ent 0.062380); local optimisation over all PVMs and
    pure states from perturbed starts never exceeds these (local-optimality evidence).  No proof.
(F) Evidence toward all d (NUMERICAL): tracial level 2 tight for d<=12 (d=11: -1.4e-9, d=12: -9.3e-9, gap 5e-9); DKZ a non-degenerate local optimum
    on |Phi_d> for d<=30 with Hessian gap 2 sec(pi/(2d))/d; I_ME(d) -> 32G/pi^2 (G Catalan).
    Obstacle to an all-d SOS: cross-power coupling is essential; adjacent-pair word set fails at d=10.

## 10. How to re-check
  python verify_tracial.py 3 4 5 6 7 8 9     (max-ent bound; exact)
  python verify_rigidity.py 3 4 5 6 7 8 9    (support conditions for rigidity)
  python verify_bip.py 3 4 5 6 7             (all-states Tsirelson bound; exact)
  python sanity_tracial.py <d> <D> ; sanity_bip.py <d> <D>   (numerical end-to-end)
Certificate files: cert_tracial_d{3..9}.pkl, cert_bip_d{3..7}.pkl (pickled Python ints/tuples:
block elements (words + cyclotomic coefficients), exact kernel bases N_b, exact Hermitian Gram Y_b,
exact bound lam; bipartite files also carry the monic polynomial p of lambda* over Q(zeta_{4d})).
Builders: exact_tracial.py / cert_tracial.py (mod), exact_bip.py / cert_bip.py, with cyclo.py,
extfield.py, modsolve.py/modsolve2.py (multi-modular exact solver), ipm.py (own interior-point SDP).

- verify_tracial.py 9 -> VERIFIED (497 classes, 131713 terms, min certified pivot 5.87e-5).
  Max-ent optimality + rigidity of DKZ now proven for d = 3..9 (I_ME(9) = 2.93650953168142802141797546004).
- Full level-2 tracial relaxation at d=11 (nvar 1510, 1650 s): tight to -1.4e-9 (gap 8e-10).
- Certificate sizes (SOS blocks / total Gram rank = number of independent SOS polynomials):
  max-ent  d=3..9: 12/48, 16/105, 20/184, 24/285, 28/264, 32/357, 36/464 (d>=7: adjacent word set);
  all-state d=3..8: 12/16, 16/33, 20/56, 24/85, 28/120, 32/161.
- Exact values d=9 (exact_values_9.txt): s_ME minimal polynomial of degree 6 (see table);
  lambda_max(K_9) of degree 15.
- facets_223.py: direct qhull enumeration of the (2,2,3) local polytope was stopped (memory > 2 GB,
  no output after ~15 min). The d=3 all-inequality noise statement in (D) therefore rests on the
  published facet classification (positivity, lifted CHSH, CGLMP_3), not on our own enumeration.
- verify_bip.py was refactored to allow parallel per-block work (-jN); the serial original is kept as
  verify_bip_serial.py. Re-verification with the current version: d=3,4 (-j4) and d=5,6,7 (-j8) all
  VERIFIED again (verify_bip_par_5to7.log).
- bipartite d=8: modsolve3 (parallel multi-modular solver, 8 workers) succeeded with 505 primes
  (15655-bit rationals); exact equations satisfied; exact LDL in progress.
- Corollary of (B): if d does not divide D, the max-ent value on |Phi_D> is STRICTLY below I_ME(d).
  Illustration (maxent_primal.py, d=3, labels k mod 3): best S found 1.34530 (D=4), 1.22771 (D=5),
  1.12707 = S_ME (D=6, DKZ (x) 1_2)  [S_ME(3) = 1.12706595].
- sparsify3.py: for d=3 a max-ent certificate supported on only 5 of the 12 symmetry blocks exists
  (blocks (0,0),(0,2),(0,3),(2,2),(2,3); kernel dims 2,4,4,3,1) -- numerical feasibility.
- NPA 1+AB numerics (npa_bip.py + ipm.py): d=8 diff -1.2e-8, d=9 diff +1.1e-10 -> tight numerically
  for d = 2..9 (ipm_bip_89.log).
- Full level-2 tracial relaxation at d=12 (nvar 1985, 2355 s): -9.3e-9 (IPM gap 5.1e-9) -> tight to solver accuracy.
- Checker speed-up: verify_tracial.CF now uses gmpy2.mpq (exact rationals; falls back to fractions).
  FINAL re-verification of all final certificate files with the final checker code:
  FINAL_verify_tracial.log (d=3..9 VERIFIED), FINAL_verify_bip.log (d=3..7 VERIFIED),
  FINAL_verify_rigidity.log (d=3..9 support conditions verified).
- bipartite d=8: cert_bip_d8.pkl (8.4 MB; builder 10550 s incl. parallel multi-modular solve with 505
  primes / 15655-bit rationals; min LDL pivot 7.1e-8) VERIFIED by verify_bip.py -j10 (132 classes,
  25145 expanded terms): I_8* = 2(lambda_max(K_8)-1)/7 = 3.10128058790546856685503528485,
  lambda_max(K_8) algebraic of degree 32.  => all-states Tsirelson bound proven for d = 3..8.
- d=10 max-ent exact certificate: exact data built (40 blocks, nvar 1118, kernels 16..27; reduced SDP
  bound 4.77089866014 vs S_DKZ 4.77089865933), but the 1118x1118 multi-modular solve (16 embeddings,
  ~6 min per prime per worker with the current numpy elimination) was stopped after ~1 h to free the
  workstation. d=10 therefore remains NUMERICAL (level-2 tight to -1.6e-9). A BLAS-blocked modular LU
  (P < 2^22, float64 matmul) would make d=10..12 feasible; not done.

## 11. FINAL CHECK (2026-09-29 07:57)
All final certificate files re-verified with the final checker code:
  FINAL_verify_tracial.log   max-ent bound, d = 3,4,5,6,7,8,9: VERIFIED
  FINAL_verify_bip.log       all-states Tsirelson bound, d = 3,4,5,6,7,8: VERIFIED
  FINAL_verify_rigidity.log  rigidity support conditions, d = 3..9: verified
