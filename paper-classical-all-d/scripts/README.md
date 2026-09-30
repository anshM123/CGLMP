# Numerical cross-checks

These scripts check the formulas and inequalities of the paper numerically. They are **not** part of any proof: the
proofs in the paper are complete without computer assistance. Everything is floating point, except
`check_constants.py`, which uses interval arithmetic.

Requirements: Python 3 with `numpy` and `scipy`, plus `mpmath` for `check_constants.py`. The recorded runs used
Python 3.14.3, numpy 2.4.6, scipy 1.17.1 and mpmath 1.3.0. Every script runs in a few seconds.

Run everything from this folder:

```bash
bash run_checks.sh                      # or: PYTHON=/path/to/python bash run_checks.sh
```

This writes one log per script to `../logs/`. Each log starts with the command and the interpreter version. The
individual commands are:

```bash
python verify_writeup.py 9
python verify_junction.py 6 7 8 12 16 24 40 64 101
python verify_chain.py 24,6,6,500 40,10,10,500 48,12,12,300 80,20,20,100 30,4,11,500 50,7,19,300
python cot_rearr.py 12,3,3 16,4,4 20,5,5 24,6,6 16,3,6 20,4,7
python swap_test.py 9,3,3 10,3,3 12,3,4 12,4,4
python classical_sa.py 12
python check_reduction.py
python check_constants.py
```

## What each script checks

Notation as in the paper: `N = 4d` in the application, `k(u) = cot(pi u/N)` (called kappa in the paper),
`<A,B> = sum_{x in A, y in B} k(x-y)`, and `C` is the complement of `A u B`.

| Script | Checks | Paper |
|---|---|---|
| `verify_writeup.py DMAX` | **(T)** For d = 2..DMAX it enumerates all 4^d classical configurations of the clock model. It checks that max F = F_DKZ with exactly 4d maximisers, all of them intervals, and the identity F(S) = <S,S-d>/2 - d/2. **(L1)** cyclic symmetry, **(L2)** exchange identity, **(L3)** smearing identity Q = Q_c + E, each on random colourings. **(L5)** the partial-fraction formula e(m) = (N/pi) sum_j eps(m+jN). | Theorem B, Lemma 3.1, Lemmas 4.1, 4.2, 4.4, Lemma 4.9(2) |
| `verify_junction.py N...` | For each N: signs e(m) <= 0 on [1, N/2], the lower bound (N/pi) eps(m) <= e(m), the bound on \|e(1)\|, and the sizes of sc and fw in units of N/pi. Also the margin -\|e(1)\| + 2 sc + fw for winding 2. | Lemma 4.9(3),(4), Lemma 4.10 |
| `verify_chain.py N,a,b,trials...` | Random 3-colourings are driven by random out-of-order exchanges until they are cyclically ordered. Each exchange is checked against the exchange identity. At the end point it checks Q_c <= Q_c(block), E <= E(block) and <A,B> <= Phi_N(a,b). | Lemma 4.2, Proposition 4.6, Lemma 4.10, Theorem A |
| `cot_rearr.py N,a,b...` | Exhaustive maximum of <A,B> over all disjoint A, B of sizes a, b. For each B the best A is its top-a set, by the bathtub principle. Also checks the identity F(S) = <S,S-d>/2 - d/2 on random transversals. | Theorem A, Lemma 3.1 |
| `swap_test.py N,a,b...` | Enumerates all colourings and all out-of-order adjacent pairs (A,B), (B,C), (C,A). Checks that no exchange decreases <A,B>. | Lemma 4.2, Corollary 4.3 |
| `classical_sa.py DMAX` | Sherali-Adams linear programs for the classical clock model (HiGHS solver), with pair marginals and with three-site marginals of the differences. | Section 6.2 |
| `check_reduction.py` | Computes CGLMP values from outcome probabilities with the expression of Collins et al. **(R1)** the chain form I_d = 4 - 2S/(d-1). **(R2)** the trace formula for S on random projective strategies on max-ent states. **(R3)** I_d = 4F(V)/(d(d-1)) for covariant strategies with random non-commuting V. **(R4)** twirl plus Stone-von Neumann reduction on random strategies. **(R5)** equal links <=> commuting reduced family. **(R6)** all clock-phase strategies for d <= 7. **(R7)** the CGLMP/DKZ measurements give I_ME(d) and equal G R^(0) G^*. **(R8)** projection form. **(R9)** Ising form. | Section 2, Appendix A, Section 6 |
| `check_constants.py` | Interval-arithmetic enclosures (60 digits) of the constants in Lemmas 4.9 and 4.10: \|eps(1)\| - \|eps(5)\| >= 0.38493; the tail sum < 0.0773 < 4/45; c_sc <= 0.07785 and c_fw <= 0.11481; the final margin 0.49974 - 0.30708 w < 0 for w >= 2. Also checks the series for eps against its closed form. | Lemmas 4.9, 4.10 |

`verify_writeup.py` and `verify_chain.py` import the functions `kbar` and `eps` from `verify_junction.py`. The
other scripts are self-contained.

## Outcomes of the recorded runs (`../logs/`)

- `verify_writeup.log`
  - d = 2..9: max F = F_DKZ to 3e-14, exactly 4d maximisers, all intervals.
  - Lemma 3.1 identity error <= 8e-14.
  - Lemma identities to <= 4e-13; partial-fraction formula to <= 2e-7 (quadrature accuracy).
- `verify_junction.log`: all signs and bounds hold for N = 6..101.
  - \|e(1)\|/(N/pi) lies in 0.38536..0.38630.
  - sc/(N/pi) <= 0.0733 and fw/(N/pi) <= 0.1097.
  - The winding-2 margin is between -0.30 and -0.13 (units of N/pi).
- `verify_chain.log`
  - The exchange identity holds to < 3e-13.
  - At every end point Q_c - Q_c(block) <= 0, E - E(block) <= 0 and Q - Q_block <= 1e-14. Here 1e-14 is rounding, for flows that end at the block itself.
- `cot_rearr.log`: in all six cases the exhaustive maximum equals Phi_N(a,b) and is attained by adjacent intervals.
- `swap_test.log`: 290,710 out-of-order exchanges; none decreases <A,B>.
- `classical_sa.log`: the three-site relaxation is exact for d = 3..8 but exceeds F_DKZ for d = 9..12, by +0.535, +2.52, +5.23 and +8.53.
- `check_reduction.log`
  - All identities hold to <= 2e-14.
  - Random strategies have unequal links (difference >= 1.29) and non-commuting reduced families (commutator >= 1.63).
  - The clock-phase exhaustion for d = 2..7 gives max I_d = I_ME(d), attained by exactly 4d one-step configurations.
- `check_constants.log`: all stated roundings are valid.
