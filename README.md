# CGLMP on maximally entangled states: exact optimality, rigidity and Tsirelson bounds

**Authors:** Ansh Mishra, Aryan Senthilkumar. **License:** MIT.

This repository contains computer-assisted proofs for the CGLMP Bell inequality, together with the certificates, independent checkers and verification logs. The work addresses part B of IQOQI Vienna Open Quantum Problem 27, [*The power of CGLMP inequalities*](https://oqp.iqoqi.oeaw.ac.at/the-power-of-cglmp-inequalities). The live site was unavailable when we checked; see the [archived copy](http://web.archive.org/web/20231029013544/https://oqp.iqoqi.oeaw.ac.at/the-power-of-cglmp-inequalities).

Every result stated as proved rests on exact data. Numbers are elements of Q(ζ₄d) or Q, positivity is certified by exact LDLᵀ pivots, and the KL bounds use interval arithmetic. The checkers use no floating-point solver, and they are independent of the programs that built the certificates. Results that are only numerical are labelled as such.

## Results

Notation:
- I_d is the CGLMP value.
- DKZ is the Durt–Kaszlikowski–Zukowski family of measurements: computational-basis measurements transformed by the discrete Fourier transform and diagonal phases.
- Φ_D is the maximally entangled state of local dimension D.
- The max-ent value of the DKZ strategy is

  I_ME(d) = 4/(d(d−1)) · Σ_{j=1}^{d−1} (d−j) sec(πj/(2d)).

### 1. DKZ is optimal on maximally entangled states (d = 3, …, 12)

For every local dimension D and all projective measurements on Φ_D, the CGLMP_d value is at most I_ME(d). DKZ attains this value.
- **Proof:** an exact sum-of-squares certificate in the tracial (maximally entangled) relaxation, one for each d.
- **Data:** `certificates/maxent/cert_tracial_d{3..12}.pkl`.
- **Checker:** `verify/verify_tracial.py`.
- **Values:** I_ME = 2.8729340512 (d=3), 2.8962432185 (d=4), 2.9105448081 (d=5), …, 2.9447540269 (d=12). The full digits are in the logs.

### 2. Rigidity and self-testing within maximally entangled strategies (d = 3, …, 12)

Suppose projective measurements on Φ_D attain I_ME(d). Then d divides D, and up to a local unitary the strategy is DKZ ⊗ 1 on C^d ⊗ C^{D/d}.
- **Certificate side:** the support conditions on the certificate, `verify/check_support.py` and `verify/verify_rigidity.py`.
- **Algebraic side:** holds for all d. It consists of a Jordan-type lemma, an orthogonality induction and imprimitivity.
  - It is written up in `rigidity/RIGIDITY.md`.
  - Its core is machine-checked in Lean 4 with Mathlib in `lean/CGLMPRigidity/`, using only the standard axioms.
- The write-up states "unconditional for d = 3..9". The d = 10, 11, 12 certificates were verified afterwards, see `logs/independent_verify_27B_d10to12.log`.

### 3. Exact Tsirelson bounds over all states (d = 3, …, 8)

The maximum quantum CGLMP value over all states and measurements is

I*(d) = 2(λ_max(K_d) − 1)/(d − 1),  where K_d = [sec(π(k−l)/(2d))]_{k,l}.

- **Data and checker:** `certificates/allstates/cert_bip_d{3..8}.pkl`, checked by `verify/verify_bip.py`.
- **Prior work:** d = 3, 4 were known (Ioannou–Rosset, 2021). As far as we found, the exact values for d = 5..8 are new.

### 4. The Kullback–Leibler clause is false for d = 4, …, 9

On Φ_d, explicit competitor strategies give certified statistical strength (KL divergence to the local set) strictly larger than the DKZ measurements.
- **Certification:** exact rational certificates plus interval arithmetic, in `kl_and_noise/`, log `kl_and_noise/FINAL_certificates.log`.
- **Example:** at d = 4 the competitor exceeds DKZ by at least 0.0255 bits.
- **d = 3:** our numerics agree that DKZ is KL-optimal.
- **Prior work:** a d = 4 counterexample was found independently, and published earlier, by Y. Zhang, *CGLMP measurements need not maximize statistical strength* (Zenodo, 28 September 2026, [doi:10.5281/zenodo.23022433](https://doi.org/10.5281/zenodo.23022433)), with a separation above 0.022 bits. Our d = 4 certificate is a second, independent one. As far as we found, the counterexamples for d = 5..9 are new.

### 5. Noise robustness

- **CGLMP witness:** the critical visibility of Φ_d, detected with the CGLMP inequality, is at least 2/I_ME(d), with equality iff DKZ. This follows from result 1 and holds for d = 3..12.
- **All Bell inequalities:** proved for d = 3 using the complete (2,2,3) facet list. It is numerical only for d = 4, 5.

### 6. Reductions toward all d

`reductions/` and `research-log/reductions_LOG.md` contain two exact reductions of the all-d problem:
- the clock-model reduction;
- the Fourier multiplier Ĵ(l) = χ₄(l)(d − l/2) − 1/2 of the classical Ising-ring picture.

They also include local-optimality checks.

### 7. Every d, commuting (equal-link) strategies: analytic proof

Paper: [`paper-classical-all-d/main.pdf`](paper-classical-all-d/main.pdf), *A rearrangement inequality for the discrete Hilbert transform, with an application to optimal CGLMP measurements* (18 pages). The proofs are by hand, with no computer assistance.

- **Theorem A (discrete Hilbert-transform rearrangement).** Let A and B be disjoint subsets of Z_N with |A| = a and |B| = b. Then Σ_{x∈A, y∈B} cot(π(x−y)/N) is largest when A and B are adjacent intervals, with B immediately before A. When a, b ≥ 1 and a + b < N, only rotations of that pair attain the maximum.
  - The proof has three ingredients: an exact exchange identity; the circle version of Stein–Weiss/Laeng universality combined with the bathtub principle; and an explicit estimate of the junction terms.
  - We have not found Theorem A, or an equivalent statement, in the literature.
- **Theorem B (clock model, commutative sector, every d ≥ 2).** F(V) ≤ F_DKZ for pairwise commuting unitaries V_0..V_{d−1} with V_k^4 = 1. Equality holds exactly for the one-step configurations.
- **Corollary C (physics).** Take any projective strategy on a maximally entangled state, of any local dimension, whose four link operators coincide. Its reduced family then commutes.
  - It satisfies I_d ≤ I_ME(d).
  - Equality holds if and only if the strategy is DKZ ⊗ 1, up to a local unitary u ⊗ ū.
  - This covers, among others, all 4^d clock-phase strategies (Fourier bases dressed by diagonal phases) and their direct sums.
- **Theorem D (a non-commuting class, every d).** F(V) ≤ F_DKZ also holds when each unitary switches between the values of two optimal configurations, V_k = i^{−a_k}(1 − P_k) + i^{−a'_k}P_k, with arbitrary, non-commuting projections P_k.
  - Equality holds only when the P_k are nested.
  - The proof: on this class F is exactly affine in the overlaps τ(P_j P_k), with coefficients at least 2 sec(π/2d) > 0.
- **Checks.** `paper-classical-all-d/scripts/` holds the numerical cross-checks (interval arithmetic for all constants), with outputs in `paper-classical-all-d/logs/`. They are not part of the proofs.

## Status

The problem statement asks for more than is proved here:

| Claim of OQP 27B | Status here |
|---|---|
| DKZ optimal on maximally entangled states | proved for d = 3..12 (projective measurements); proved for **every d** among equal-link (commuting) strategies and among the non-commuting two-window families (result 7); open in general for d ≥ 13 |
| DKZ unique (rigidity) | proved for d = 3..12; for every d among equal-link strategies (result 7); algebraic part proved for all d |
| Tsirelson bound, all states | exact for d = 3..8 |
| KL optimality of DKZ | **false** for d = 4..9 (d = 4 first shown by Y. Zhang, 2026) |
| Noise robustness against all Bell inequalities | proved for d = 3; numerical for d = 4, 5 |

Part 27A (whether all facets of the (2,2,d) local polytope are CGLMP-type) was already answered negatively by Bancal, Gisin and Pironio (J. Phys. A 43, 385303 (2010)). The literature check is `kl_and_noise/LITERATURE_AUDIT.md`.

## Layout

| Path | Contents |
|---|---|
| `certificates/maxent/` | exact max-ent SOS certificates, d = 3..12 |
| `certificates/allstates/` | exact all-state SOS certificates, d = 3..8 |
| `verify/` | independent exact checkers |
| `builders/` | the programs that produced the certificates (not needed for checking) |
| `values/` | exact and high-precision values of I_ME(d) and I*(d) |
| `rigidity/` | the rigidity proof and its support checks |
| `lean/` | Lean 4 formalisation of the algebraic core of rigidity |
| `kl_and_noise/` | KL counterexample certificates, noise checks, literature audit |
| `reductions/` | exact reductions and local-optimality checks toward all d |
| `paper-classical-all-d/` | paper for result 7 (LaTeX source and PDF), numerical cross-check scripts and logs |
| `research-log/` | dated working logs |
| `logs/` | verification runs, including the independent re-checks |

## How to check

Requirements:
- Python 3 with `mpmath`; `gmpy2` is optional and makes the checks faster.
- The Lean part needs Lean 4 v4.33.1 (see `lean/lean-toolchain`) and Mathlib at the pinned revision in `lean/lakefile.toml`.

The max-ent certificates must be checked from their own folder:

```bash
cd certificates/maxent && python ../../verify/verify_tracial.py 3 4 5 6 7 8 9 10 11 12
```

Each d ends with the line `VERIFIED: max CGLMP value over maximally entangled states (any dimension) = I_ME`. Run times are seconds for small d and about three hours each for d = 11, 12.

The all-state Tsirelson certificates:

```bash
cd certificates/allstates && python ../../verify/verify_bip.py 3 4 5 6 7 8
```

The rigidity support conditions, run from the repository root:

```bash
python verify/check_support.py 3 4 5 6 7 8 9
```

The Lean core of rigidity, run from `lean/`:

```bash
cd lean && lake exe cache get && bash CGLMPRigidity/check.sh
```

## Citation

See `CITATION.cff`. Please cite the repository and the authors above.
