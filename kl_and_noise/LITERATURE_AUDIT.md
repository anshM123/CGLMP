# Literature audit for OQP 27B (workstream D), 2026-09-29

Sources were read in full where marked (PDF text extracted locally); others from abstracts or targeted
sections. Citation lists were taken from the Semantic Scholar API (citations of arXiv:2112.10803,
quant-ph/0506225 and 1607.04578). Earlier novelty check: ../../oqp27B/LOG.md section 1, which covers the
arXiv full-text search for "CGLMP" from 2015 to 2026.

## A. Origin and exact meaning of the problem

| Ref | Date | Exact content relevant to 27B | Read |
|---|---|---|---|
| R. Gill, "Better Bell inequalities (passion at a distance)", IMS Lecture Notes-Monograph Series 55, 135-148 (2007), arXiv:math/0610115 (OQP ref [4]) | Oct 2006 / Sep 2007 | Conjecture (a): all faces of the 2x2xr polytope are nonnegativity or (lifted) CGLMP. Gill writes "certainly true for r = 2,3,4 and 5", which is false for r = 4 (BGP 2010). Conjecture (b), the claim behind 27B: "numerical search for optimal experiments using the maximally entangled state has only uncovered the CGLMP measurements ... true both using Euclidean and relative entropy distances". Here noise resistance means mixing q with UNIFORM outcomes, and violation means violation of local realism, i.e. any facet. He also asks "why are the CGLMP measurements optimal for the CGLMP inequality?". | full |
| van Dam, Gruenwald, Gill, IEEE TIT 51, 2812 (2005), quant-ph/0307125 (OQP ref [3]) | Jul 2003 / Jun 2004 | Defines S^UNI <= S^UC <= S^COR (sup over setting distributions of inf over LHV of KL). CHSH gives 0.0462738469 bits. Conjecture 4: the strongest 2x2x3 proof has strength 0.077, with state ~0.6475/0.6475/0.4019. The max-ent state gives only 0.058. | full (relevant sections) |
| Acin, Gill, Gisin, PRL 95, 210402 (2005), quant-ph/0506225 (OQP ref [5]) | Jun 2005 | Numerical, rank-one PVMs, setting distribution free (found uniform). d=3: DKZ is optimal for Phi_3 (0.058 bits). Global optimum: DKZ with delta ~ 0.642 (0.077 bits). d=4: "optimal measurements are those maximizing the Bell violation for Phi_4", 0.098 bits. For d > 4 they ASSUME DKZ measurements and a dual inequality of Gill/CGLMP form. | full |

## B. Max-ent optimality, Tsirelson bounds, SOS

| Ref | Date | Statement | Rigor |
|---|---|---|---|
| Kaszlikowski, Gnacinski, Zukowski, Miklaszewski, Zeilinger, PRL 85, 4418 | 2000 | max-ent qudits, unbiased multiports (DFT + phases), LP against all LHV models; thresholds decrease with d | numerical, restricted family |
| Durt, Kaszlikowski, Zukowski, PRA 64, 024101 (OQP [2]) | 2001 | the specific phases (DKZ) up to N = 16 | numerical |
| Collins, Gisin, Linden, Massar, Popescu, PRL 88, 040404 (OQP [1]) | 2002 | CGLMP inequality; reproduces the thresholds of [2] | analytic inequality; optimality numerical |
| Acin, Durt, Gisin, Latorre, PRA 65, 052325 | 2002 | non-max-ent state maximises CGLMP (d=3: 1+sqrt(11/3)); more noise robust; remark on rank-2 states | numerical (d=3 value analytic for DKZ) |
| Navascues, Pironio, Acin, PRL 98, 010401 | 2007 | NPA reproduces the ADGL values | numerical |
| Zohren, Gill, PRL 100, 120406 (quant-ph/0612020) | 2006/2008 | large d; "best measurements are the same as conjectured best measurements for the maximally entangled state" | numerical |
| Zohren, Reska, Gill, Westra, EPL 90, 10002 (1003.0616) | 2010 | d -> infinity limit | analytic (limit only) |
| Lang, Vertesi, Navascues, J. Phys. A 47, 424029 (1402.2850) (OQP [7]) | 2014 | tracial relaxation Q^2_+ of max-ent correlations of ANY dimension, with POVMs (positivity constraints, not projector constraints); CGLMP_3: DKZ value to 1e-10 | numerical, d = 3 only |
| Ioannou, Rosset, arXiv:2112.10803 (OQP [10]) | Dec 2021 | analytic SOS for the CGLMP Tsirelson bound (all states), d = 3,4 | analytic |
| Citations of 2112.10803 (S2 API, 11 entries incl. 2601.02893, 2606.19033, 2608.16512, 2609.05555) | to Sep 2026 | none gives an exact CGLMP value for d >= 5 or max-ent optimality | -- |
| Coccia, Padovan, Vallone, arXiv:2606.21626 | Jun 2026 | systematic SOS Tsirelson bounds for max-ent states; generalised SATWAP (2 inputs) built on the CGLMP measurements; the CGLMP functional itself is NOT treated | analytic, different expressions |

## C. Self-testing and rigidity

| Ref | Date | Statement |
|---|---|---|
| Yang, Vertesi, Bancal, Scarani, Navascues, PRL 113, 040401; Bancal et al., PRA 91, 022115 | 2014/2015 | self-test of the (non-max-ent) CGLMP_3 optimum, numerical SWAP method |
| Salavrakos, Augusiak, Tura, Wittek, Acin, Pironio, PRL 119, 040402 (1607.04578) | 2016/2017 | SATWAP: Bell expressions (m settings, d outcomes) maximised by Phi_d with CGLMP measurements; SOS for all d |
| Kaniewski, Supic, Tura, Baccari, Salavrakos, Augusiak, Quantum 3, 198 (1807.03332) | 2018/2019 | d-setting MUB inequalities; self-test of Phi_3 with 3 MUBs |
| Sarkar, Saha, Kaniewski, Augusiak, npj QI 7, 151 (1909.12722) | 2019/2021 | DEVICE-INDEPENDENT self-test of Phi_d and the CGLMP measurements (2 settings, all d) via SATWAP (m=2) |
| Pandya, Sarkar, Augusiak, arXiv:2602.08469 | Feb 2026 | (2,2,3) "CGLMP-style" SOS inequalities self-testing Phi_3 and a class of measurements; SATWAP is a special case; not CGLMP itself |

No paper was found that proves CGLMP_d optimality on maximally entangled states for any d exactly, or rigidity within max-ent strategies for the CGLMP functional. That includes 2024-2026 arXiv, the citations of 2112.10803 and the citations of 1607.04578 with CGLMP / qudit / max-ent / SOS / self-testing in the title.

## D. Facets and noise robustness

| Ref | Date | Statement |
|---|---|---|
| Masanes, QIC 3, 345 (OQP [6]) | 2003 | CGLMP is a facet |
| Collins, Gisin, J. Phys. A 37, 1775 | 2004 | (2,2,3): only positivity, CHSH (lifted), CGLMP |
| Bancal, Gisin, Pironio, J. Phys. A 43, 385303 (1004.4146) | 2010 | (2,2,4): 8 symmetric facet classes (S1 = CGLMP, 7 new), which answers Gill's conjecture (a), i.e. 27A, negatively |
| Gruca, Laskowski, Zukowski, PRA 85, 022118 (1111.3955) | 2012 | exhaustive LP numerics, 2 settings, unrestricted U(d): max-ent thresholds 0.6962 / 0.6906 / 0.6871 (d = 3,4,5), equal to DKZ. Schmidt-rank-2 states are MORE robust: 0.6821 / 0.6442 / 0.6071 (the "embedding" effect; see also ADGL 2002) |
| Jesus, Zambrini Cruzeiro, PRA 108, 052220 (2212.03212) | 2022/2023 | (2,2,4,4): 34 classes (10 non-lifted) from polytope slicing, completeness NOT proven; per-facet noise with the state maximally violating that facet (not max-ent). No ancillary data on arXiv |

## E. Statistical strength follow-ups
Citations of AGG 2005 (Semantic Scholar, 81 entries) include 1108.2468 ("Asymptotically optimal data analysis for rejecting local realism", test factors / PBR), 1510.07233 ("(Nearly) optimal P values for all Bell inequalities"), 1812.06236 (hypothesis testing of physical theories), 1805.09451 (Fonseca et al., probability of violation), 1802.09982 (Lipinska et al., random measurements), and a 2011 item without arXiv id, "On the optimization of nonlocality proofs in quantum statistics", which could not be located. None of the located papers revisits the KL optimisation of 2 x 2 x d tests. We found no later work on KL-optimal measurements for CGLMP / 2 x 2 x d beyond vDGG 2005, AGG 2005 and Gill 2007.

## F. Classification of our results

| Result | Status | Classification |
|---|---|---|
| Max-ent optimality of DKZ for CGLMP_d (any local dimension, PVMs), d = 3..9, exact SOS certificates (../../oqp27B) | proved (computer-assisted, exact) | APPARENTLY NEW as a proof. Known numerically (LVN 2014 for d = 3 including POVMs; Gill/AGG/ZG numerics). Scope: PVMs. Our POVM seesaw (d = 3,4) finds no improvement. |
| Rigidity within max-ent strategies (DKZ (x) 1 unique), d = 3..9; algebraic reduction for all d (C_rigidity) | proved | APPARENTLY NEW for CGLMP. KNOWN UNDER DIFFERENT ASSUMPTIONS: DI self-testing of Phi_d + CGLMP measurements via SATWAP (Sarkar et al. 2021, all d; Kaniewski et al. 2019 d=3 via MUB inequalities) |
| All-state Tsirelson bounds I*(d) = 2(lambda_max(K_d)-1)/(d-1), d = 5..8 exact | proved (computer-assisted) | NEW for d >= 5 (d = 3,4: Ioannou-Rosset 2021; values numerically known since ADGL 2002 / NPA 2007) |
| Noise (N) on Phi_d against ALL Bell inequalities | numerical d = 3,4,5 (ours, and 34 + BGP facet checks at d=4); proof d = 3 | numerics KNOWN (GLZ 2012); the d = 3 PROOF (facet list + max-ent theorem + lifted-CHSH bound) is new but elementary given the max-ent theorem |
| CGLMP-witness noise corollary (v_c >= 2/I_ME, equality iff DKZ) | proved d = 3..9, conditional all d | new corollary, trivial given the max-ent theorem and rigidity |
| (N) with the state also optimised is FALSE for d >= 3 (Schmidt-rank-2 embedding) | rigorous (explicit, LP-checked) | KNOWN (ADGL 2002, GLZ 2012) |
| KL clause on Phi_d FALSE for d = 4..9 (certified intervals) | rigorous (exact rational certificates + interval arithmetic) | APPARENTLY NEW. It contradicts Gill 2007's numerical claim. AGG 2005 made no max-ent claim for d >= 4 |
| KL clause with the state optimised (AGG's global problem) FALSE for d = 4, 6, 7, 8 | numerical (competitor side certified; DKZ side = best of many state starts) | APPARENTLY NEW. It contradicts AGG 2005's d = 4 claim, by 6.2e-5 bits |
| DKZ remains KL-optimal on Phi_3, and DKZ + state remains optimal at d = 3,5 | numerical | consistent with AGG and vDGG Conjecture 4 (0.0768501461 bits, Schmidt 0.64747/0.64747/0.40194) |
| AGG Gill-form dual ansatz exact only for d <= 5 at the max-ent DKZ point | numerical (1e-14) | APPARENTLY NEW observation |
| DKZ is a non-degenerate local maximum of KL in the Fourier-diagonal family, d = 3..6 | numerical (Hessian) | new, numerical |
