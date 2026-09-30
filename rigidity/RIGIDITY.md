# Rigidity of the DKZ strategy for CGLMP on maximally entangled states

IQOQI Open Quantum Problem 27B, first sentence: *the optimal CGLMP observables on a maximally entangled
state are necessarily of the Durt-Kaszlikowski-Zukowski (DKZ) form* (computational-basis measurements
transformed only by the discrete Fourier transform and diagonal unitaries).

Workstream C, 2026-09-29. Status: **complete proof, conditional on an explicitly specified certificate
property (holds for every d for which such a certificate exists); unconditional for d = 3, 4, 5, 6, 7, 8, 9.**
The algebraic core (Section 3, Proposition 2.5, $d\mid D$, the core of Lemma 7.1, and non-vacuity) is
machine-checked in Lean 4 + Mathlib for all $d$ (Section 10); the certificate side is checked by two independent
exact programs (Section 6).

Contents
0. Results
1. Setting, conventions, the DKZ strategy
2. From a certificate to equality conditions (faithful-trace step)
3. Algebraic rigidity (Jordan-type lemma, orthogonality induction, imprimitivity)
4. The self-testing theorem
5. Consequence: d not dividing D gives a strictly smaller value
6. Verified instances d = 3..9
7. What an all-d certificate must provide
8. Scope and limitations
9. Changes with respect to the original sketch (lab LOG section 6)
10. Formalisation status

---------------------------------------------------------------------------------------------------------

## 0. Results

Throughout, $d\ge 3$, $z:=e^{2\pi i/(4d)}$, $w:=z^4=e^{2\pi i/d}$ (so $z^d=i$), $m_d:=\lfloor d/2\rfloor+1$.
A *strategy of local dimension D* is a 4-tuple $R=(R_1,R_2,R_3,R_4)$ of $D\times D$ unitaries with
$R_k^d=1$; it encodes projective measurements on $|\Phi_D\rangle$ via $R_1=A_1,\ R_2=B_1^{T},\ R_3=A_2,\
R_4=B_2^{T}$ (Section 1). $S(R)$ is the CGLMP functional in Gill form (CGLMP value $I_d=4-2S/(d-1)$;
maximising $I_d$ = minimising $S$). The four *links* are
$$L^{(0)}_n=R_1^nR_2^{-n},\quad L^{(1)}_n=R_2^nR_3^{-n},\quad L^{(2)}_n=R_3^nR_4^{-n},\quad L^{(3)}_n=w^{-n}R_4^nR_1^{-n}.$$

**(OPT_d)** $S(R)\ge S(R^{\rm DKZ})$ for every strategy of every local dimension $D$.

**(RIG_d)** every strategy $R$ (any $D$) with $S(R)=S(R^{\rm DKZ})$ satisfies

* (E_1') $R_1R_2^{-1}=R_2R_3^{-1}=R_3R_4^{-1}$, and
* (T_n) $(U_n-z^{-n})(U_n-iz^{-n})=0$ for $U_n:=R_1^nR_2^{-n}$ and $n=1,\dots,m_d$.

**Theorem A (algebraic rigidity; Theorems 3.8, 3.9).** Let $R_1,\dots,R_4$ be unitaries with $R_1^d=R_2^d=1$
in a unital $*$-algebra in which $x^*x=0\Rightarrow x=0$ (e.g. $M_D(\mathbb C)$, any C*-algebra), satisfying
(E_1') and (T_n) for $n=1,\dots,m_d$. Then there is a unital $*$-homomorphism $\pi:M_d(\mathbb C)\to\mathcal A$
with $R_k=\pi(R^{c}_k)$, $k=1,\dots,4$, where $R^c$ is the canonical DKZ tuple (1.4). If $\mathcal A=M_D(\mathbb C)$,
then $d\mid D$ and $R_k=W(R^c_k\otimes 1_M)W^*$ for a unitary $W$, $M=D/d$.

**Theorem B (certificate $\Rightarrow$ OPT and RIG; Theorem 2.6).** If a symmetric SOS certificate for
$S\ge S(R^{\rm DKZ})$ exists whose support (Definition 2.1) contains the chain modes $M_{n,t}$ and the
two-eigenvalue elements $T_n$ for $1\le n\le m_d$ (or, alternatively, $M_{1,t}$ and the first-link elements
$V_n$), then (OPT_d) and (RIG_d) hold.

**Theorem C (self-testing form; Theorem 4.1).** Assume (OPT_d) and (RIG_d). If projective measurements on
$|\Phi_D\rangle$ attain the maximal max-ent CGLMP value $I_{\rm ME}(d)$, then $d\mid D$ and there is a unitary
$V:\mathbb C^D\to\mathbb C^d\otimes\mathbb C^M$ with $(V\otimes\bar V)|\Phi_D\rangle=|\Phi_d\rangle\otimes|\Phi_M\rangle$,
$VA_x^aV^*=A_x^{a,\rm DKZ}\otimes 1_M$, $\bar VB_y^b\bar V^*=B_y^{b,\rm DKZ}\otimes 1_M$. $V$ is unique up to
$V\mapsto(1_d\otimes u)V$, $u\in U(M)$.

**Corollary D (Corollary 5.1).** Assume (OPT_d), (RIG_d). If $d\nmid D$, the maximum of $I_d$ over projective
measurements on $|\Phi_D\rangle$ is *strictly* smaller than $I_{\rm ME}(d)$ (and the gap tends to 0 along
$D=kd+r$, $k\to\infty$, so there is no uniform gap).

**Corollary E (Section 6).** For $d=3,4,5,6,7,8,9$ the certificates `cert_tracial_d{d}.pkl` satisfy the
hypotheses of Theorem B (exact, independently re-checked). Hence Theorems A-C and Corollary D hold
unconditionally for these $d$: *every* optimal maximally entangled strategy is DKZ $\otimes\,1_M$ up to local
unitaries; the DKZ form is necessary.

**Bridges (Section 7.1).** In the chain-mode coordinates $W(n,t)$ of workstream A, a strategy that is
*grade-sharp* ($W(n,1)=W(n,2)=0$) for $n\le\lfloor d/2\rfloor+1$ is DKZ$\,\otimes1_M$ (Lemma 7.1, all $d$); in the
reduced variables $V_k$ of workstream B, (E_1') is cyclic adjacent commutation of the $V_k$ and (T_n) are B's
relations $g_{jk}$ (Lemma 7.2). Also: face-reduced certificates with positive definite Gram blocks carry the
rigidity support conditions automatically (Lemma 1.9, Section 7).

---------------------------------------------------------------------------------------------------------

## 1. Setting, conventions, the DKZ strategy

### 1.1 The CGLMP functional

Outcomes $a,b\in\mathbb Z_d$, inputs $x,y\in\{1,2\}$; $X_x,Y_y$ denote the outcomes. For $k\in\mathbb Z$ let
$m(k)\in\{0,\dots,d-1\}$ be the residue of $k$ mod $d$. The Gill (chain) form is
$$S:=\mathbb E\,m(X_1-Y_1)+\mathbb E\,m(Y_1-X_2)+\mathbb E\,m(X_2-Y_2)+\mathbb E\,m(Y_2-X_1-1).$$
The CGLMP value is $I_d=4-2S/(d-1)$ (lab LOG sections 2-3; established there numerically on random boxes,
`conventions.py`, `conventions2.py`): the literal Collins-Gisin-Linden-Massar-Popescu expression coincides with
this after exchanging the two settings of *each* party ($X_1\leftrightarrow X_2$, $Y_1\leftrightarrow Y_2$). That
relabelling permutes the four measurements and is irrelevant for every statement below. Since
$S\mapsto 4-2S/(d-1)$ is strictly decreasing, *maximisers of $I_d$ are exactly minimisers of $S$*.

**Lemma 1.1 (Fourier expansion).** For $k\in\mathbb Z$: $m(k)=\frac{d-1}{2}+\sum_{n=1}^{d-1}c_n w^{nk}$ with
$c_n=-1/(1-w^{-n})$.

*Proof.* $m$ is $d$-periodic, so $m(k)=\sum_{n=0}^{d-1}\hat m_n w^{nk}$ with
$\hat m_n=\frac1d\sum_{j=0}^{d-1}j\,x^{j}$, $x=w^{-n}$. For $n=0$ this is $(d-1)/2$. For $x^d=1\ne x$:
$(1-x)\sum_{j=0}^{d-1}jx^j=\sum_{j=1}^{d-1}x^j-(d-1)x^d=-1-(d-1)=-d$, so $\hat m_n=-1/(1-x)$. $\square$

### 1.2 Projective strategies on the maximally entangled state

**Definition 1.2.** A *maximally entangled projective strategy of local dimension $D$* consists of
$H_A=H_B=\mathbb C^D$ with the computational basis $\{|k\rangle\}_{k=1}^D$, the state
$|\Phi_D\rangle=D^{-1/2}\sum_k|k\rangle|k\rangle$, and projection-valued measurements $\{A_x^a\}_{a\in\mathbb Z_d}$,
$\{B_y^b\}_{b\in\mathbb Z_d}$ ($A^a=(A^a)^*=(A^a)^2$, $\sum_aA^a=1$; some $A^a$ may be $0$), with
$p(a,b|x,y)=\langle\Phi_D|A_x^a\otimes B_y^b|\Phi_D\rangle$.

The generalised observable $A_x:=\sum_a w^aA_x^a$ is unitary with $A_x^d=1$, $A_x^n=\sum_aw^{na}A_x^a$ and
$A_x^a=\frac1d\sum_{n=0}^{d-1}w^{-na}A_x^n$. Conversely every unitary $U$ with $U^d=1$ is normal with spectrum in
$\{w^a\}$, and its spectral projections (spectral theorem for normal matrices) form a $d$-outcome PVM. So
PVMs $\leftrightarrow$ unitaries of order dividing $d$, bijectively.

**Lemma 1.3 (transposition convention).** For $X,Y\in M_D(\mathbb C)$:
$\langle\Phi_D|X\otimes Y|\Phi_D\rangle=\mathrm{tr}(XY^T)$, where $\mathrm{tr}=\frac1D\mathrm{Tr}$ and $Y^T$ is the
transpose *in the computational basis used to define $|\Phi_D\rangle$*.

*Proof.* $\frac1D\sum_{k,l}\langle k|X|l\rangle\langle k|Y|l\rangle=\frac1D\sum_{k,l}X_{kl}(Y^T)_{lk}$. $\square$

**Definition 1.4.** $R_1:=A_1,\ R_2:=B_1^T,\ R_3:=A_2,\ R_4:=B_2^T$. The map $(A_1,A_2,B_1,B_2)\mapsto R$ is a
bijection onto $\mathcal U_D:=\{R\in U(D)^4:\ R_k^d=1\}$ (transposition preserves unitarity and order and is
invertible). Elements of $\mathcal U_D$ are called *strategies of local dimension $D$*.

**Lemma 1.5.** $S=S(R):=2(d-1)+\sum_{n=1}^{d-1}c_n\,\mathrm{tr}\big(L_n^{(0)}+L_n^{(1)}+L_n^{(2)}+L_n^{(3)}\big)$.

*Proof.* By Lemma 1.1, $\mathbb E\,m(X_x-Y_y)=\frac{d-1}2+\sum_nc_n\,\mathbb E\,w^{n(X_x-Y_y)}$ and
$\mathbb E\,w^{n(X_x-Y_y)}=\sum_{a,b}w^{n(a-b)}p(a,b|x,y)=\langle\Phi_D|A_x^n\otimes B_y^{-n}|\Phi_D\rangle
=\mathrm{tr}(A_x^n(B_y^T)^{-n})$ (Lemma 1.3). Hence $\mathbb E\,w^{n(X_1-Y_1)}=\mathrm{tr}(R_1^nR_2^{-n})$,
$\mathbb E\,w^{n(Y_1-X_2)}=\mathrm{tr}(R_3^{-n}R_2^{n})=\mathrm{tr}(R_2^nR_3^{-n})$,
$\mathbb E\,w^{n(X_2-Y_2)}=\mathrm{tr}(R_3^nR_4^{-n})$,
$\mathbb E\,w^{n(Y_2-X_1-1)}=w^{-n}\mathrm{tr}(R_1^{-n}R_4^n)=\mathrm{tr}(L_n^{(3)})$. $\square$

**Lemma 1.6 (faithfulness).** $\mathrm{tr}(X^*X)=\frac1D\sum_{k,l}|X_{kl}|^2\ge0$, with equality iff $X=0$.
(Physically: both reduced states of $|\Phi_D\rangle$ equal $1/D$, which has full rank; this is the only place
where maximal entanglement *with full-rank marginals* enters.)

**Optimality.** $\lambda_d:=\inf\{S(R):D\ge1,R\in\mathcal U_D\}$; $R$ is *optimal* if $S(R)=\lambda_d$.
$I_{\rm ME}(d):=4-2\lambda_d/(d-1)$ whenever (OPT_d) holds (then $\lambda_d=S(R^{\rm DKZ})$ and
$I_{\rm ME}(d)=\frac{4}{d(d-1)}\sum_{j=1}^{d-1}(d-j)\sec\frac{\pi j}{2d}$, lab LOG section 8).

### 1.3 Algebraic framework

$\Gamma_d:=\mathbb Z_d*\mathbb Z_d*\mathbb Z_d*\mathbb Z_d$ with generators $r_1,\dots,r_4$ ($r_k^d=e$);
$\mathcal A_d:=\mathbb C[\Gamma_d]$ with $r^*=r^{-1}$. Each $R\in\mathcal U_D$ gives the unital $*$-representation
$\pi_R:\mathcal A_d\to M_D(\mathbb C)$, $r_k\mapsto R_k$, and conversely. Put
$\ell^{(0)}_n=r_1^nr_2^{-n}$, $\ell^{(1)}_n=r_2^nr_3^{-n}$, $\ell^{(2)}_n=r_3^nr_4^{-n}$,
$\ell^{(3)}_n=w^{-n}r_4^nr_1^{-n}$ and
$S_{\rm op}:=2(d-1)\cdot1+\sum_{n=1}^{d-1}c_n\sum_{k=0}^{3}\ell^{(k)}_n$, so $S(R)=\mathrm{tr}\,\pi_R(S_{\rm op})$.
Target polynomials (all in $\mathcal A_d$; $n=1,\dots,d-1$, $t=1,2,3$):
$$M_{n,t}:=\sum_{k=0}^3 i^{-tk}\ell^{(k)}_n,\qquad
T_n:=\sum_{k=0}^3\Big[(\ell^{(k)}_n)^{-1}-(1-i)z^n-iz^{2n}\ell^{(k)}_n\Big],\qquad
V_n:=(\ell^{(0)}_n)^{-1}-(1-i)z^n-iz^{2n}\ell^{(0)}_n .$$

### 1.4 The DKZ strategy and the canonical model

DKZ bases (lab LOG section 3): Alice $|a\rangle_x=d^{-1/2}\sum_{k=0}^{d-1}w^{k(a+\alpha_x)}|k\rangle$,
$\alpha=(\tfrac12,0)$; Bob $|b\rangle_y=d^{-1/2}\sum_kw^{-k(b+\beta_y)}|k\rangle$, $\beta=(\tfrac14,-\tfrac14)$.
Let $G:=\mathrm{diag}(z^k)_{k=0}^{d-1}$, $X|k\rangle=|k+1\bmod d\rangle$, $Z:=\mathrm{diag}(w^k)=G^4$,
$F|a\rangle=d^{-1/2}\sum_kw^{ka}|k\rangle$ (DFT). Then $XF=FZ^{-1}$, i.e. $FZF^*=X^{-1}$;
$|a\rangle_1=G^2F|a\rangle$, $|a\rangle_2=F|a\rangle$, $\overline{|b\rangle_y}=G^{4\beta_y}F|b\rangle$ and
$(|v\rangle\langle v|)^T=|\bar v\rangle\langle\bar v|$. Hence
$$R^{\rm DKZ}_k=G^{3-k}X^{-1}G^{k-3}\qquad(k=1,2,3,4).$$
(Each measurement basis is the Fourier basis $F|a\rangle$ or $\bar F|b\rangle$ dressed by a diagonal unitary:
this is the DKZ form.)

**Canonical model.** $R^c_1:=X$, $\Pi^c:=|0\rangle\langle0|$, $D^c:=1-(1+i)\Pi^c=\mathrm{diag}(-i,1,\dots,1)$,
$R^c_{k+1}:=zD^cR^c_k$ ($k=1,2,3$).

**Lemma 1.7 (properties of DKZ).**
(a) Each $R^{\rm DKZ}_k$ is a weighted cyclic shift with $(R^{\rm DKZ}_k)^d=1$.
(b) $G^{-1}R^{\rm DKZ}_kG=R^{\rm DKZ}_{k+1}$ ($k=1,2,3$) and $G^{-1}R^{\rm DKZ}_4G=wR^{\rm DKZ}_1$.
(c) For every $n$ all four links of $R^{\rm DKZ}$ coincide.
(d) $R^{\rm DKZ}_2R^{\rm DKZ\,-1}_1=z\,(1-(1+i)|d-1\rangle\langle d-1|)$.
(e) $W_0|j\rangle:=z^{-2j}|d-1-j\rangle$ ($j=0,\dots,d-1$) is unitary and $W_0^*R^{\rm DKZ}_kW_0=R^c_k$ ($k=1..4$).
(f) $R^c$ (hence $R^{\rm DKZ}$) satisfies (E_n) and (T_n) for every $n$, and $R^c_k$ generate $M_d(\mathbb C)$.

*Proof.* (a) $R^{\rm DKZ}_k|j\rangle=z^{-a j+a((j-1)\bmod d)}|j-1\rangle$ with $a=3-k$: a monomial matrix whose
$d$-th power is $z^{a\sum_j[(j-1)\bmod d-j]}=z^{a(d(d-1)/2-d(d-1)/2)}=1$ times the identity.
(b) The first claims are immediate; the last is equivalent to $G^{-4}X^{-1}G^{4}=wX^{-1}$, i.e.
$Z^{-1}XZ=w^{-1}X$, which follows from $ZXZ^{-1}=wX$.
(c) By (b), $\varphi(x):=G^{-1}xG$ maps $R_1\mapsto R_2\mapsto R_3\mapsto R_4\mapsto wR_1$, hence
$\varphi(L^{(k)}_n)=L^{(k+1)}_n$ ($k=0,1,2$). $L^{(0)}_n=R_1^nG^{-1}R_1^{-n}G$ is diagonal (a monomial matrix
conjugates diagonals to diagonals), so it commutes with $G$, and $L^{(k)}_n=G^{-k}L^{(0)}_nG^{k}=L^{(0)}_n$.
(d) $R_2R_1^{-1}=GX^{-1}G^{-1}G^2XG^{-2}=G(X^{-1}GX)G^{-2}$ and $X^{-1}GX=\mathrm{diag}(z^{(j+1)\bmod d})
=zG\,\mathrm{diag}(1,\dots,1,z^{-d})$ with $z^{-d}=-i=1-(1+i)$.
(e) $W_0$ is monomial with unimodular entries. From (a), $R^{\rm DKZ}_1|j\rangle=z^{-2}|j-1\rangle$ ($j\ge1$),
$R^{\rm DKZ}_1|0\rangle=z^{2(d-1)}|d-1\rangle$; hence $W_0^*R_1^{\rm DKZ}W_0=X$, and $W_0^*|d-1\rangle\langle
d-1|W_0=|0\rangle\langle0|$. By (c),(d): $R^{\rm DKZ}_{k+1}=zD^{\rm DKZ}R^{\rm DKZ}_k$ with
$D^{\rm DKZ}=1-(1+i)|d-1\rangle\langle d-1|$, and conjugating by $W_0$ gives the recursion of $R^c$.
(Also verified exactly by `check_support.py`, check (S5).)
(f) Unitary conjugation preserves polynomial identities, so (E_n) follows from (c). For (T_n) use Lemma 3.3
below in the canonical frame: $P_j=X^j|0\rangle\langle0|X^{-j}=|j\rangle\langle j|$ are orthogonal, so
$z^nU_n=\mathcal Q(P_{n-1})\cdots\mathcal Q(P_0)=\mathcal Q(\sum_{j<n}|j\rangle\langle j|)$, which is two-valued
(Lemma 3.1). Generation: the algebra contains $\Pi^c=(zU_1-1)/(i-1)$ (Lemma 3.1), all
$X^j\Pi^cX^{-j}=|j\rangle\langle j|$ and $X|j\rangle\langle j|=|j+1\rangle\langle j|$, i.e. all matrix units. $\square$

Every weighted cyclic shift $Y=\sum_jc_j|j+1\rangle\langle j|$ with $Y^d=1$ ($\prod_jc_j=1$) equals $\Delta X\Delta^{-1}$
for a diagonal unitary $\Delta$ ($\delta_0=1,\ \delta_{j+1}=c_j\delta_j$), so its eigenbasis is $\{\Delta F|a\rangle\}$.
Thus "unitarily equivalent to $R^c$ by a monomial unitary" = "DKZ form" in the sense of OQP 27B.

### 1.5 Symmetries

Let $\rho,\sigma$ be the unital $*$-automorphisms of $\mathcal A_d$ with $\rho:r_1\mapsto r_2\mapsto r_3\mapsto
r_4\mapsto wr_1$ and $\sigma:r_k\mapsto r_{5-k}^{-1}$, and $\tau$ the linear $*$-anti-automorphism with
$\tau(r_k)=r_k$ (word reversal). ($\rho$ is well defined since $(wr_1)^d=e$, and invertible.) Let $G$ be the
group they generate.

**Lemma 1.8.** (i) Every $\alpha\in G$ is an automorphism or anti-automorphism with $\alpha(r_k)=\lambda_kr_{j_k}^{s_k}$,
$\lambda_k\in\{w^m\}$, $s_k=\pm1$; hence $G$ is finite (only finiteness is used; $|G|=16d$ for
$d=3,\dots,9$, computed in `check_support.py`).
(ii) For $R\in\mathcal U_D$ and $\alpha\in G$ put $R^\alpha_k:=\pi_R(\alpha(r_k))$ if $\alpha$ is an automorphism and
$R^\alpha_k:=\pi_R(\alpha(r_k))^T$ otherwise. Then $R^\alpha\in\mathcal U_D$ and
$\mathrm{tr}\,\pi_{R^\alpha}(x)=\mathrm{tr}\,\pi_R(\alpha(x))$ for all $x\in\mathcal A_d$; more precisely
$\pi_{R^\alpha}=\pi_R\circ\alpha$ resp. $\pi_{R^\alpha}(x)=\pi_R(\alpha(x))^T$.
(iii) $\alpha(S_{\rm op})-S_{\rm op}\in[\mathcal A_d,\mathcal A_d]$ (span of commutators) for all $\alpha\in G$.

*Proof.* (i) True for the generators, and the set of such maps is closed under composition and finite.
(ii) $\lambda R_j^{\pm1}$ is unitary with $d$-th power $1$; $\pi_R\circ\alpha$ (resp. $x\mapsto\pi_R(\alpha(x))^T$,
a composition of two anti-homomorphisms and a homomorphism) is a unital homomorphism agreeing with
$\pi_{R^\alpha}$ on generators; $\mathrm{tr}(Y^T)=\mathrm{tr}\,Y$.
(iii) It suffices to check generators, since automorphisms and anti-automorphisms preserve
$[\mathcal A_d,\mathcal A_d]$. $\rho$ permutes the links exactly: $\ell^{(k)}_n\mapsto\ell^{(k+1)}_n$ and
$\ell^{(3)}_n\mapsto w^{-n}(wr_1)^nr_2^{-n}=\ell^{(0)}_n$. $\sigma$: $\ell^{(0)}_n\mapsto r_4^{-n}r_3^n\sim\ell^{(2)}_n$,
$\ell^{(1)}_n\mapsto r_3^{-n}r_2^n\sim\ell^{(1)}_n$, $\ell^{(2)}_n\mapsto r_2^{-n}r_1^n\sim\ell^{(0)}_n$,
$\ell^{(3)}_n\mapsto w^{-n}r_1^{-n}r_4^n\sim\ell^{(3)}_n$, where $u\sim v$ means $u-v\in[\mathcal A_d,\mathcal A_d]$
($ab\sim ba$). $\tau$ maps each $\ell^{(k)}_n$ to the reversed word, which is a cyclic rotation of it. $\square$

**Lemma 1.9 ($G$ preserves the DKZ class).** Call $R\in\mathcal U_d$ *DKZ-structured* if there is a projection $\Pi$
with (a) $R_{k+1}=zDR_k$ ($k=1,2,3$), $D=1-(1+i)\Pi$, and (b) $P_j=R_1^j\Pi R_1^{-j}$, $0\le j<d$, mutually
orthogonal with $\sum_jP_j=1$. Then: $R$ is DKZ-structured iff $R$ is unitarily equivalent to $R^c$ (iff to
$R^{\rm DKZ}$); and for every $\alpha\in G$ the strategy $(R^{\rm DKZ})^\alpha$ of Lemma 1.8(ii) is unitarily equivalent
to $R^{\rm DKZ}$. In particular every $G$-image of $R^{\rm DKZ}$ satisfies (E_n) and (T_n) for all $n$.

*Proof.* "If": $R^c$ is DKZ-structured with $\Pi^c$, and the conditions are unitarily invariant. "Only if":
Lemma 3.7 and Theorem 3.9 below ($M=1$ since $D=d$). For the second claim it suffices, by Lemma 1.8(ii)
and the invariance of unitary equivalence under $R\mapsto R^\alpha$ ($(WRW^*)^\alpha=WR^\alpha W^*$ for automorphisms,
$=\bar WR^\alpha W^T$ for $\tau$), to show that each generator maps DKZ-structured tuples to DKZ-structured
tuples. Let $R$ be DKZ-structured; recall $wR_1=zDR_4$ (Remark 3.10) and that the $P_j$ commute with each other.
$\rho$: $(R_2,R_3,R_4,wR_1)$ satisfies (a) with the same $\Pi$; for (b), $R_2^j=z^jD_0\cdots D_{j-1}R_1^j$
(Lemma 3.3) with $D_l=1-(1+i)P_l$ commuting with $P_j$, so $R_2^j\Pi R_2^{-j}=P_j$. Iterating (the family
$\{P_j\}$ is the same for $R_2$), $R_k^j\Pi R_k^{-j}=P_j$ for $k=1,2,3,4$ and all $j\ge0$; since $R_k^d=1$ this holds
for all $j\in\mathbb Z$ with $P_{j}:=P_{j\bmod d}$. $\sigma$: $R'=(R_4^{-1},R_3^{-1},R_2^{-1},R_1^{-1})$.
From $R_{k+1}=zDR_k$: $R_k^{-1}=R_{k+1}^{-1}zD$, so $R'_{k+1}=zR'_kD=zD'_kR'_k$ with $D'_k:=R'_kDR'^{-1}_k$;
$D'_k$ does not depend on $k$ because $R'^{-1}_{1}R'_{k}=(zD)^{k-1}$ commutes with $D$; so (a) holds
with $\Pi'=R_4^{-1}\Pi R_4$, and $R'^j_1\Pi'R'^{-j}_1=R_4^{-(j+1)}\Pi R_4^{j+1}=P_{-(j+1)}$ by the previous
sentence: (b) holds. $\tau$: $R'=(R_1^T,\dots,R_4^T)$;
$R'_{k+1}=zR'_kD^T=zD'R'_k$ with $D'=R'_kD^TR'^{-1}_k$, independent of $k$ since $R'^{-1}_kR'_{k'}$ is a power of
$(zD)^T$; $\Pi'=R_1^T\Pi^TR_1^{-T}$ and $R'^j_1\Pi'R'^{-j}_1=(P_{-(j+1)})^T$: (b) holds. The last sentence follows
from Lemma 1.7(f). (Numerically confirmed for all $8d$ distinct images, $d\le30$.) $\square$

---------------------------------------------------------------------------------------------------------

## 2. From a certificate to equality conditions (faithful-trace step)

**Definition 2.1 (symmetric SOS certificate).** A certificate for $d$ consists of $\lambda\in\mathbb R$ and finitely
many blocks $b$, each with elements $e_{b,1},\dots,e_{b,n_b}\in\mathcal A_d$ and a Hermitian positive semidefinite
$X_b\in M_{n_b}(\mathbb C)$, such that

(Id) $L(S_{\rm op}-\lambda1)=\sum_b\sum_{i,j}(X_b)_{ij}\,L(e_{b,j}^*e_{b,i})$ for every linear functional
$L:\mathcal A_d\to\mathbb C$ that is tracial ($L(xy)=L(yx)$) and $G$-invariant ($L\circ\alpha=L$, $\alpha\in G$).

Its *support* is $\mathcal S:=\sum_b\mathcal S_b$, $\mathcal S_b:=\{\sum_iu_ie_{b,i}:u\in\mathrm{ran}\,X_b\}$, and
$\mathcal I$ denotes the two-sided ideal of $\mathcal A_d$ generated by $\bigcup_{\alpha\in G}\alpha(\mathcal S)$.

(The lab certificates have $X_b=N_bY_bN_b^*$ with $N_b$ of full column rank and $Y_b$ positive definite, so
$\mathrm{ran}\,X_b=\mathrm{ran}\,N_b$. (Id) is verified by reducing every word to a canonical representative of
its class under cyclic rotation and $G$ (with phases), a word being dropped when its class forces $L=0$;
equality of the reduced coefficient vectors implies (Id) for every tracial $G$-invariant $L$.)

**Lemma 2.2 (orbit direct sum).** For $R\in\mathcal U_D$ let $R^\oplus:=\bigoplus_{\alpha\in G}R^\alpha\in
\mathcal U_{D|G|}$. Then $L_R:=\mathrm{tr}\circ\pi_{R^\oplus}=|G|^{-1}\sum_{\alpha\in G}\mathrm{tr}\circ\pi_R\circ\alpha$
is a tracial, $G$-invariant state on $\mathcal A_d$, and $S(R^\oplus)=S(R^\alpha)=S(R)$ for all $\alpha$.

*Proof.* The normalised trace of a block-diagonal matrix with $|G|$ blocks of equal size is the average of the
block traces; use Lemma 1.8(ii). Traciality is that of $\mathrm{tr}$. For $\beta\in G$,
$L_R\circ\beta=|G|^{-1}\sum_\alpha\mathrm{tr}\circ\pi_R\circ(\alpha\circ\beta)=L_R$ because $\alpha\mapsto\alpha\circ\beta$
is a bijection of $G$. Finally $S(R^\alpha)=\mathrm{tr}\,\pi_R(\alpha(S_{\rm op}))=\mathrm{tr}\,\pi_R(S_{\rm op})=S(R)$
by Lemma 1.8(iii), since traces vanish on commutators. $\square$

**Lemma 2.3.** If $X=\sum_ku_ku_k^*$ ($u_k\in\mathbb C^n$), then $\sum_{i,j}X_{ij}e_j^*e_i=\sum_kf_k^*f_k$ with
$f_k=\sum_i(u_k)_ie_i$. For $X\succeq0$ take the spectral decomposition; then $\{u_k\}$ spans $\mathrm{ran}\,X$.

*Proof.* $X_{ij}=\sum_k(u_k)_i\overline{(u_k)_j}$, and $\big(\sum_j(u_k)_je_j\big)^*=\sum_j\overline{(u_k)_j}e_j^*$. $\square$

**Proposition 2.4 (bound and faithful-trace step).** Given a certificate:
(a) $S(R)\ge\lambda$ for every $D$ and $R\in\mathcal U_D$;
(b) if $S(R)=\lambda$, then $\pi_R(g)=0$ for every $g\in\mathcal I$.

*Proof.* Let $L=L_R$ (Lemma 2.2). By (Id) and Lemma 2.3,
$$S(R)-\lambda=S(R^\oplus)-\lambda=\sum_{b,k}L(f_{b,k}^*f_{b,k})=\sum_{b,k}\mathrm{tr}\big(F_{b,k}^*F_{b,k}\big),\qquad
F_{b,k}:=\pi_{R^\oplus}(f_{b,k}),$$
a sum of non-negative terms: (a). If $S(R)=\lambda$, all terms vanish, so $F_{b,k}=0$ (Lemma 1.6, which is
exactly where faithfulness of the trace, i.e. full-rank marginals of $|\Phi_D\rangle$, is used). As
$\pi_{R^\oplus}=\bigoplus_\alpha\pi_{R^\alpha}$, we get $\pi_{R^\alpha}(f_{b,k})=0$ for all $\alpha$, hence
$\pi_R(\alpha(f_{b,k}))=0$ by Lemma 1.8(ii). The $f_{b,k}$ span $\mathcal S_b$, so $\pi_R(\alpha(f))=0$ for all
$f\in\mathcal S$, $\alpha\in G$; and $\pi_R(xgy)=\pi_R(x)\pi_R(g)\pi_R(y)$. $\square$

Note: the identity (Id) holds only modulo the symmetry relations, so it cannot be evaluated on an arbitrary
strategy directly; the orbit direct sum $R^\oplus$ is what makes the evaluation legitimate. This justifies the
phrase "faithful trace $\Rightarrow f(R)=0$ for every $f$ in the support" of the original sketch.

**Proposition 2.5 (links).** Let $R\in\mathcal U_D$ and $1\le n\le d-1$.
(a) If $\pi_R(M_{n,t})=0$ for $t=1,2,3$, then (E_n): $L^{(0)}_n=L^{(1)}_n=L^{(2)}_n=L^{(3)}_n$.
(b) If $\pi_R(V_n)=0$, or if (E_n) holds and $\pi_R(T_n)=0$, then (T_n) holds:
$(U_n-z^{-n})(U_n-iz^{-n})=0$ with $U_n=R_1^nR_2^{-n}$.

*Proof.* (a) Put $\ell_k:=\pi_R(\ell^{(k)}_n)$ and $\hat\ell_t:=\sum_ki^{-tk}\ell_k$. Since $\sum_{t=0}^3i^{t(k-k')}=4\delta_{kk'}$,
$\ell_k=\frac14\sum_{t=0}^3i^{tk}\hat\ell_t=\frac14\hat\ell_0$ for every $k$. (b) In $\mathcal A_d$ (with
$u=\ell^{(0)}_n$): $V_n=-iz^{2n}u^{-1}(u-z^{-n})(u-iz^{-n})$, by expanding
$-iz^{2n}u^{-1}(u^2-(1+i)z^{-n}u+iz^{-2n})=-iz^{2n}u-(1-i)z^n+u^{-1}$. As $-iz^{2n}U_n^{-1}$ is invertible,
$\pi_R(V_n)=0\iff$ (T_n). Under (E_n), $\pi_R((\ell^{(k)}_n)^{-1})=\ell_k^{-1}=U_n^{-1}$ and $\ell_k=U_n$, so
$\pi_R(T_n)=4\,\pi_R(V_n)$. $\square$

**Theorem 2.6 (certificate $\Rightarrow$ (OPT_d) and (RIG_d)).** Suppose a certificate exists with
$\lambda=S(R^{\rm DKZ})$ and one of
(C1) $M_{n,t}\in\mathcal S$ and $T_n\in\mathcal S$ for all $1\le n\le m_d$, $t=1,2,3$; or
(C2) $M_{1,t}\in\mathcal I$ ($t=1,2,3$) and $V_n\in\mathcal I$ for $1\le n\le m_d$.
Then (OPT_d) and (RIG_d) hold.

*Proof.* Proposition 2.4(a) gives $S(R)\ge S(R^{\rm DKZ})$: (OPT_d); so optimal means $S(R)=\lambda$, and
Proposition 2.4(b) applies. Under (C1), Proposition 2.5(a) gives (E_n) and then 2.5(b) gives (T_n) for
$n\le m_d$; (E_1) contains (E_1'). Under (C2), 2.5(a) with $n=1$ and 2.5(b) directly. $\square$

---------------------------------------------------------------------------------------------------------

## 3. Algebraic rigidity

Throughout this section $\mathcal A$ is a unital $*$-algebra over $\mathbb C$ with

(C*) $x^*x=0\Rightarrow x=0$.

This holds in $M_D(\mathbb C)$ ($\mathrm{Tr}(x^*x)=\sum|x_{ij}|^2$), in every C*-algebra ($\|x^*x\|=\|x\|^2$)
and hence in every von Neumann algebra. A projection is $P=P^*=P^2$. Write $a:=i-1$ and, for a projection
$P$, $\mathcal Q(P):=1+(i-1)P$ ("quarter-phase unitary": eigenvalue $1$ on $\ker P$, $i$ on $\mathrm{ran}\,P$).

**Facts.** $\mathcal Q(P)^*=1-(1+i)P=\mathcal Q(P)^{-1}$ (since $(1+aP)(1+\bar aP)=1+(a+\bar a+|a|^2)P=1$ as
$a+\bar a=-2$, $|a|^2=2$); for $P\perp Q$ (i.e. $PQ=0$, hence $QP=(PQ)^*=0$), $P+Q$ is a projection and
$\mathcal Q(P)\mathcal Q(Q)=\mathcal Q(P+Q)$; $R\mathcal Q(P)R^{-1}=\mathcal Q(RPR^{-1})$ for unitary $R$.

**Lemma 3.1 (two-valued unitaries).** For a unitary $A\in\mathcal A$ the following are equivalent:
(i) $(A-1)(A-i)=0$; (ii) $A+iA^*=(1+i)1$; (iii) $A=\mathcal Q(E)$ for a projection $E$.
In that case $E=(A-1)/(i-1)$ (unique).

*Proof.* (i)$\Leftrightarrow$(ii): multiply $A^2-(1+i)A+i=0$ on the left by $A^*=A^{-1}$, resp. (ii) by $A$.
(i)$\Rightarrow$(iii): let $E:=(A-1)/(i-1)$. Since $(A-1)^2=(A-1)(A-i)+(i-1)(A-1)=(i-1)(A-1)$, $E^2=E$. By
(ii), $A^*=(1-i)+iA$, so $E^*=(A^*-1)/(-1-i)=i(A-1)/(-1-i)=(A-1)\cdot\frac{-1-i}{2}=E$, because
$\frac{1}{i-1}=\frac{-1-i}{2}$. Clearly $A=1+(i-1)E$. (iii)$\Rightarrow$(i):
$(\mathcal Q(E)-1)(\mathcal Q(E)-i)=aE\,(1-i+aE)=aE(1-i+a)=0$. $\square$

**Lemma 3.2 (Jordan-type lemma).** Let $P,Q\in\mathcal A$ be projections and
$A:=\mathcal Q(P)\mathcal Q(Q)=(1+(i-1)P)(1+(i-1)Q)$. Then $(A-1)(A-i)=0$ if and only if $PQ=0$.

*Proof.* ($\Leftarrow$) If $PQ=0$ then $A=\mathcal Q(P+Q)$ with $P+Q$ a projection; Lemma 3.1.
($\Rightarrow$) $A$ is unitary. With $a=i-1$, $\bar a=-1-i$, $a^2=-2i$, $\bar a^2=2i$, $i\bar a=1-i$:
$$A=1+a(P+Q)+a^2PQ,\qquad iA^*=i+(1-i)(P+Q)-2\,QP,$$
$$A+iA^*=(1+i)+\big(a+1-i\big)(P+Q)-2i\,PQ-2\,QP=(1+i)-2\,(iPQ+QP).$$
By Lemma 3.1(ii), $iPQ+QP=0$. Multiplying by $P$ on both sides: $(1+i)PQP=0$, so $PQP=0$. Then
$(QP)^*(QP)=PQ^2P=PQP=0$, hence $QP=0$ by (C*), and $PQ=(QP)^*=0$. $\square$

*Remark 3.2' (the spectral picture, including all non-generic configurations).* The proof above treats all
configurations at once and needs no decomposition. For orientation, in finite dimension Jordan's lemma
(C. Jordan 1875; P. R. Halmos, "Two subspaces", Trans. AMS 144 (1969)) splits $\mathbb C^D$ into $P,Q$-invariant
pieces of dimension 1 and 2. On a line, $(P,Q)\in\{0,1\}^2$ and $A=1,i,i,-1$ respectively: the eigenvalue $-1$
occurs exactly on $\mathrm{ran}P\cap\mathrm{ran}Q$ (the non-generic case $P=Q$ gives $A=1-2P$). On a
2-dimensional piece in generic position ($P=|p\rangle\langle p|$, $Q=|q\rangle\langle q|$, $0<c:=|\langle p|q\rangle|<1$):
$\det A=i^2=-1$ and $\mathrm{tr}A=2+2a+a^2c^2=2i(1-c^2)$, so the eigenvalues are $e^{i\phi},e^{i(\pi-\phi)}$ with
$\sin\phi=1-c^2\in(0,1)$, never in $\{1,i\}$. Hence $\sigma(A)\subseteq\{1,i\}$ iff there is no common vector and no
generic piece, i.e. iff $\mathrm{ran}P\perp\mathrm{ran}Q$ - consistent with Lemma 3.2.

**Lemma 3.3 (product formula).** Let $R_1,R_2\in\mathcal A$ be unitaries and $\Pi$ a projection with
$R_2=zDR_1$, $D:=1-(1+i)\Pi=\mathcal Q(\Pi)^{-1}$. Put $P_j:=R_1^j\Pi R_1^{-j}$ ($j\in\mathbb Z$). Then for $n\ge1$
$$z^nR_1^nR_2^{-n}=\mathcal Q(P_{n-1})\,\mathcal Q(P_{n-2})\cdots\mathcal Q(P_0).$$

*Proof.* With $D_j:=R_1^jDR_1^{-j}=\mathcal Q(P_j)^{-1}$, induction gives $(DR_1)^n=D_0D_1\cdots D_{n-1}R_1^n$:
$(DR_1)^{n+1}=D_0\cdots D_{n-1}R_1^nDR_1=D_0\cdots D_{n-1}(R_1^nDR_1^{-n})R_1^{n+1}$. Hence
$R_2^n=z^nD_0\cdots D_{n-1}R_1^n$ and $R_1^nR_2^{-n}=z^{-n}D_{n-1}^{-1}\cdots D_0^{-1}$. $\square$

**Lemma 3.4 (orthogonality induction).** In the situation of Lemma 3.3 let $K\ge1$ and assume
$(z^nU_n-1)(z^nU_n-i)=0$ for $U_n:=R_1^nR_2^{-n}$ and $n=2,\dots,K+1$. Then $P_0,P_1,\dots,P_K$ are mutually
orthogonal.

*Proof.* Induction on $k\le K$; $H(k)$: $P_0,\dots,P_k$ mutually orthogonal ($H(0)$ is empty). Assume $H(k-1)$.
Then $E_k:=P_0+\dots+P_{k-1}$ is a projection and $\mathcal Q(P_{k-1})\cdots\mathcal Q(P_0)=\mathcal Q(E_k)$
(repeatedly $\mathcal Q(P_j)\mathcal Q(E_j)=\mathcal Q(E_{j+1})$ since $P_j\perp E_j$). By Lemma 3.3,
$z^{k+1}U_{k+1}=\mathcal Q(P_k)\mathcal Q(E_k)$, and by hypothesis ($n=k+1$) and Lemma 3.2, $P_kE_k=0$. For
$j<k$: $P_j=E_kP_j$, so $P_kP_j=P_kE_kP_j=0$. $\square$

*What is used here:* two-valuedness for $n=2,\dots,K+1$ only, and each application of Lemma 3.2 is to a
product of exactly two quarter-phase unitaries, the second one being $\mathcal Q(E_k)$ with $E_k$ the projection
accumulated so far. The products have the stated form $\mathcal Q(P_{n-1})\cdots\mathcal Q(P_0)$ (Lemma 3.3), leftmost
factor $\mathcal Q(P_{n-1})$.

**Lemma 3.5 (cyclic transport: $K=\lfloor d/2\rfloor$ suffices).** If moreover $R_1^d=1$ and $K\ge\lfloor d/2\rfloor$,
then $P_{j+d}=P_j$ for all $j$, and $P_jP_k=0$ whenever $j\not\equiv k\pmod d$.

*Proof.* $P_{j+d}=R_1^jR_1^d\Pi R_1^{-d}R_1^{-j}=P_j$, and $R_1^jP_lR_1^{-j}=P_{j+l}$. Let $j\not\equiv k$ and
$r:=(k-j)\bmod d\in\{1,\dots,d-1\}$. If $r\le K$: $P_jP_k=R_1^j(P_0P_r)R_1^{-j}=0$ by Lemma 3.4. If $r>K$, then
$1\le d-r\le d-K-1\le\lceil d/2\rceil-1\le\lfloor d/2\rfloor\le K$, and $j\equiv k+(d-r)$, so
$P_kP_j=R_1^k(P_0P_{d-r})R_1^{-k}=0$ and $P_jP_k=(P_kP_j)^*=0$. $\square$

**Lemma 3.6 (resolution of the identity).** If, in addition, $R_2^d=1$, then $\sum_{j=0}^{d-1}P_j=1$.

*Proof.* $U_d=R_1^dR_2^{-d}=1$ and $z^d=i$. By Lemma 3.3 and Lemma 3.5,
$i\cdot1=z^dU_d=\mathcal Q(P_{d-1})\cdots\mathcal Q(P_0)=\mathcal Q(E)$ with $E:=\sum_{j<d}P_j$ a projection. So
$(i-1)E=(i-1)1$, i.e. $E=1$. $\square$

(The original sketch obtained the last projection from $U_d=1$ after using all $T_n$, $n\le d-1$; Lemma 3.5
shows that the half-range $n\le\lfloor d/2\rfloor+1$ already gives all orthogonality relations.)

**Lemma 3.7 (imprimitivity $\Rightarrow$ matrix units).** Let $R\in\mathcal A$ be unitary with $R^d=1$ and $\Pi$ a
projection such that $P_j:=R^j\Pi R^{-j}$ ($0\le j\le d-1$) are mutually orthogonal and $\sum_jP_j=1$. Put
$e_{jk}:=R^j\Pi R^{-k}$ ($0\le j,k\le d-1$). Then
$$e_{jk}e_{lm}=\delta_{kl}e_{jm},\qquad e_{jk}^*=e_{kj},\qquad \textstyle\sum_je_{jj}=1,\qquad
e_{00}=\Pi,\qquad R=\sum_{j=0}^{d-1}e_{j+1\bmod d,\;j}.$$
Consequently $\pi(|j\rangle\langle k|):=e_{jk}$ extends linearly to a unital $*$-homomorphism
$\pi:M_d(\mathbb C)\to\mathcal A$ with $\pi(X)=R$, $\pi(|0\rangle\langle0|)=\Pi$.

*Proof.* $e_{jk}^*=R^k\Pi R^{-j}=e_{kj}$. $e_{jk}e_{lm}=R^j(\Pi R^{l-k}\Pi)R^{-m}$ and
$\Pi R^{l-k}\Pi=(P_0P_{l-k})R^{l-k}$, which is $\Pi$ if $l=k$ and $0$ otherwise ($l-k\not\equiv0$). $e_{jj}=P_j$.
$R=R\sum_je_{jj}=\sum_jR^{j+1}\Pi R^{-j}$, and $R^d\Pi R^{-(d-1)}=e_{0,d-1}$ as $R^d=1$. Linearity and the
matrix-unit relations give multiplicativity, $*$-compatibility and unitality of $\pi$. $\square$

**Theorem 3.8 (algebraic rigidity).** Let $d\ge3$ and let $R_1,\dots,R_4\in\mathcal A$ be unitaries with
$R_1^d=R_2^d=1$ such that

* (E_1') $R_1R_2^{-1}=R_2R_3^{-1}=R_3R_4^{-1}$;
* (T_n) $(U_n-z^{-n})(U_n-iz^{-n})=0$ for $U_n=R_1^nR_2^{-n}$, $n=1,\dots,m_d=\lfloor d/2\rfloor+1$.

Then there is a unital $*$-homomorphism $\pi:M_d(\mathbb C)\to\mathcal A$ with $R_k=\pi(R^c_k)$, $k=1,2,3,4$.
In particular $R_3^d=R_4^d=1$ and $w^{-1}R_4R_1^{-1}=R_1R_2^{-1}$ hold automatically, and $R$ satisfies (E_n)
and (T_n) for every $n$.

*Proof.* (T_1) says $(zU_1-1)(zU_1-i)=0$; by Lemma 3.1, $zU_1=\mathcal Q(\Pi)$ for a projection $\Pi$, i.e.
$R_2=zDR_1$ with $D=\mathcal Q(\Pi)^{-1}=1-(1+i)\Pi$. By (E_1'), $R_3=U_1^{-1}R_2=zDR_2$ and $R_4=zDR_3$.
(T_n) for $2\le n\le m_d$ is the hypothesis of Lemma 3.4 with $K=\lfloor d/2\rfloor$ (note $m_d\le d-1$ for $d\ge3$);
Lemmas 3.5, 3.6 show that $P_j=R_1^j\Pi R_1^{-j}$, $0\le j<d$, are mutually orthogonal with sum $1$. Lemma 3.7
gives $\pi$ with $\pi(X)=R_1$, $\pi(|0\rangle\langle0|)=\Pi$, hence $\pi(D^c)=D$ and inductively
$\pi(R^c_{k+1})=z\pi(D^c)\pi(R^c_k)=zDR_k=R_{k+1}$. The last sentence follows from Lemma 1.7(f), since unital
$*$-homomorphisms preserve polynomial identities. $\square$

**Theorem 3.9 (finite-dimensional form).** If $\mathcal A=M_D(\mathbb C)$ in Theorem 3.8, then $d\mid D$ and there is
a unitary $W:\mathbb C^d\otimes\mathbb C^M\to\mathbb C^D$, $M=D/d$, with $R_k=W(R^c_k\otimes1_M)W^*$ for $k=1,\dots,4$
and $\Pi=W(|0\rangle\langle0|\otimes1_M)W^*$.

*Proof.* Let $K:=\mathrm{ran}\,e_{00}$, $M:=\dim K$, and identify $K\cong\mathbb C^M$ by an orthonormal basis. Define
$W(|j\rangle\otimes v):=e_{j0}v$. Isometry: $\langle e_{j0}v,e_{k0}u\rangle=\langle v,e_{0j}e_{k0}u\rangle=\delta_{jk}\langle
v,u\rangle$. Surjectivity: $x=\sum_je_{jj}x=\sum_je_{j0}(e_{0j}x)$ with $e_{0j}x\in K$. So $W$ is unitary and
$D=dM$. On $W(|l\rangle\otimes v)$ both $W(|j\rangle\langle k|\otimes1)W^*$ and $e_{jk}$ act as
$\delta_{kl}\,e_{j0}v$, so $\pi(x)=W(x\otimes1_M)W^*$. $\square$

*Remark 3.10 (the four equalities in (E_1)).* Only three of the four equalities of (E_1) are used; the fourth,
$w^{-1}R_4R_1^{-1}=U_1$, is automatic: $w^{-1}(zD)^3R_1R_1^{-1}=z^{-1}D^{-1}\cdot(z^4D^4)w^{-1}=z^{-1}D^{-1}=U_1$
because $D^4=1$ ($D=1+(\mu-1)\Pi$ with $\mu=-i$, so $D^k=1+(\mu^k-1)\Pi$). Likewise $R_3^d=R_4^d=1$ are not
needed as hypotheses.

*Remark 3.11 (d = 2).* The argument also works for $d=2$ (then $m_2=2=d$ and (T_2) is automatic since
$z^2U_2=i$), but no certificate for $d=2$ is used here.

---------------------------------------------------------------------------------------------------------

## 4. The self-testing theorem

**Theorem 4.1.** Let $d\ge3$ and assume (OPT_d) and (RIG_d) (e.g. $d\in\{3,\dots,9\}$, Section 6). Let
$\{A_x^a\},\{B_y^b\}$ be $d$-outcome PVMs on $\mathbb C^D$ whose correlations on $|\Phi_D\rangle$ attain
$I_d=I_{\rm ME}(d)$ (equivalently, maximise $I_d$ over all maximally entangled projective strategies of all local
dimensions). Then:

1. $D=dM$ for an integer $M\ge1$;
2. there is a unitary $V:\mathbb C^D\to\mathbb C^d\otimes\mathbb C^M$ such that, with $\bar V$ the entrywise complex
   conjugate (w.r.t. the computational bases $\{|k\rangle\}$ and $\{|j\rangle\otimes|\mu\rangle\}$),
   * $(V\otimes\bar V)|\Phi_D\rangle=|\Phi_d\rangle\otimes|\Phi_M\rangle$ (after reordering
     $\mathbb C^d\otimes\mathbb C^M\otimes\mathbb C^d\otimes\mathbb C^M\to(\mathbb C^d\otimes\mathbb C^d)\otimes(\mathbb C^M\otimes\mathbb C^M)$);
   * $VA_x^aV^*=A_x^{a,\rm DKZ}\otimes1_M$ and $\bar VB_y^b\bar V^{*}=B_y^{b,\rm DKZ}\otimes1_M$ for all $x,y,a,b$;
3. $V$ is unique up to $V\mapsto(1_d\otimes u)V$, $u\in U(M)$.

Conversely, for every $M$, $(A^{\rm DKZ}\otimes1_M,B^{\rm DKZ}\otimes1_M)$ on $|\Phi_{dM}\rangle$ attains $I_{\rm ME}(d)$.

*Proof.* Let $R=(A_1,B_1^T,A_2,B_2^T)\in\mathcal U_D$. By (OPT_d) and the strict monotonicity of
$I_d=4-2S/(d-1)$, $R$ is optimal, so (RIG_d) gives (E_1') and (T_n), $n\le m_d$. Theorem 3.9 gives $W$ with
$R_k=W(R^c_k\otimes1_M)W^*$; with $W_0$ of Lemma 1.7(e), $R^c_k=W_0^*R^{\rm DKZ}_kW_0$. Put
$V:=(W_0\otimes1_M)W^*$. Then $VR_kV^*=R^{\rm DKZ}_k\otimes1_M$. For Alice ($k=1,3$) this is
$VA_xV^*=A_x^{\rm DKZ}\otimes1_M$, and the PVM elements follow from $A^a=\frac1d\sum_nw^{-na}A^n$. For Bob,
$VB_y^TV^*=(B_y^{\rm DKZ})^T\otimes1_M=(B_y^{\rm DKZ}\otimes1_M)^T$; transposing, $\bar VB_yV^T=B_y^{\rm DKZ}\otimes1_M$,
and $V^T=\bar V^{*}$. State: $(V\otimes\bar V)\sum_k|k\rangle|k\rangle=\sum_{\alpha,\beta}(\sum_kV_{\alpha k}\overline{V_{\beta k}})
|\alpha\rangle|\beta\rangle=\sum_{\alpha,\beta}(VV^*)_{\alpha\beta}|\alpha\rangle|\beta\rangle=\sum_\alpha|\alpha\rangle|\alpha\rangle$,
with $\alpha$ running over $\{|j\rangle\otimes|\mu\rangle\}$; normalise with $D=dM$. Uniqueness: if $V'$ is another
such unitary, $C:=V'V^*$ commutes with all $R^{\rm DKZ}_k\otimes1_M$; these generate $M_d(\mathbb C)\otimes1_M$
(Lemma 1.7(f)), whose commutant is $1_d\otimes M_M(\mathbb C)$, so $C=1_d\otimes u$. Converse:
$\mathrm{tr}_{dM}(x\otimes1_M)=\mathrm{tr}_d(x)$, so $S(R^{\rm DKZ}\otimes1_M)=S(R^{\rm DKZ})$, and
$(B\otimes1_M)^T=B^T\otimes1_M$. $\square$

**What "uniqueness" means precisely.**

* *Up to local unitaries:* Alice's frame is changed by $V$, Bob's by $\bar V$. Such pairs are exactly the local
  unitaries mapping $|\Phi_D\rangle$ to a maximally entangled state of the product form above; no other
  freedom (e.g. no local isometry into a larger space) is needed.
* *On the support of the state:* the reduced states of $|\Phi_D\rangle$ are $1/D$, of full rank, so the support is
  all of $\mathbb C^D$: the operators themselves (not only their action on the state) are determined. This is
  where the full-rank hypothesis is used (Lemma 1.6).
* *Self-testing form:* measurements $=\,$DKZ$\,\otimes1_M$, junk state $|\Phi_M\rangle$ (maximally entangled;
  forced, since the total state is); the multiplicity space $\mathbb C^M$ carries no information about the
  outcomes.
* *Dimension:* $d\mid D$ (Corollary 5.1 below).
* *DKZ form:* in the frames $V,\bar V$ every measurement basis is (Fourier basis dressed by a diagonal unitary)
  $\otimes$ (arbitrary orthonormal basis of $\mathbb C^M$); the PVM elements $A^{a,\rm DKZ}\otimes1_M$ have rank $M$.
* *Rigidity of the frame:* $V$ is unique up to $1_d\otimes u$; for $M=1$ up to a global phase. Alice's two
  observables alone already generate $M_d(\mathbb C)\otimes1_M$: in the canonical frame
  $R^c_3R^{c\,-1}_1=(zD^c)^2=z^2\,\mathrm{diag}(-1,1,\dots,1)=z^2(1-2\Pi^c)$, so the algebra generated by $R^c_1=X$
  and $R^c_3$ contains $\Pi^c$ and hence all matrix units. The same holds for Bob's pair $R^c_2,R^c_4$
  ($R^c_4R^{c\,-1}_2=(zD^c)^2$).

---------------------------------------------------------------------------------------------------------

## 5. Consequence: $d\nmid D$ gives a strictly smaller value

**Corollary 5.1.** Assume (OPT_d) and (RIG_d). For every $D$ with $d\nmid D$,
$$\max\{\,I_d(A,B;\Phi_D):\ A_1,A_2,B_1,B_2\ \text{$d$-outcome PVMs on }\mathbb C^D\,\}\;<\;I_{\rm ME}(d).$$

*Proof.* The set $\mathcal P_D$ of $d$-outcome PVMs on $\mathbb C^D$ is a closed (defined by the polynomial
equations $P_a=P_a^*=P_a^2$, $\sum_aP_a=1$) and bounded ($\|P_a\|\le1$) subset of $M_D(\mathbb C)^d$, hence
compact, and $I_d$ is a continuous (polynomial) function on $\mathcal P_D^4$. So the maximum is attained by some
strategy. If it equalled $I_{\rm ME}(d)$, that strategy would be optimal and Theorem 4.1 would give $d\mid D$. $\square$

*Remark 5.2 (no uniform gap).* For $D=kd+r$, $0<r<d$, take $(R^{\rm DKZ}\otimes1_k)\oplus R'$ with any $R'\in\mathcal U_r$.
Since $\mathrm{tr}_D(x\oplus y)=\frac{kd}D\mathrm{tr}_{kd}(x)+\frac rD\mathrm{tr}_r(y)$ and $S$ is affine in the
traces, $S=\frac{kd}D\lambda_d+\frac rDS(R')$, so $I_{\rm ME}(d)-\max_D\le\frac rD\cdot\frac{2(S(R')-\lambda_d)}{d-1}\to0$.
Thus $\sup_{d\nmid D}\max_D=I_{\rm ME}(d)$ is not attained; a quantitative gap for fixed $D$ would require a
robust version of the rigidity theorem (not attempted). Numerical illustration (lab LOG, `maxent_primal.py`,
$d=3$): best $S$ found $1.34530$ ($D=4$), $1.22771$ ($D=5$) versus $S_{\rm ME}(3)=1.12707$ ($D=6$, DKZ$\otimes1_2$).

---------------------------------------------------------------------------------------------------------

## 6. Verified instances $d=3,\dots,9$

The certificates `../../oqp27B/cert_tracial_d{d}.pkl` (level-2 tracial moment relaxation; full word set for
$d\le6$, adjacent-pair word set for $d=7,8,9$) are symmetric SOS certificates in the sense of Definition 2.1
(blocks labelled by $\mathbb Z_d$-charge $q$ and $\rho$-eigenvalue $\mu=\zeta_{4d}^{q+dt}$; elements
$e=\sum_{k=0}^3\mu^{-k}\rho^k(u)$ for words $u$ of length $\le2$; $X_b=N_bY_bN_b^*$).

| check | lab tool (log) | independent re-check `check_support.py` (log `check_support_d3to9.log`) |
|---|---|---|
| (Id) modulo commutators and $G$ | `verify_tracial.py` (FINAL_verify_tracial.log) | (P1): own canonicalisation with $G$ enumerated explicitly ($|G|=16d$); class counts 23, 59, 122, 220, 229, 344, 497 agree |
| $S_{\rm op}$ $G$-invariant mod commutators | `verify_tracial.py` | (P1), all $\alpha\in G$ |
| $Y_b\succ0$ | exact LDL + interval signs | (S1) interval Cholesky of real embedding, 320 bits |
| $\lambda=S(R^{\rm DKZ})$ | `verify_tracial.py` | (S0) own exact evaluation |
| $M_{n,t}\in\mathcal S$, all $n\le d-1$ | `verify_rigidity.py` (whole blocks $(0,t)$, $t\ne0$, in support) | (S2) exact membership with explicit coefficients, re-verified |
| $T_n\in\mathcal S$, all $n\le d-1$ | `verify_rigidity.py` | (S3) |
| $V_n\in\mathcal S$ (non-symmetrised) | - | (S3') |
| consistency | - | (S4) all targets and support polynomials vanish on DKZ; (S5) $W_0$; negative controls rejected |

Hence conditions (C1) and (C2) of Theorem 2.6 hold for $d=3,\dots,9$ (for all $n\le d-1$, more than the needed
$n\le m_d$), and Theorems 4.1 and Corollary 5.1 hold unconditionally for these $d$.

Re-run: `python -B check_support.py 3 4 5 6 7 8 9 --identity` (about one minute in total).

---------------------------------------------------------------------------------------------------------

## 7. What an all-$d$ certificate must provide (target for workstreams A and B)

For all-$d$ rigidity it suffices to establish, for each $d$:

1. **(OPT_d)** the bound $S\ge S(R^{\rm DKZ})$ over all strategies of all local dimensions, by any method; and
2. **(RIG_d)**: optimal strategies satisfy (E_1') and (T_n) for $n\le\lfloor d/2\rfloor+1$ on the first link.

For an SOS-type proof $S_{\rm op}-\lambda=\sum_kf_k^*f_k$ (mod commutators, and mod $G$-relations if
symmetry-reduced), (RIG_d) follows (Proposition 2.4, 2.5) as soon as the two-sided ideal generated by the
$G$-images of the $f_k$ contains either
* (C1) $M_{n,t}$ ($t=1,2,3$) and $T_n$ for $1\le n\le\lfloor d/2\rfloor+1$, or
* (C2) the two link differences $\ell^{(0)}_1-\ell^{(1)}_1$, $\ell^{(1)}_1-\ell^{(2)}_1$ and
  $V_n=(\ell_n^{(0)})^{-1}-(1-i)z^n-iz^{2n}\ell^{(0)}_n$ for $1\le n\le\lfloor d/2\rfloor+1$.

For an operator-inequality proof (workstream B), the same holds if the equality case of each inequality
used forces these relations.

**Face-reduced certificates give rigidity for free.** Suppose a certificate is built as in the lab: blocks of
symmetry-adapted elements, $X_b=N_bY_bN_b^*$ with $N_b$ a basis of the kernel of the block's moment matrix
evaluated on *all* $G$-images of $R^{\rm DKZ}$, and $Y_b\succ0$. Then $\mathcal S_b$ is exactly the set of block
polynomials vanishing on all $G$-images of $R^{\rm DKZ}$. By Lemma 1.7(c),(f) and Lemma 1.9, $M_{n,t}$, $T_n$
and every $\rho$-Fourier component $\sum_ki^{-tk}\rho^k(V_n)$ of $V_n$ vanish on all these images, for every $d$
and $n$. Hence, as soon as these polynomials lie in the linear span of the block elements (true whenever the
word set contains $1$, $r_k^nr_{k+1}^{-n}$ and $r_{k+1}^nr_k^{-n}$ for the four adjacent pairs, as both lab word
sets do), (C1) and (C2) hold automatically: for face-reduced certificates with positive definite $Y_b$ the
rigidity conclusion costs nothing beyond the bound itself. If instead a certificate uses a Gram matrix of
lower rank on the face, the membership conditions must be checked (e.g. with `check_support.py`).

### 7.1 Bridges to the formulations of workstreams A and B

**A's chain modes.** For $R\in\mathcal U_D$ write (0-indexed as in A) $R_0,\dots,R_3$ for $R_1,\dots,R_4$,
$R_4:=wR_0$, $c:=z^n$ and $W_t:=W(n,t):=\frac14\sum_{k=0}^3(c\,i^t)^{-k}R_k^n$ ($t\in\mathbb Z_4$). $R$ is
*grade-sharp at $n$* if $W(n,1)=W(n,2)=0$.

**Lemma 7.1 (grade-sharp $\Rightarrow$ (E_n), (T_n)).** If $R$ is grade-sharp at $n$, then (E_n) and (T_n) hold,
with $\Pi^{(n)}=W_3W_3^*$. Consequently, if $R$ is grade-sharp at $n=1,\dots,\lfloor d/2\rfloor+1$, then $R$ is
unitarily equivalent to $R^{\rm DKZ}\otimes1_M$ (Theorems 3.8, 3.9), unconditionally in $d$.

*Proof.* (i) Inversion: $R_k^n=\sum_t(c\,i^t)^kW_t$ for $k=0,\dots,4$ (DFT on $\mathbb Z_4$; for $k=4$ use $c^4=w^n$).
(ii) $\sum_tW_tW_{t+r}^*=\delta_{r0}1$ and $\sum_tW_t^*W_{t+r}=\delta_{r0}1$: e.g.
$\sum_tW_tW_{t+r}^*=\frac1{16}\sum_{k,l}c^{l-k}i^{rl}\big(\sum_ti^{t(l-k)}\big)R_k^nR_l^{-n}=\frac14\sum_{k=0}^3i^{rk}=\delta_{r0}$
($|c|=1$, $R_l^{-n}=(R_l^n)^*$); the other identity is the same computation with $R_k^{-n}R_l^{n}$ (A-I1).
(iii) With $W_1=W_2=0$: $r=3$ gives $W_0W_3^*=0$ and $W_0^*W_3=0$; $r=0$ gives $W_0W_0^*+W_3W_3^*=1$ and
$W_0^*W_0+W_3^*W_3=1$. Hence $W_3W_3^*W_3=W_3(1-W_0^*W_0)=W_3-(W_0W_3^*)^*W_0=W_3$, so $E:=W_3W_3^*$ is a
projection and $W_0W_0^*=1-E$. (iv) By (i), $R_k^n=c^k(W_0+(-i)^kW_3)$, so for $k=0,\dots,3$
$$L^{(k)}_n=R_k^nR_{k+1}^{-n}=c^{-1}\big[W_0W_0^*+i^{k+1}W_0W_3^*+(-i)^kW_3W_0^*+(-i)^ki^{k+1}W_3W_3^*\big]
=c^{-1}\big[(1-E)+iE\big],$$
independent of $k$ (this is (E_n); note $L^{(3)}_n=R_3^n(wR_0)^{-n}$ is our fourth link), and
$z^nU_n=(1-E)+iE=\mathcal Q(E)$ is two-valued (Lemma 3.1): (T_n). $\square$

Conversely $R^{\rm DKZ}$ is grade-sharp at every $n$ (A's "counter walk" description; checked numerically here
for $d=3,\dots,8$, together with $W_0W_3^*=0$). So, given (OPT_d), **(RIG_d) is equivalent to: every optimal
strategy is grade-sharp at $n=1,\dots,\lfloor d/2\rfloor+1$**, and A's Lemma A-L1 (grade-sharp $\Rightarrow S=\lambda_d$)
is consistent with Lemma 7.1: grade-sharp strategies are exactly the DKZ$\,\otimes1_M$ strategies.

**B's reduced variables.** B-T1 represents a (twirled) max-ent strategy as $R_1=X\otimes1_M$,
$R_{i+1}=WR_iW^*$, $W=\sum_k|k\rangle\langle k|\otimes W_k$, $W_k=z^kV_k$, $V_k^4=1$ ($k\in\mathbb Z_d$).

**Lemma 7.2 (translation of (E_1'), (T_n)).** In this representation:
(a) $R_1R_2^{-1}=R_2R_3^{-1}=R_3R_4^{-1}$ iff $V_kV_{k-1}=V_{k-1}V_k$ for all $k\in\mathbb Z_d$ (cyclically adjacent
$V$'s commute);
(b) (T_n) holds iff for every $k$ the unitary $V_{k-n}V_k^*$ ($k\ge n$), resp. $i\,V_{k-n+d}V_k^*$ ($k<n$), has spectrum
in $\{1,i\}$; in B's notation $E_k=V_k^*$ these are the relations $g_{jk}=E_j^*E_k+iE_k^*E_j-(1+i)=0$ for the pairs
$j<k$ with $k-j\in\{n,d-n\}$ (Lemma 3.1(ii)).

*Proof.* $L^{(0)}=R_1R_2^{-1}=\hat XW\hat X^*W^*=\sum_k|k\rangle\langle k|\otimes W_{k-1}W_k^*$ ($\hat X=X\otimes1$), and
$L^{(j+1)}=WL^{(j)}W^*$. Hence $L^{(0)}=L^{(1)}$ iff $[W,L^{(0)}]=0$ iff $W_k$ commutes with $W_{k-1}W_k^*$ for all
$k$ iff $W_kW_{k-1}=W_{k-1}W_k$; and $L^{(1)}=L^{(2)}$ is the $W$-conjugate of the same condition. (b)
$U_n=\hat X^nW\hat X^{-n}W^*=\sum_k|k\rangle\langle k|\otimes W_{k-n}W_k^*$ and $z^nW_{k-n}W_k^*=V_{k-n}V_k^*$ for $k\ge n$,
$=z^dV_{k-n+d}V_k^*=iV_{k-n+d}V_k^*$ for $k<n$. For $k<n$ the condition says $V_{k-n+d}V_k^*$ has spectrum in
$\{1,-i\}$, i.e. $V_kV_{k-n+d}^*=E_k^*E_{k-n+d}$ has spectrum in $\{1,i\}$. $\square$

Since B's twirled strategy $\tilde R=\bigoplus_{k=0}^{4d-1}\rho^k(R)$ contains $R$ as a direct summand, relations
that hold for $\tilde R$ hold for $R$. **Target for B:** (RIG_d) follows if every optimal family $(V_k)$ (any $M$)
satisfies (a) cyclic adjacent commutation and (b) the relations $g_{jk}$ for $k-j\in\{n,d-n\}$,
$1\le n\le\lfloor d/2\rfloor+1$. All of these vanish on B's classical one-step optima, so for a face-reduced
certificate of the reduced problem with positive definite Gram blocks they lie in the support as soon as the
words $E_jE_k$, $E_kE_j$, $E_j^*E_k$ are in the word set (commutators included).

---------------------------------------------------------------------------------------------------------

## 8. Scope and limitations

* **Projective measurements, full-rank maximally entangled state.** The theorem concerns PVMs on
  $\mathbb C^D\otimes\mathbb C^D$ and $|\Phi_D\rangle$. POVMs are not covered (the certificates bound
  $S$ only over unitaries of order $d$); nor are states that are maximally entangled on a proper subspace of
  larger local spaces (their compressed measurements are POVMs).
* **Tracial von Neumann algebras (remark).** Theorem 3.8 holds in any C*-algebra. Proposition 2.4 holds for every
  tracial state; if $(\mathcal N,\tau)$ is a von Neumann algebra with faithful normal tracial state and
  $R_k\in\mathcal N$ attain $\lambda$ (use $\mathcal N^{\rm op}$ in place of the transpose for anti-automorphisms, and
  the faithful trace $|G|^{-1}\sum\tau$ on the orbit direct sum), Theorem 3.8 yields matrix units $e_{jk}\in\mathcal N$
  and $\mathcal N\cong M_d(\mathbb C)\otimes e_{00}\mathcal Ne_{00}$ via $x\mapsto[e_{0j}xe_{k0}]_{jk}$, under which
  $R_k\mapsto R^c_k\otimes1$. So the rigidity statement extends to this commuting-operator analogue of
  maximally entangled strategies.
* **Exact, not robust.** No statement is made about near-optimal strategies.
* **All $d$ is conditional** on (OPT_d) + (RIG_d), i.e. on the existence of an all-$d$ certificate with the
  support properties of Section 7. Theorem A (the algebraic part) is proved for all $d$.

---------------------------------------------------------------------------------------------------------

## 9. Changes with respect to the original sketch (lab LOG section 6)

1. Faithful-trace step made rigorous: the SOS identity holds only modulo symmetry relations, so it is evaluated
   on the orbit direct sum $R^\oplus$ (Lemma 2.2, Proposition 2.4); conclusion for $R$ itself and for all $G$-images.
2. "(E_n)": precisely, vanishing of the three nontrivial chain modes (Proposition 2.5(a)). (E_n) for $n\ge2$ is used
   only to de-symmetrise $T_n$, and is not needed at all when $V_n\in\mathcal S$ (verified for $d=3..9$).
3. Only three of the four equalities of (E_1) are needed (Remark 3.10).
4. Jordan-type lemma: purely algebraic proof valid in any C*-algebra, covering every degenerate configuration at
   once (Lemma 3.2; spectral picture in Remark 3.2').
5. Two-valuedness is needed only for $n\le\lfloor d/2\rfloor+1$, not for all $n\le d-1$ (Lemma 3.5); the products are
   exactly $\mathcal Q(P_{n-1})\cdots\mathcal Q(P_0)$ (Lemma 3.3), and each induction step applies the lemma to a
   product of two quarter-phase unitaries $\mathcal Q(P_k)\mathcal Q(E_k)$.
6. Imprimitivity via matrix units (Lemma 3.7): works in any unital $*$-algebra; the finite-dimensional unitary
   and $d\mid D$ follow (Theorem 3.9).
7. Precise self-testing statement including Bob's conjugate frame $\bar V$, the junk state, and the uniqueness of $V$
   (Theorem 4.1).
8. New: the symmetry group maps the DKZ class to itself (Lemma 1.9), which shows that face-reduced certificates
   with positive definite Gram blocks automatically carry the rigidity support conditions (Section 7).

---------------------------------------------------------------------------------------------------------

## 10. Formalisation status (Lean 4 + Mathlib, `formal-conjectures/CGLMPRigidity/`)

Machine-checked, no `sorry`/`admit`/`axiom`/`native_decide`; `#print axioms` of every theorem below gives only
`[propext, Classical.choice, Quot.sound]` (`Axioms.lean`). Check: from `formal-conjectures/`,
`bash CGLMPRigidity/check.sh` (about 3 minutes; files `Rigidity`, `Links`, `Model`, `Bridge`, `Axioms`;
log `CGLMPRigidity/check_all.log`).

Setting of the formal statements: `A` any unital `*`-algebra over `ℂ` (`Ring`, `StarRing`, `Algebra ℂ`,
`StarModule ℂ`) with the hypothesis `∀ x, star x * x = 0 → x = 0`; unitaries as `star U * U = 1 ∧ U * star U = 1`;
projections as Mathlib's `IsStarProjection`; `R⁻¹` written `star R`.

| paper | Lean (file) | content |
|---|---|---|
| Lemma 3.1 | `exists_proj_of_two_valued`, `two_valued_ii` (Rigidity), `qp_two_valued` (Model) | two-valued unitary = quarter-phase unitary of a projection |
| Lemma 3.2 | `jordan` (Rigidity) | $(X-1)(X-i)=0$ for $X=\mathcal Q(P)\mathcal Q(Q)$ implies $PQ=0$ |
| Lemma 3.3 | `prod_formula`, `zpow_link` (Rigidity) | product formula |
| Lemma 3.4 | `orth_upto` (Rigidity) | orthogonality induction from $T_2,\dots,T_{K+1}$ |
| Lemma 3.5 | `orth_all` (Rigidity) | cyclic transport, $K=\lfloor d/2\rfloor$ suffices |
| Lemma 3.6-3.7 | inside `rigidity`; `mu_mul`, `mu_star`, `R_eq_sum` | resolution of identity; matrix units, $R_1$ = cyclic shift |
| Theorem 3.8 | `rigidity` (Rigidity) | hypotheses (E_1'), (T_n) for $n\le d/2+1$ $\Rightarrow$ projection $p$, $R_{k+1}=z(1-(1+i)p)R_k$, $d\times d$ matrix units $R_1^jpR_1^{\star k}$ summing to 1, $R_1=\sum_j e_{j+1,j}$ |
| Theorem 3.9 ($d\mid D$) | `dvd_of_rigidity` (Rigidity) | for `Matrix (Fin D) (Fin D) ℂ`: $d\mid D$ (trace of $p$ = rank) |
| Proposition 2.5 | `eq_of_modes`, `two_valued_iff_linear`, `V_eq` (Links) | chain modes $\Rightarrow$ equal links; $V_n=0\iff(T_n)$ |
| Theorem 2.6 (C2) $\Rightarrow$ 3.8 | `rigidity_of_support` (Links) | vanishing of $M_{1,t}$ and $V_n$ ($n\le d/2+1$) $\Rightarrow$ DKZ structure |
| non-vacuity | `canonical_model` (Model) | for every $d\ge1$ the canonical DKZ tuple in $M_d(\mathbb C)$ with $z=e^{i\pi/(2d)}$ satisfies all hypotheses of `rigidity` (for all $n\le d$) |
| Lemma 7.1 (core) | `grade_sharp_core` (Bridge) | if $u_k=Y_0+(-i)^kY_3$ ($k=0..3$) are unitary, $E=Y_3Y_3^\star$ is a projection and all four links $u_ku_{k+1}^\star$ equal $\mathcal Q(E)$ |

Not formalised: the certificate side (Definition 2.1, Lemma 2.2, Proposition 2.4: group algebra of
$\mathbb Z_d^{*4}$, symmetry group, faithful trace), which is covered by the exact computer checks of Section 6;
the explicit unitary $W$ of Theorem 3.9 (only its consequence $d\mid D$ is formal); the translation to PVMs and
$V\otimes\bar V$ in Theorem 4.1; Corollary 5.1 (compactness); Lemma 1.9.
