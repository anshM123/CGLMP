import CGLMPRigidity.Model

/-!
# Grade-sharp strategies satisfy the rigidity hypotheses (Lemma 7.1 of `RIGIDITY.md`)

In the chain-mode coordinates `W(n, t)` of workstream A, a strategy is grade-sharp at `n` if
`W(n, 1) = W(n, 2) = 0`; then `R_k ^ n = z ^ (n k) (Y₀ + (-i) ^ k Y₃)` for `k = 0, …, 4` with
`Y₀ = W(n, 0)`, `Y₃ = W(n, 3)` (and `R₄ = w R₀`). `CGLMPRigidity.grade_sharp_core` shows: if the four
elements `u_k = Y₀ + (-i) ^ k Y₃` are unitary, then `E = Y₃ Y₃⋆` is a star projection and every link
`u_k u_{k+1}⋆` (indices mod 4) equals `qp E = 1 + (i - 1) E`. Since the phases `z ^ (n k)` contribute the
common factor `z ^ (-n)`, this is exactly: all four links are equal (`(E_n)`) and `z ^ n U_n` is two-valued
(`(T_n)`).
-/

set_option linter.unusedSectionVars false

namespace CGLMPRigidity

open Complex

variable {A : Type*} [Ring A] [StarRing A] [Algebra ℂ A] [StarModule ℂ A]

lemma expand_mul_star (Y₀ Y₃ : A) (a : ℂ) (ha : a * star a = 1) :
    (Y₀ + a • Y₃) * star (Y₀ + a • Y₃)
      = Y₀ * star Y₀ + star a • (Y₀ * star Y₃) + a • (Y₃ * star Y₀) + Y₃ * star Y₃ := by
  have ha' : star a * a = 1 := by rw [mul_comm]; exact ha
  simp only [star_add, star_smul, mul_add, add_mul, smul_mul_assoc, mul_smul_comm, smul_add,
    smul_smul, ha', one_smul]
  abel

lemma expand_star_mul (Y₀ Y₃ : A) (a : ℂ) (ha : star a * a = 1) :
    star (Y₀ + a • Y₃) * (Y₀ + a • Y₃)
      = star Y₀ * Y₀ + a • (star Y₀ * Y₃) + star a • (star Y₃ * Y₀) + star Y₃ * Y₃ := by
  have ha' : a * star a = 1 := by rw [mul_comm]; exact ha
  simp only [star_add, star_smul, mul_add, add_mul, smul_mul_assoc, mul_smul_comm, smul_add,
    smul_smul, ha', one_smul]
  abel

lemma expand_link (Y₀ Y₃ : A) (a b : ℂ) :
    (Y₀ + a • Y₃) * star (Y₀ + b • Y₃)
      = Y₀ * star Y₀ + star b • (Y₀ * star Y₃) + a • (Y₃ * star Y₀) + (a * star b) • (Y₃ * star Y₃) := by
  simp only [star_add, star_smul, mul_add, add_mul, smul_mul_assoc, mul_smul_comm, smul_add,
    smul_smul]
  module

/-- **Lemma 7.1 (core).** If `u_k = Y₀ + (-i) ^ k Y₃` (`k = 0, 1, 2, 3`) are unitary, then `E = Y₃ Y₃⋆` is a
star projection and all four cyclic links `u_k u_{k+1}⋆` equal `qp E`. -/
theorem grade_sharp_core {Y₀ Y₃ : A}
    (h0 : IsU (Y₀ + (1 : ℂ) • Y₃)) (h1 : IsU (Y₀ + (-I) • Y₃))
    (h2 : IsU (Y₀ + (-1 : ℂ) • Y₃)) (h3 : IsU (Y₀ + I • Y₃)) :
    IsStarProjection (Y₃ * star Y₃) ∧
      (Y₀ + (1 : ℂ) • Y₃) * star (Y₀ + (-I) • Y₃) = qp (Y₃ * star Y₃) ∧
      (Y₀ + (-I) • Y₃) * star (Y₀ + (-1 : ℂ) • Y₃) = qp (Y₃ * star Y₃) ∧
      (Y₀ + (-1 : ℂ) • Y₃) * star (Y₀ + I • Y₃) = qp (Y₃ * star Y₃) ∧
      (Y₀ + I • Y₃) * star (Y₀ + (1 : ℂ) • Y₃) = qp (Y₃ * star Y₃) := by
  have s1 : (1 : ℂ) * star (1 : ℂ) = 1 ∧ star (1 : ℂ) * 1 = 1 := by constructor <;> cring
  have sI : (-I) * star (-I) = 1 ∧ star (-I) * (-I) = 1 := by constructor <;> cring
  have sm : (-1 : ℂ) * star (-1 : ℂ) = 1 ∧ star (-1 : ℂ) * (-1) = 1 := by constructor <;> cring
  have sJ : I * star I = 1 ∧ star I * I = 1 := by constructor <;> cring
  -- the four relations `u u⋆ = 1`
  have k0 := h0.2
  have k1 := h1.2
  have k2 := h2.2
  have k3 := h3.2
  rw [expand_mul_star _ _ _ s1.1] at k0
  rw [expand_mul_star _ _ _ sI.1] at k1
  rw [expand_mul_star _ _ _ sm.1] at k2
  rw [expand_mul_star _ _ _ sJ.1] at k3
  -- the four relations `u⋆ u = 1`
  have l0 := h0.1
  have l1 := h1.1
  have l2 := h2.1
  have l3 := h3.1
  rw [expand_star_mul _ _ _ s1.2] at l0
  rw [expand_star_mul _ _ _ sI.2] at l1
  rw [expand_star_mul _ _ _ sm.2] at l2
  rw [expand_star_mul _ _ _ sJ.2] at l3
  have cancel4 : ∀ x : A, (4 : ℂ) • x = 0 → x = 0 := by
    intro x hx
    rcases smul_eq_zero.1 hx with h | h
    · norm_num at h
    · exact h
  -- `Y₀ Y₃⋆ = 0`
  have hA : Y₀ * star Y₃ = 0 := by
    apply cancel4
    linear_combination (norm := skip) k0 - k2 - I • k1 + I • k3
    match_scalars <;> cring
  have hB : Y₃ * star Y₀ = 0 := by
    have := congrArg star hA
    rwa [star_mul, star_star, star_zero] at this
  -- `Y₀⋆ Y₃ = 0`
  have hA' : star Y₀ * Y₃ = 0 := by
    apply cancel4
    linear_combination (norm := skip) l0 - l2 + I • l1 - I • l3
    match_scalars <;> cring
  have hB' : star Y₃ * Y₀ = 0 := by
    have := congrArg star hA'
    rwa [star_mul, star_star, star_zero] at this
  have hP : Y₀ * star Y₀ + Y₃ * star Y₃ = 1 := by
    rw [hA, hB, smul_zero, smul_zero, add_zero, add_zero] at k0
    exact k0
  have hP' : star Y₀ * Y₀ + star Y₃ * Y₃ = 1 := by
    rw [hA', hB', smul_zero, smul_zero, add_zero, add_zero] at l0
    exact l0
  -- `E = Y₃ Y₃⋆` is a star projection
  have hE : IsStarProjection (Y₃ * star Y₃) := by
    rw [isStarProjection_iff']
    constructor
    · have e : star Y₃ * Y₃ = 1 - star Y₀ * Y₀ := by rw [← hP']; abel
      calc Y₃ * star Y₃ * (Y₃ * star Y₃) = Y₃ * (star Y₃ * Y₃) * star Y₃ := by
            simp only [mul_assoc]
        _ = Y₃ * star Y₃ - (Y₃ * star Y₀) * (Y₀ * star Y₃) := by
            rw [e, mul_sub, sub_mul, mul_one]; simp only [mul_assoc]
        _ = Y₃ * star Y₃ := by rw [hB, zero_mul, sub_zero]
    · rw [star_mul, star_star]
  -- the links
  have hY0 : Y₀ * star Y₀ = 1 - Y₃ * star Y₃ := by rw [← hP]; abel
  have link : ∀ a b : ℂ, a * star b = I →
      (Y₀ + a • Y₃) * star (Y₀ + b • Y₃) = qp (Y₃ * star Y₃) := by
    intro a b hab
    rw [expand_link, hA, hB, smul_zero, smul_zero, add_zero, add_zero, hab, hY0, qp]
    module
  refine ⟨hE, link _ _ (by cring), link _ _ (by cring), link _ _ (by cring), link _ _ (by cring)⟩

end CGLMPRigidity
