import NewenSpec.Float
import NewenSpec.Nonbonded
import NewenTheory.Interactions.SoftCore

/-!
# Model of the perturbed (soft-core) nonbonded kernels of `crates/newen-forces`

Rust, all in `src/nonbonded/perturbed.rs` (spec keys in `spec/SPEC.toml`):

* `nonbonded.perturbed_lambda_params`: `PerturbedLambdaParams::from_lambda`, the λ table of one
  pair of energy groups (gromosXX `set_lambda`), with `powf` for `(1 − λ)ⁿ`, `λⁿ`;
* `nonbonded.perturbed_lj_crf_pair`: `perturbed_lj_crf_pair`, one perturbed LJ + CRF pair;
* `nonbonded.rf_soft_interaction`: `rf_soft_interaction`, the reaction field of an excluded
  perturbed pair (a `PERTATOMPAIR` state);
* `nonbonded.perturbed_rf_excluded_pair`: `perturbed_rf_excluded_pair`, the same term as the
  excluded-pair correction rounds it;
* `nonbonded.perturbed_rf_self_term`: `perturbed_rf_self_term`, the self term of a perturbed atom;
* `nonbonded.lj_crf_soft_interaction_ext`, `nonbonded.rf_soft_interaction_ext`: the unweighted
  end-state energies and their λ derivatives that PRECALCLAM books at every λ_k.

Theory: `NewenTheory/Interactions/SoftCore.lean`.

As in `NewenSpec/Nonbonded.lean`, each kernel is written once over any type `α` with the
operations the Rust uses, operation by operation in the Rust order, and used at `ℝ` (the
`*_real_eq` theorems: equal to the theory) and at `Float` (`lake exe spec-vectors`, compared bit
for bit by `crates/newen-forces/tests/softcore_spec.rs`). The kernels use `+ − × ÷`, `sqrt` and,
in the λ table only, `pow` (the `Powf` class of `NewenSpec/Float.lean`).

The central result is `perturbedLjCrfPairModel_dvdl_is_deriv`: over ℝ, the `dV/dλ` the pair
kernel returns is the derivative in λ of the energy it returns, with the λ table built by the
table kernel from individual λ that are differentiable in λ.

What is not proved: the rounding error between the two instantiations.
-/

namespace NewenSpec.SoftCore

open NewenSpec.Nonbonded NewenTheory.Interactions.SoftCore
open NewenTheory.Interactions.LennardJones NewenTheory.Interactions.ReactionField

/-- `PerturbedLambdaParams`: the λ table one pair reads. -/
structure LambdaParams (α : Type) where
  /-- `(1 − λ_LJ)ⁿ` -/
  aLjN : α
  /-- `λ_LJⁿ` -/
  bLjN : α
  /-- `(1 − λ_CRF)ⁿ` -/
  aCrfN : α
  /-- `λ_CRFⁿ` -/
  bCrfN : α
  /-- `(1 − λ_LJ)ⁿ⁻¹ · λ_LJ′` -/
  aLjN1 : α
  /-- `λ_LJⁿ⁻¹ · λ_LJ′` -/
  bLjN1 : α
  /-- `(1 − λ_CRF)ⁿ⁻¹ · λ_CRF′` -/
  aCrfN1 : α
  /-- `λ_CRFⁿ⁻¹ · λ_CRF′` -/
  bCrfN1 : α
  /-- `(1 − λ_LJs) · λ_LJs′` -/
  aLjs : α
  /-- `λ_LJs · λ_LJs′` -/
  bLjs : α
  /-- `(1 − λ_LJs)²` -/
  aLjs2 : α
  /-- `λ_LJs²` -/
  bLjs2 : α
  /-- `(1 − λ_CRFs) · λ_CRFs′` -/
  aCrfs : α
  /-- `λ_CRFs · λ_CRFs′` -/
  bCrfs : α
  /-- `(1 − λ_CRFs)²` -/
  aCrfs2 : α
  /-- `λ_CRFs²` -/
  bCrfs2 : α
  /-- NLAM as a number (`lambda_exp as f64`). -/
  nExp : α

/-! ## 1. Model (mirrors the Rust functions) -/

section Model

variable (α : Type) [Add α] [Sub α] [Mul α] [Div α] [Neg α] [LT α] [DecidableRel (α := α) (· < ·)]
  [BEq α] [Sqrt α] [Powf α] [Fabs α]
  [OfNat α 0] [OfNat α 1] [OfNat α 2] [OfNat α 3] [OfNat α 6] [OfNat α 12]

/-- `PerturbedLambdaParams::from_lambda(lj, ljs, crf, crfs, dlj, dljs, dcrf, dcrfs, n)`, with
`n as f64` and `(n − 1) as f64` given as `nF`, `nm1F`. -/
def lambdaParamsModel (lj ljs crf crfs dlj dljs dcrf dcrfs nF nm1F : α) : LambdaParams α where
  aLjN := Powf.powf (1 - lj) nF
  aCrfN := Powf.powf (1 - crf) nF
  bLjN := Powf.powf lj nF
  bCrfN := Powf.powf crf nF
  aLjN1 := Powf.powf (1 - lj) nm1F * dlj
  aCrfN1 := Powf.powf (1 - crf) nm1F * dcrf
  bLjN1 := Powf.powf lj nm1F * dlj
  bCrfN1 := Powf.powf crf nm1F * dcrf
  aLjs := (1 - ljs) * dljs
  aCrfs := (1 - crfs) * dcrfs
  bLjs := ljs * dljs
  bCrfs := crfs * dcrfs
  aLjs2 := (1 - ljs) * (1 - ljs)
  aCrfs2 := (1 - crfs) * (1 - crfs)
  bLjs2 := ljs * ljs
  bCrfs2 := crfs * crfs
  nExp := nF

/-- `C12/C6` for the soft-core length, `0` when `|C6| ≤ guard` (`1e-10` in the pair kernel). -/
def c126Model (guard c6 c12 : α) : α :=
  if guard < Fabs.fabs c6 then c12 / c6 else 0

/-- `crf/2 / sqrt(cut2soft³)` as the soft kernels form it: `crf_2cut3i · cut³ / sqrt(…)`. -/
def softCrf2Cut3iModel (p : CrfParams α) (cut2soft : α) : α :=
  p.crf2Cut3i * (p.cutoffSq * Sqrt.sqrt p.cutoffSq) / Sqrt.sqrt (cut2soft * cut2soft * cut2soft)

/-- `perturbed_lj_crf_pair(r2, a_c6, a_c12, b_c6, b_c12, a_q, b_q, alpha_lj, alpha_crf, lp, crf,
four_pi_eps_i, coulomb_scaling)`: `(force, e_lj, e_crf, de_lj, de_crf)`, with the distance guard
and the `|C6|` guard (both `1e-10`) as `guard`, `c6guard`. -/
def perturbedLjCrfPairModel (guard c6guard r2 aC6 aC12 bC6 bC12 aQ bQ αlj αc : α)
    (lp : LambdaParams α) (p : CrfParams α) (fpe cs : α) : α × α × α × α × α :=
  if r2 < guard then (0, 0, 0, 0, 0)
  else
    let r6 := r2 * r2 * r2
    let r4 := r2 * r2
    let aC126 := c126Model α c6guard aC6 aC12
    let bC126 := c126Model α c6guard bC6 bC12
    let aDist2soft := r2 + αc * lp.bCrfs2
    let aDistisoft := 1 / Sqrt.sqrt aDist2soft
    let aDist3isoft := aDistisoft / aDist2soft
    let aDist6soft := r6 + αlj * lp.bLjs2 * aC126
    let aDist6isoft := 1 / aDist6soft
    let bDist2soft := r2 + αc * lp.aCrfs2
    let bDistisoft := 1 / Sqrt.sqrt bDist2soft
    let bDist3isoft := bDistisoft / bDist2soft
    let bDist6soft := r6 + αlj * lp.aLjs2 * bC126
    let bDist6isoft := 1 / bDist6soft
    let cut2 := p.cutoffSq
    let aCut2soft := cut2 + αc * lp.bCrfs2
    let bCut2soft := cut2 + αc * lp.aCrfs2
    let aCut2soft3 := aCut2soft * aCut2soft * aCut2soft
    let bCut2soft3 := bCut2soft * bCut2soft * bCut2soft
    let cut3 := cut2 * Sqrt.sqrt cut2
    let aCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt aCut2soft3
    let bCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt bCut2soft3
    let aCrfCut3i := 2 * aCrf2cut3i
    let bCrfCut3i := 2 * bCrf2cut3i
    let aCrfPert := 3 * aCrf2cut3i / aCut2soft
    let bCrfPert := 3 * bCrf2cut3i / bCut2soft
    let force0 := (lp.aCrfN * aQ * (cs * aDist3isoft + aCrfCut3i) +
      lp.bCrfN * bQ * (cs * bDist3isoft + bCrfCut3i)) * fpe
    let force1 := force0 + (-6) * (lp.aLjN * aC6 * aDist6isoft * aDist6isoft +
      lp.bLjN * bC6 * bDist6isoft * bDist6isoft) * r4
    let force := force1 + 12 * (lp.aLjN * aC12 * aDist6isoft * aDist6isoft * aDist6isoft +
      lp.bLjN * bC12 * bDist6isoft * bDist6isoft * bDist6isoft) * r4
    let aELj := (aC12 * aDist6isoft - aC6) * aDist6isoft
    let bELj := (bC12 * bDist6isoft - bC6) * bDist6isoft
    let aECrf := aQ * (cs * aDistisoft - aCrf2cut3i * r2 - p.crfCut)
    let bECrf := bQ * (cs * bDistisoft - bCrf2cut3i * r2 - p.crfCut)
    let eLj := lp.aLjN * aELj + lp.bLjN * bELj
    let eCrf := (lp.aCrfN * aECrf + lp.bCrfN * bECrf) * fpe
    let deLj := (-2) * αlj *
        (lp.aLjN * lp.bLjs * aC126 * aDist6isoft * aDist6isoft * (2 * aC12 * aDist6isoft - aC6) -
          lp.bLjN * lp.aLjs * bC126 * bDist6isoft * bDist6isoft * (2 * bC12 * bDist6isoft - bC6)) +
      lp.nExp * (lp.bLjN1 * bELj - lp.aLjN1 * aELj)
    let deCrf := (-(lp.aCrfN * aQ * lp.bCrfs * (cs * aDist3isoft - aCrfPert * r2) -
          lp.bCrfN * bQ * lp.aCrfs * (cs * bDist3isoft - bCrfPert * r2))) * fpe * αc +
      lp.nExp * (lp.bCrfN1 * bECrf - lp.aCrfN1 * aECrf) * fpe
    (force, eLj, eCrf, deLj, deCrf)

/-- `rf_soft_interaction(r2, a_q, b_q, alpha_crf, lp, crf, four_pi_eps_i)`:
`(force, e_rf, de_rf)`. -/
def rfSoftModel (r2 aQ bQ αc : α) (lp : LambdaParams α) (p : CrfParams α) (fpe : α) :
    α × α × α :=
  let cut2 := p.cutoffSq
  let aCut2soft := cut2 + αc * lp.bCrfs2
  let bCut2soft := cut2 + αc * lp.aCrfs2
  let cut3 := cut2 * Sqrt.sqrt cut2
  let aCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt (aCut2soft * aCut2soft * aCut2soft)
  let bCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt (bCut2soft * bCut2soft * bCut2soft)
  let aCrfPert := 3 * aCrf2cut3i / aCut2soft
  let bCrfPert := 3 * bCrf2cut3i / bCut2soft
  let force := (lp.aCrfN * aQ * 2 * aCrf2cut3i + lp.bCrfN * bQ * 2 * bCrf2cut3i) * fpe
  let aERf := aQ * ((-aCrf2cut3i) * r2 - p.crfCut)
  let bERf := bQ * ((-bCrf2cut3i) * r2 - p.crfCut)
  let eRf := (lp.aCrfN * aERf + lp.bCrfN * bERf) * fpe
  let deRf := ((lp.aCrfN * aQ * lp.bCrfs * aCrfPert - lp.bCrfN * bQ * lp.aCrfs * bCrfPert) *
      r2 * αc + (lp.bCrfN1 * bERf - lp.aCrfN1 * aERf) * lp.nExp) * fpe
  (force, eRf, deRf)

/-- `x.powi(3)`: compiler-rt's `__powidf2` and LLVM's expansion both form `x · (x · x)`, which is
`x · x · x` bit for bit (rounded multiplication is commutative). -/
def powi3 (x : α) : α := x * x * x

/-- `perturbed_rf_excluded_pair(r2, a_q, b_q, alpha_crf, lp, crf, four_pi_eps_i)`:
`(force, e_crf, de_crf)`. -/
def perturbedRfExcludedModel (r2 aQ bQ αc : α) (lp : LambdaParams α) (p : CrfParams α)
    (fpe : α) : α × α × α :=
  let cut2 := p.cutoffSq
  let aCut2soft := cut2 + αc * lp.bCrfs2
  let bCut2soft := cut2 + αc * lp.aCrfs2
  let aCrf2cut3i := p.crf2Cut3i * Sqrt.sqrt (powi3 α (cut2 / aCut2soft))
  let bCrf2cut3i := p.crf2Cut3i * Sqrt.sqrt (powi3 α (cut2 / bCut2soft))
  let aCrfCut3i := 2 * aCrf2cut3i
  let bCrfCut3i := 2 * bCrf2cut3i
  let aCrfPert := 3 * aCrf2cut3i / aCut2soft
  let bCrfPert := 3 * bCrf2cut3i / bCut2soft
  let aECrf := aQ * ((-aCrf2cut3i) * r2 - p.crfCut)
  let bECrf := bQ * ((-bCrf2cut3i) * r2 - p.crfCut)
  let eCrf := (lp.aCrfN * aECrf + lp.bCrfN * bECrf) * fpe
  let force := (lp.aCrfN * aQ * aCrfCut3i + lp.bCrfN * bQ * bCrfCut3i) * fpe
  let deCrf := ((lp.aCrfN * aQ * lp.bCrfs * aCrfPert - lp.bCrfN * bQ * lp.aCrfs * bCrfPert) *
      r2 * αc + lp.nExp * (lp.bCrfN1 * bECrf - lp.aCrfN1 * aECrf)) * fpe
  (force, eCrf, deCrf)

/-- `perturbed_rf_self_term(qa2, qb2, lp, crf, four_pi_eps_i)`: `(e_rf, de_rf)`. -/
def perturbedRfSelfModel (qa2 qb2 : α) (lp : LambdaParams α) (p : CrfParams α) (fpe : α) :
    α × α :=
  ((-(lp.aCrfN * qa2) - lp.bCrfN * qb2 + qa2) * p.crfCut * fpe,
    lp.nExp * (lp.aCrfN1 * qa2 - lp.bCrfN1 * qb2) * p.crfCut * fpe)

/-- `lj_crf_soft_interaction_ext(r2, a_c6, a_c12, b_c6, b_c12, a_q, b_q, alpha_lj, alpha_crf,
crf, four_pi_eps_i, lam)`: `([A_e_lj, B_e_lj, A_e_crf, B_e_crf], [A_de_lj, …])`. -/
def ljCrfSoftExtModel (r2 aC6 aC12 bC6 bC12 aQ bQ αlj αc : α) (p : CrfParams α) (fpe lam : α) :
    (α × α × α × α) × (α × α × α × α) :=
  let aLjs2 := (1 - lam) * (1 - lam)
  let bLjs2 := lam * lam
  let aC126 := if aC6 != 0 then aC12 / aC6 else 0
  let bC126 := if bC6 != 0 then bC12 / bC6 else 0
  let aDist2soft := r2 + αc * bLjs2
  let bDist2soft := r2 + αc * aLjs2
  let aDistisoft := 1 / Sqrt.sqrt aDist2soft
  let bDistisoft := 1 / Sqrt.sqrt bDist2soft
  let aDist3isoft := aDistisoft / aDist2soft
  let bDist3isoft := bDistisoft / bDist2soft
  let dist6 := r2 * r2 * r2
  let aDist6isoft := 1 / (dist6 + αlj * bLjs2 * aC126)
  let bDist6isoft := 1 / (dist6 + αlj * aLjs2 * bC126)
  let cut2 := p.cutoffSq
  let aCut2soft := cut2 + αc * bLjs2
  let bCut2soft := cut2 + αc * aLjs2
  let cut3 := cut2 * Sqrt.sqrt cut2
  let aCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt (aCut2soft * aCut2soft * aCut2soft)
  let bCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt (bCut2soft * bCut2soft * bCut2soft)
  let aCrfPert := 3 * aCrf2cut3i / aCut2soft
  let bCrfPert := 3 * bCrf2cut3i / bCut2soft
  (((aC12 * aDist6isoft - aC6) * aDist6isoft,
    (bC12 * bDist6isoft - bC6) * bDist6isoft,
    aQ * (aDistisoft - aCrf2cut3i * r2 - p.crfCut) * fpe,
    bQ * (bDistisoft - bCrf2cut3i * r2 - p.crfCut) * fpe),
   ((-2) * αlj * lam * aC126 * aDist6isoft * aDist6isoft * (2 * aC12 * aDist6isoft - aC6),
    (-2) * αlj * (1 - lam) * bC126 * bDist6isoft * bDist6isoft * (2 * bC12 * bDist6isoft - bC6),
    (-aQ) * lam * (aDist3isoft - aCrfPert * r2) * fpe * αc,
    (-bQ) * (1 - lam) * (bDist3isoft - bCrfPert * r2) * fpe * αc))

/-- `rf_soft_interaction_ext(r2, a_q, b_q, alpha_crf, crf, four_pi_eps_i, lam)`:
`([A_e_rf, B_e_rf], [A_de_rf, B_de_rf])`. -/
def rfSoftExtModel (r2 aQ bQ αc : α) (p : CrfParams α) (fpe lam : α) : (α × α) × (α × α) :=
  let cut2 := p.cutoffSq
  let aCut2soft := cut2 + αc * lam * lam
  let bCut2soft := cut2 + αc * (1 - lam) * (1 - lam)
  let cut3 := cut2 * Sqrt.sqrt cut2
  let aCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt (aCut2soft * aCut2soft * aCut2soft)
  let bCrf2cut3i := p.crf2Cut3i * cut3 / Sqrt.sqrt (bCut2soft * bCut2soft * bCut2soft)
  let aCrfPert := 3 * aCrf2cut3i / aCut2soft
  let bCrfPert := 3 * bCrf2cut3i / bCut2soft
  ((aQ * ((-aCrf2cut3i) * r2 - p.crfCut) * fpe, bQ * ((-bCrf2cut3i) * r2 - p.crfCut) * fpe),
   (aQ * lam * aCrfPert * r2 * fpe * αc, bQ * (1 - lam) * bCrfPert * r2 * fpe * αc))

end Model

/-! ## 2. The `Float` instantiation (what `lake exe spec-vectors` runs) -/

/-- The `1e-10` of the pair kernel's guards (`r2 < 1e-10`, `|C6| > 1e-10`), as `f64` bits; the
same double as `ZERO_DISTANCE2`. -/
def guardFloat : Float := Float.ofBits 0x3ddb7cdfd9d7bdbb

def lambdaParamsFloat (lj ljs crf crfs dlj dljs dcrf dcrfs nF nm1F : Float) : LambdaParams Float :=
  lambdaParamsModel Float lj ljs crf crfs dlj dljs dcrf dcrfs nF nm1F

def perturbedLjCrfPairFloat (r2 aC6 aC12 bC6 bC12 aQ bQ αlj αc : Float) (lp : LambdaParams Float)
    (p : CrfParams Float) (fpe cs : Float) : Float × Float × Float × Float × Float :=
  perturbedLjCrfPairModel Float guardFloat guardFloat r2 aC6 aC12 bC6 bC12 aQ bQ αlj αc lp p fpe cs

def rfSoftFloat (r2 aQ bQ αc : Float) (lp : LambdaParams Float) (p : CrfParams Float)
    (fpe : Float) : Float × Float × Float :=
  rfSoftModel Float r2 aQ bQ αc lp p fpe

def perturbedRfExcludedFloat (r2 aQ bQ αc : Float) (lp : LambdaParams Float)
    (p : CrfParams Float) (fpe : Float) : Float × Float × Float :=
  perturbedRfExcludedModel Float r2 aQ bQ αc lp p fpe

def perturbedRfSelfFloat (qa2 qb2 : Float) (lp : LambdaParams Float) (p : CrfParams Float)
    (fpe : Float) : Float × Float :=
  perturbedRfSelfModel Float qa2 qb2 lp p fpe

def ljCrfSoftExtFloat (r2 aC6 aC12 bC6 bC12 aQ bQ αlj αc : Float) (p : CrfParams Float)
    (fpe lam : Float) : (Float × Float × Float × Float) × (Float × Float × Float × Float) :=
  ljCrfSoftExtModel Float r2 aC6 aC12 bC6 bC12 aQ bQ αlj αc p fpe lam

def rfSoftExtFloat (r2 aQ bQ αc : Float) (p : CrfParams Float) (fpe lam : Float) :
    (Float × Float) × (Float × Float) :=
  rfSoftExtModel Float r2 aQ bQ αc p fpe lam

/-! ## 3. The ℝ instantiation equals the theory -/

/-- The λ table over ℝ, from λ values and NLAM `n ≥ 1`: `powf` is `Real.rpow`, which at the
natural exponents `n`, `n − 1` is the power. -/
theorem lambdaParamsModel_real (lj ljs crf crfs dlj dljs dcrf dcrfs : ℝ) (n : ℕ) :
    lambdaParamsModel ℝ lj ljs crf crfs dlj dljs dcrf dcrfs n ((n - 1 : ℕ) : ℝ) =
      { aLjN := (1 - lj) ^ n, bLjN := lj ^ n, aCrfN := (1 - crf) ^ n, bCrfN := crf ^ n,
        aLjN1 := (1 - lj) ^ (n - 1) * dlj, bLjN1 := lj ^ (n - 1) * dlj,
        aCrfN1 := (1 - crf) ^ (n - 1) * dcrf, bCrfN1 := crf ^ (n - 1) * dcrf,
        aLjs := (1 - ljs) * dljs, bLjs := ljs * dljs,
        aLjs2 := (1 - ljs) * (1 - ljs), bLjs2 := ljs * ljs,
        aCrfs := (1 - crfs) * dcrfs, bCrfs := crfs * dcrfs,
        aCrfs2 := (1 - crfs) * (1 - crfs), bCrfs2 := crfs * crfs, nExp := n } := by
  simp only [lambdaParamsModel, Powf.powf, Real.rpow_eq_pow, Real.rpow_natCast]

/-- The table of individual λ `(lj, ljs, crf, crfs)` with derivatives `(dlj, dljs, dcrf, dcrfs)`
and NLAM `n`, over ℝ. -/
noncomputable abbrev tableR (n : ℕ) (lj ljs crf crfs dlj dljs dcrf dcrfs : ℝ) : LambdaParams ℝ :=
  lambdaParamsModel ℝ lj ljs crf crfs dlj dljs dcrf dcrfs n ((n - 1 : ℕ) : ℝ)

/-- The Lennard-Jones state the pair kernel sees: `C6`, `C12` and its guarded `C12/C6`. -/
noncomputable abbrev ljStateOf (c6guard c6 c12 : ℝ) : LjState := ⟨c6, c12, c126Model ℝ c6guard c6 c12⟩

/-- The softened reaction-field coefficient over ℝ: `crf_2cut3i · cut³` is `crf/2` for `R > 0`. -/
theorem crf2Cut3i_cut3_real (crf : ℝ) {R : ℝ} (hR : 0 < R) :
    (crfDerivedModel ℝ crf R).crf2Cut3i *
        ((crfDerivedModel ℝ crf R).cutoffSq * Real.sqrt (crfDerivedModel ℝ crf R).cutoffSq) =
      crf / 2 := by
  simp only [crfDerivedModel, Real.sqrt_mul_self hR.le]
  have := hR.ne'
  field_simp

/-- **The pair kernel over ℝ is the theory.** For `r > 0` outside the guard, `R > 0` and NLAM
`n`, the model returns the soft-core force, energies and `dV/dλ` of
`NewenTheory/Interactions/SoftCore.lean`, with the individual λ `(lj, ljs, crf, crfs)` and their
derivatives as the table kernel receives them. -/
theorem perturbedLjCrfPairModel_real_eq (n : ℕ) {guard c6guard r R : ℝ} (hR : 0 < R)
    (hg : guard ≤ r ^ 2) (aC6 aC12 bC6 bC12 qA qB αlj αc fpe cs crf : ℝ)
    (lj ljs lc lcs dlj dljs dlc dlcs : ℝ) :
    perturbedLjCrfPairModel ℝ guard c6guard (r ^ 2) aC6 aC12 bC6 bC12 qA qB αlj αc
        (tableR n lj ljs lc lcs dlj dljs dlc dlcs) (crfDerivedModel ℝ crf R) fpe cs =
      (pertLjForce n αlj (ljStateOf c6guard aC6 aC12) (ljStateOf c6guard bC6 bC12) r lj ljs +
          pertCrfForce n αc cs fpe qA qB crf R r lc lcs,
        pertLjEnergy n αlj (ljStateOf c6guard aC6 aC12) (ljStateOf c6guard bC6 bC12) r lj ljs,
        pertCrfEnergy n αc cs fpe qA qB crf R r lc lcs,
        pertLjDvdl n αlj (ljStateOf c6guard aC6 aC12) (ljStateOf c6guard bC6 bC12) r lj ljs
          dlj dljs,
        pertCrfDvdl n αc cs fpe qA qB crf R r lc lcs dlc dlcs) := by
  have hc := crf2Cut3i_cut3_real crf hR
  simp only [perturbedLjCrfPairModel, tableR, lambdaParamsModel_real, not_lt.mpr hg, ite_false,
    Sqrt.sqrt, hc]
  simp only [crfDerivedModel, pertLjForce, pertCrfForce, pertLjEnergy, pertCrfEnergy, pertLjDvdl,
    pertCrfDvdl, scLjEnergy, scLjDs, scLjForce, scCrfEnergy, scCrfDs, scCrfForce, Prod.mk.injEq]
  simp only [div_eq_mul_inv, ← inv_pow]
  ring_nf
  simp only [and_self]

/-- **The pair kernel's `dV/dλ` is the derivative of its energy.** Over ℝ, with the λ table
built by `lambdaParamsModel` from individual λ `Λ`, `Λₛ`, `Λ_c`, `Λ_cs` that are differentiable
in λ (derivatives `dΛ`, … at `lam`, the inputs the table receives), the `de_lj + de_crf` the
pair kernel returns at `lam` is the derivative at `lam` of the `e_lj + e_crf` it returns, wherever
the softened distances are positive. -/
theorem perturbedLjCrfPairModel_dvdl_is_deriv (n : ℕ) {guard c6guard r R : ℝ} (hR : 0 < R)
    (hg : guard ≤ r ^ 2) (aC6 aC12 bC6 bC12 qA qB αlj αc fpe cs crf : ℝ)
    {Λ Λs Λc Λcs : ℝ → ℝ} {dΛ dΛs dΛc dΛcs lam : ℝ}
    (hΛ : HasDerivAt Λ dΛ lam) (hΛs : HasDerivAt Λs dΛs lam)
    (hΛc : HasDerivAt Λc dΛc lam) (hΛcs : HasDerivAt Λcs dΛcs lam)
    (hA : r ^ 6 + αlj * Λs lam ^ 2 * c126Model ℝ c6guard aC6 aC12 ≠ 0)
    (hB : r ^ 6 + αlj * (1 - Λs lam) ^ 2 * c126Model ℝ c6guard bC6 bC12 ≠ 0)
    (hXA : 0 < r ^ 2 + αc * Λcs lam ^ 2) (hYA : 0 < R ^ 2 + αc * Λcs lam ^ 2)
    (hXB : 0 < r ^ 2 + αc * (1 - Λcs lam) ^ 2) (hYB : 0 < R ^ 2 + αc * (1 - Λcs lam) ^ 2) :
    HasDerivAt
      (fun l =>
        let o := perturbedLjCrfPairModel ℝ guard c6guard (r ^ 2) aC6 aC12 bC6 bC12 qA qB αlj αc
          (tableR n (Λ l) (Λs l) (Λc l) (Λcs l) dΛ dΛs dΛc dΛcs) (crfDerivedModel ℝ crf R) fpe cs
        o.2.1 + o.2.2.1)
      (let o := perturbedLjCrfPairModel ℝ guard c6guard (r ^ 2) aC6 aC12 bC6 bC12 qA qB αlj αc
          (tableR n (Λ lam) (Λs lam) (Λc lam) (Λcs lam) dΛ dΛs dΛc dΛcs)
          (crfDerivedModel ℝ crf R) fpe cs
       o.2.2.2.1 + o.2.2.2.2)
      lam := by
  simp only [perturbedLjCrfPairModel_real_eq n hR hg]
  exact softcore_dvdl_is_deriv n αlj αc cs fpe _ _ qA qB crf R r hΛ hΛs hΛc hΛcs hA hB hXA hYA
    hXB hYB

/-- The excluded-pair kernel `rf_soft_interaction` over ℝ is the theory's excluded soft-core
reaction field. -/
theorem rfSoftModel_real_eq (n : ℕ) {R : ℝ} (hR : 0 < R) (r qA qB αc fpe crf : ℝ)
    (lj ljs lc lcs dlj dljs dlc dlcs : ℝ) :
    rfSoftModel ℝ (r ^ 2) qA qB αc (tableR n lj ljs lc lcs dlj dljs dlc dlcs)
        (crfDerivedModel ℝ crf R) fpe =
      (pertRfExclForce n αc fpe qA qB crf R lc lcs, pertRfExclEnergy n αc fpe qA qB crf R r lc lcs,
        pertRfExclDvdl n αc fpe qA qB crf R r lc lcs dlc dlcs) := by
  have hc := crf2Cut3i_cut3_real crf hR
  simp only [rfSoftModel, tableR, lambdaParamsModel_real, Sqrt.sqrt, hc]
  simp only [crfDerivedModel, pertRfExclForce, pertRfExclEnergy, pertRfExclDvdl, scRfExclEnergy,
    scRfExclDs, Prod.mk.injEq]
  simp only [div_eq_mul_inv]
  ring_nf
  simp only [and_self]

/-- The excluded-pair correction's coefficient over ℝ: `crf_2cut3i · √((R²/Y)³)` is
`crf/2 / √(Y³)` for `R > 0`, `Y ≥ 0`. -/
theorem crf2Cut3i_ratio_real (crf : ℝ) {R Y : ℝ} (hR : 0 < R) (hY : 0 ≤ Y) :
    (crfDerivedModel ℝ crf R).crf2Cut3i *
        Real.sqrt ((crfDerivedModel ℝ crf R).cutoffSq / Y * ((crfDerivedModel ℝ crf R).cutoffSq / Y) *
          ((crfDerivedModel ℝ crf R).cutoffSq / Y)) =
      crf / 2 / Real.sqrt (Y ^ 3) := by
  have h : R * R / Y * (R * R / Y) * (R * R / Y) = (R ^ 3) ^ 2 / Y ^ 3 := by ring
  simp only [crfDerivedModel, h, Real.sqrt_div' _ (pow_nonneg hY 3),
    Real.sqrt_sq (pow_nonneg hR.le 3)]
  have := hR.ne'
  field_simp

/-- `perturbed_rf_excluded_pair` over ℝ is the same theory as `rf_soft_interaction` (it differs
only in rounding), where both softened cutoffs are positive. -/
theorem perturbedRfExcludedModel_real_eq (n : ℕ) {R : ℝ} (hR : 0 < R) (r qA qB αc fpe crf : ℝ)
    (lj ljs lc lcs dlj dljs dlc dlcs : ℝ) (hYA : 0 < R ^ 2 + αc * lcs ^ 2)
    (hYB : 0 < R ^ 2 + αc * (1 - lcs) ^ 2) :
    perturbedRfExcludedModel ℝ (r ^ 2) qA qB αc (tableR n lj ljs lc lcs dlj dljs dlc dlcs)
        (crfDerivedModel ℝ crf R) fpe =
      (pertRfExclForce n αc fpe qA qB crf R lc lcs, pertRfExclEnergy n αc fpe qA qB crf R r lc lcs,
        pertRfExclDvdl n αc fpe qA qB crf R r lc lcs dlc dlcs) := by
  have hA := crf2Cut3i_ratio_real (Y := R * R + αc * (lcs * lcs)) crf hR (by nlinarith)
  have hB := crf2Cut3i_ratio_real (Y := R * R + αc * ((1 - lcs) * (1 - lcs))) crf hR
    (by nlinarith)
  simp only [perturbedRfExcludedModel, powi3, tableR, lambdaParamsModel_real, Sqrt.sqrt]
  simp only [crfDerivedModel] at hA hB ⊢
  rw [hA, hB]
  simp only [pertRfExclForce, pertRfExclEnergy, pertRfExclDvdl, scRfExclEnergy, scRfExclDs,
    Prod.mk.injEq]
  simp only [div_eq_mul_inv]
  ring_nf
  simp only [and_self]

/-- The self-term kernel over ℝ: `(rfSelfChange, its λ derivative)` with
`c = crf_cut · fpe`. -/
theorem perturbedRfSelfModel_real_eq (n : ℕ) (qa2 qb2 fpe : ℝ) (p : CrfParams ℝ)
    (lj ljs lc lcs dlj dljs dlc dlcs : ℝ) :
    perturbedRfSelfModel ℝ qa2 qb2 (tableR n lj ljs lc lcs dlj dljs dlc dlcs) p fpe =
      (rfSelfChange n qa2 qb2 (p.crfCut * fpe) lc,
        n * ((1 - lc) ^ (n - 1) * qa2 - lc ^ (n - 1) * qb2) * (p.crfCut * fpe) * dlc) := by
  simp only [perturbedRfSelfModel, tableR, lambdaParamsModel_real, rfSelfChange, Prod.mk.injEq]
  constructor <;> ring

/-- The self-term kernel's `de_rf` is the derivative of its `e_rf`. -/
theorem perturbedRfSelfModel_dvdl_is_deriv (n : ℕ) (qa2 qb2 fpe : ℝ) (p : CrfParams ℝ)
    {Λ Λs Λc Λcs : ℝ → ℝ} {dΛ dΛs dΛc dΛcs lam : ℝ} (hΛc : HasDerivAt Λc dΛc lam) :
    HasDerivAt
      (fun l => (perturbedRfSelfModel ℝ qa2 qb2
        (tableR n (Λ l) (Λs l) (Λc l) (Λcs l) dΛ dΛs dΛc dΛcs) p fpe).1)
      (perturbedRfSelfModel ℝ qa2 qb2
        (tableR n (Λ lam) (Λs lam) (Λc lam) (Λcs lam) dΛ dΛs dΛc dΛcs) p fpe).2 lam := by
  simp only [perturbedRfSelfModel_real_eq]
  exact rf_self_dvdl_is_deriv n qa2 qb2 (p.crfCut * fpe) hΛc

/-- The Lennard-Jones state the ext kernel sees: its `C12/C6` is guarded by `C6 ≠ 0`. -/
noncomputable abbrev extStateOf (c6 c12 : ℝ) : LjState := ⟨c6, c12, if c6 != 0 then c12 / c6 else 0⟩

/-- `lj_crf_soft_interaction_ext` over ℝ: the unweighted soft-core energies of the two states at
softness `λ²` and `(1 − λ)²`, and `A_de = ∂E_A/∂λ`, `B_de = −∂E_B/∂λ` written as the chain rule
`∂E/∂s · ds/dλ`. -/
theorem ljCrfSoftExtModel_real_eq {R : ℝ} (hR : 0 < R) (r aC6 aC12 bC6 bC12 qA qB αlj αc fpe crf
    lam : ℝ) :
    ljCrfSoftExtModel ℝ (r ^ 2) aC6 aC12 bC6 bC12 qA qB αlj αc (crfDerivedModel ℝ crf R) fpe lam =
      ((scLjEnergy αlj (extStateOf aC6 aC12) (lam ^ 2) r,
          scLjEnergy αlj (extStateOf bC6 bC12) ((1 - lam) ^ 2) r,
          scCrfEnergy αc 1 qA crf R (lam ^ 2) r * fpe,
          scCrfEnergy αc 1 qB crf R ((1 - lam) ^ 2) r * fpe),
        (scLjDs αlj (extStateOf aC6 aC12) (lam ^ 2) r * (2 * lam),
          -(scLjDs αlj (extStateOf bC6 bC12) ((1 - lam) ^ 2) r * (-(2 * (1 - lam)))),
          scCrfDs αc 1 qA crf R (lam ^ 2) r * (2 * lam) * fpe,
          -(scCrfDs αc 1 qB crf R ((1 - lam) ^ 2) r * (-(2 * (1 - lam)))) * fpe)) := by
  have hc := crf2Cut3i_cut3_real crf hR
  simp only [ljCrfSoftExtModel, Sqrt.sqrt, hc]
  simp only [crfDerivedModel, scLjEnergy, scLjDs, scCrfEnergy, scCrfDs, Prod.mk.injEq]
  simp only [div_eq_mul_inv, ← inv_pow]
  ring_nf
  simp only [and_self]

/-- **The ext kernel's derivatives are the derivatives of its energies**: `A_de = ∂E_A/∂λ` and
`B_de = −∂E_B/∂λ` (gromosXX's sign convention for B), for the LJ and the CRF columns, wherever
the softened distances are positive. -/
theorem ljCrfSoftExtModel_de_is_deriv {R : ℝ} (hR : 0 < R)
    (r aC6 aC12 bC6 bC12 qA qB αlj αc fpe crf lam : ℝ)
    (hA : r ^ 6 + αlj * lam ^ 2 * (extStateOf aC6 aC12).c126 ≠ 0)
    (hB : r ^ 6 + αlj * (1 - lam) ^ 2 * (extStateOf bC6 bC12).c126 ≠ 0)
    (hXA : 0 < r ^ 2 + αc * lam ^ 2) (hYA : 0 < R ^ 2 + αc * lam ^ 2)
    (hXB : 0 < r ^ 2 + αc * (1 - lam) ^ 2) (hYB : 0 < R ^ 2 + αc * (1 - lam) ^ 2) :
    let ext := fun l =>
      ljCrfSoftExtModel ℝ (r ^ 2) aC6 aC12 bC6 bC12 qA qB αlj αc (crfDerivedModel ℝ crf R) fpe l
    HasDerivAt (fun l => (ext l).1.1) (ext lam).2.1 lam ∧
      HasDerivAt (fun l => (ext l).1.2.1) (-(ext lam).2.2.1) lam ∧
      HasDerivAt (fun l => (ext l).1.2.2.1) (ext lam).2.2.2.1 lam ∧
      HasDerivAt (fun l => (ext l).1.2.2.2) (-(ext lam).2.2.2.2) lam := by
  intro ext
  simp only [ext, ljCrfSoftExtModel_real_eq hR, neg_neg]
  have sa := hasDerivAt_soft_a (hasDerivAt_id lam)
  have sb := hasDerivAt_soft_b (hasDerivAt_id lam)
  simp only [id] at sa sb
  refine ⟨?_, ?_, ?_, ?_⟩
  · have h := (scLj_hasDerivAt_s αlj _ r hA).comp (h := fun l => l ^ 2) lam sa
    refine h.congr_deriv ?_
    ring
  · have h := (scLj_hasDerivAt_s αlj _ r hB).comp lam sb
    refine h.congr_deriv ?_
    ring
  · have h := ((scCrf_hasDerivAt_s αc 1 qA crf R r hXA hYA).comp (h := fun l => l ^ 2) lam sa).mul_const fpe
    refine h.congr_deriv ?_
    ring
  · have h := ((scCrf_hasDerivAt_s αc 1 qB crf R r hXB hYB).comp lam sb).mul_const fpe
    refine h.congr_deriv ?_
    ring

/-- `rf_soft_interaction_ext` over ℝ: the unweighted excluded-pair energies and
`A_de = ∂E_A/∂λ`, `B_de = −∂E_B/∂λ`. -/
theorem rfSoftExtModel_real_eq {R : ℝ} (hR : 0 < R) (r qA qB αc fpe crf lam : ℝ) :
    rfSoftExtModel ℝ (r ^ 2) qA qB αc (crfDerivedModel ℝ crf R) fpe lam =
      ((scRfExclEnergy αc qA crf R (lam ^ 2) r * fpe,
          scRfExclEnergy αc qB crf R ((1 - lam) ^ 2) r * fpe),
        (scRfExclDs αc qA crf R (lam ^ 2) r * (2 * lam) * fpe,
          -(scRfExclDs αc qB crf R ((1 - lam) ^ 2) r * (-(2 * (1 - lam)))) * fpe)) := by
  have hc := crf2Cut3i_cut3_real crf hR
  simp only [rfSoftExtModel, Sqrt.sqrt, hc]
  simp only [crfDerivedModel, scRfExclEnergy, scRfExclDs, Prod.mk.injEq]
  simp only [div_eq_mul_inv]
  ring_nf
  simp only [and_self]

/-- The `rf_soft_interaction_ext` derivatives are the derivatives of its energies. -/
theorem rfSoftExtModel_de_is_deriv {R : ℝ} (hR : 0 < R) (r qA qB αc fpe crf lam : ℝ)
    (hYA : 0 < R ^ 2 + αc * lam ^ 2) (hYB : 0 < R ^ 2 + αc * (1 - lam) ^ 2) :
    let ext := fun l => rfSoftExtModel ℝ (r ^ 2) qA qB αc (crfDerivedModel ℝ crf R) fpe l
    HasDerivAt (fun l => (ext l).1.1) (ext lam).2.1 lam ∧
      HasDerivAt (fun l => (ext l).1.2) (-(ext lam).2.2) lam := by
  intro ext
  simp only [ext, rfSoftExtModel_real_eq hR]
  have sa := hasDerivAt_soft_a (hasDerivAt_id lam)
  have sb := hasDerivAt_soft_b (hasDerivAt_id lam)
  simp only [id] at sa sb
  refine ⟨?_, ?_⟩
  · have h := ((scRfExcl_hasDerivAt_s αc qA crf R r hYA).comp (h := fun l => l ^ 2) lam sa).mul_const fpe
    refine h.congr_deriv ?_
    ring
  · have h := ((scRfExcl_hasDerivAt_s αc qB crf R r hYB).comp lam sb).mul_const fpe
    refine h.congr_deriv ?_
    ring

end NewenSpec.SoftCore
