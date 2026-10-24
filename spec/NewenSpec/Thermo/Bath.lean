import NewenTheory.Dynamics.DegreesOfFreedom

/-!
# Model of `crates/newen-logic/src/thermo/bath.rs`

Rust: `crates/newen-logic/src/thermo/bath.rs`: `baths_of`, `groups_starting_in`, `bath_half_dof`,
`bath_half_dofs`, `ndfmin_spread`, `empty_baths` (spec keys `thermo.bath.*` in `spec/SPEC.toml`).
Theory: `NewenTheory/Dynamics/DegreesOfFreedom.lean`, section "Temperature baths" (`BathSystem`,
`comHalf`, `irHalf`, `spreadNumerators`, `spreadDenominator`, `EmptyBath`).

This file has two parts:

1. the *models*, written as the Rust is: over the DOFSET ranges of a `MULTIBATH` block, range by
   range, with the first atom of each range carried along;
2. the *bridges*: when the ranges partition the atoms (`RangesWellFormed`, what `newen-run`'s plan
   validation checks), the range-by-range count is the theory's atom-by-atom count, and the
   spread is the theory's spread.

The Aeneas refinements (Rust = model) are in `NewenSpec/Refinement.lean`.

Only core Lean is used, no Mathlib.
-/

namespace NewenSpec.Thermo.Bath

open NewenTheory.Dynamics.DegreesOfFreedom

/-! ## 1. Models (mirror the Rust functions) -/

/-- A DOFSET range (`DofRange`): the atoms up to and including `lastAtom` not in an earlier
range, with their COM and IR bath. -/
structure Range where
  lastAtom : Nat
  comBath : Nat
  irBath : Nat

/-- `baths_of`: the baths of the first range that reaches `atom`. -/
def bathsOfModel : List Range → Nat → Option (Nat × Nat)
  | [], _ => none
  | r :: rs, atom => if atom ≤ r.lastAtom then some (r.comBath, r.irBath) else bathsOfModel rs atom

/-- `groups_starting_in`: how many group firsts lie in `first ..= last`. -/
def groupsStartingInModel : List Nat → Nat → Nat → Nat
  | [], _, _ => 0
  | g :: gs, first, last =>
    (if first ≤ g ∧ g ≤ last then 1 else 0) + groupsStartingInModel gs first last

/-- The first loop of `bath_half_dof`: the ranges' contributions `(com, ir)` to bath `bath`,
the first range starting at atom `start`. A range has `lastAtom + 1 − start` atoms (`Nat`
subtraction: 0 when it is empty or backwards). -/
def rangesHalfModel (bath : Nat) (groupFirsts : List Nat) : Nat → List Range → Int × Int
  | _, [] => (0, 0)
  | start, r :: rs =>
    let atoms : Int := ((r.lastAtom + 1 - start : Nat) : Int)
    let groups : Int := (groupsStartingInModel groupFirsts start r.lastAtom : Int)
    let rest := rangesHalfModel bath groupFirsts (r.lastAtom + 1) rs
    ((if r.comBath = bath then 6 * groups else 0) + rest.1,
     (if r.irBath = bath then 6 * atoms - 6 * groups else 0) + rest.2)

/-- 1 if `atom` is reached by a range whose IR bath is `bath`. -/
def endInModel (bath : Nat) (ranges : List Range) (atom : Nat) : Int :=
  match bathsOfModel ranges atom with
  | some (_, ir) => if ir = bath then 1 else 0
  | none => 0

/-- The second loop of `bath_half_dof`: constraint ends in the IR bath `bath`. -/
def endsModel (bath : Nat) (ranges : List Range) : List (Nat × Nat) → Int
  | [] => 0
  | c :: cs => endInModel bath ranges c.1 + endInModel bath ranges c.2 + endsModel bath ranges cs

/-- `bath_half_dof`: `(com, ir)` half degrees of freedom of bath `bath`. -/
def bathHalfDofModel (bath : Nat) (ranges : List Range) (groupFirsts : List Nat)
    (constraints : List (Nat × Nat)) (rottrans : Bool) : Int × Int :=
  let r := rangesHalfModel bath groupFirsts 0 ranges
  let com := r.1
  let ir := r.2 - endsModel bath ranges constraints
  if rottrans then
    match bathsOfModel ranges 0 with
    | some (com0, ir0) => (if com0 = bath then com - 6 else com, if ir0 = bath then ir - 6 else ir)
    | none => (com, ir)
  else (com, ir)

/-- `bath_half_dofs`: `bathHalfDofModel` for the baths `0 .. n − 1`. -/
def bathHalfDofsModel (nBaths : Nat) (ranges : List Range) (groupFirsts : List Nat)
    (constraints : List (Nat × Nat)) (rottrans : Bool) : List (Int × Int) :=
  (List.range nBaths).map fun b => bathHalfDofModel b ranges groupFirsts constraints rottrans

/-- `ndfmin_spread`: `(numerators, denominator)`. -/
def ndfminSpreadModel (halfDofs : List Int) (ndfmin : Int) : List Int × Int :=
  let total := halfDofs.sum
  let fd : Int × Int := if 0 < total then (total - 2 * ndfmin, 2 * total) else (1, 2)
  (halfDofs.map (fun h => h * fd.1), fd.2)

/-- `empty_baths`: `h ≤ 0` for every bath's half count before the spread. -/
def emptyBathsModel (halves : List Int) : List Bool := halves.map fun h => decide (h ≤ 0)

/-! ## 2. Bridges: the models are the theory's counts -/

/-- The ranges, the first starting at `start`, cover the atoms `start .. n − 1` one after the
other: each range is non-empty (`start ≤ lastAtom`), the next starts after it, and the last ends
at atom `n − 1`. From `start = 0` this is "`lastAtom` strictly increasing, the last one `n − 1`". -/
def Covers : Nat → List Range → Nat → Prop
  | start, [], n => start = n
  | start, r :: rs, n => start ≤ r.lastAtom ∧ Covers (r.lastAtom + 1) rs n

/-- What `newen-run`'s plan validation guarantees for `n` atoms and `nBaths` baths: the DOFSET
ranges cover the atoms in order (`lastAtom` strictly increasing, the last one `n − 1`), every
bath index is a bath, and the temperature groups' first atoms are strictly increasing from 0.
The bridges below use only the first part. -/
def RangesWellFormed (nBaths n : Nat) (ranges : List Range) (groupFirsts : List Nat) : Prop :=
  Covers 0 ranges n ∧ (∀ r, r ∈ ranges → r.comBath < nBaths ∧ r.irBath < nBaths)
    ∧ groupFirsts.Pairwise (· < ·) ∧ groupFirsts.head? = some 0

/-- The COM bath of `atom` (0 for an atom no range reaches). -/
def comOf (ranges : List Range) (atom : Nat) : Nat :=
  match bathsOfModel ranges atom with
  | some (c, _) => c
  | none => 0

/-- The IR bath of `atom` (0 for an atom no range reaches). -/
def irOf (ranges : List Range) (atom : Nat) : Nat :=
  match bathsOfModel ranges atom with
  | some (_, i) => i
  | none => 0

/-- The theory's system for `n` atoms with the baths the ranges give. -/
def systemOf (n : Nat) (ranges : List Range) (groupFirsts : List Nat)
    (constraints : List (Nat × Nat)) (rottrans : Bool) : BathSystem where
  atoms := n
  comBath := comOf ranges
  irBath := irOf ranges
  groupFirsts := groupFirsts
  constraints := constraints
  rottrans := rottrans

theorem covers_le {start : Nat} {rs : List Range} {n : Nat} (h : Covers start rs n) :
    start ≤ n := by
  induction rs generalizing start with
  | nil => simp only [Covers] at h; omega
  | cons r rs ih => simp only [Covers] at h; have := ih h.2; omega

/-- From `start` on, an atom is reached by some range exactly when it is one of the `n` atoms. -/
theorem bathsOf_some_iff {start : Nat} {rs : List Range} {n : Nat} (h : Covers start rs n)
    (a : Nat) (ha : start ≤ a) : (bathsOfModel rs a).isSome ↔ a < n := by
  induction rs generalizing start with
  | nil => simp only [Covers] at h; simp [bathsOfModel]; omega
  | cons r rs ih =>
    simp only [Covers] at h
    have hle := covers_le h.2
    by_cases har : a ≤ r.lastAtom
    · simp [bathsOfModel, har]; omega
    · simp only [bathsOfModel, har, ite_false]
      exact ih h.2 (by omega)

theorem bathsOf_eq {rs : List Range} (a : Nat) (h : (bathsOfModel rs a).isSome) :
    bathsOfModel rs a = some (comOf rs a, irOf rs a) := by
  unfold comOf irOf
  cases e : bathsOfModel rs a with
  | none => simp [e] at h
  | some p => rfl

theorem sumOver_zero {α : Type} (xs : List α) (f : α → Int) (h : ∀ x, x ∈ xs → f x = 0) :
    sumOver xs f = 0 :=
  (sumOver_congr xs f (fun _ => 0) h).trans ((sumOver_const xs 0).trans (by simp))

theorem sumBelow_mul (n : Nat) (c : Int) (f : Nat → Int) :
    sumBelow n (fun a => c * f a) = c * sumBelow n f := by
  induction n with
  | zero => simp [sumBelow]
  | succ n ih => simp only [sumBelow, ih, Int.mul_add]

theorem sumOver_mul {α : Type} (xs : List α) (c : Int) (f : α → Int) :
    sumOver xs (fun x => c * f x) = c * sumOver xs f := by
  induction xs with
  | nil => simp [sumOver]
  | cons x xs ih => simp only [sumOver, ih, Int.mul_add]

/-- `groups_starting_in` counts the group firsts in `first ..= last`. -/
theorem groupsStartingIn_eq (gs : List Nat) (first last : Nat) :
    (groupsStartingInModel gs first last : Int)
      = sumOver gs (fun g => if first ≤ g ∧ g ≤ last then 1 else 0) := by
  induction gs with
  | nil => simp [groupsStartingInModel, sumOver]
  | cons g gs ih =>
    simp only [groupsStartingInModel, sumOver, ← ih]
    by_cases hg : first ≤ g ∧ g ≤ last <;> simp [hg]

/-- The atoms `start ..= last` of `n` atoms, `last < n`, are `last + 1 − start`. -/
theorem sumBelow_interval (n start last : Nat) (h : last < n) :
    sumBelow n (fun a => if start ≤ a ∧ a ≤ last then 1 else 0)
      = ((last + 1 - start : Nat) : Int) := by
  suffices ∀ m, sumBelow m (fun a => if start ≤ a ∧ a ≤ last then (1 : Int) else 0)
      = ((min m (last + 1) - start : Nat) : Int) by
    rw [this n, show min n (last + 1) = last + 1 by omega]
  intro m
  induction m with
  | zero => simp [sumBelow]
  | succ m ih =>
    simp only [sumBelow, ih]
    by_cases hm : start ≤ m ∧ m ≤ last
    · simp only [hm, and_self, ite_true]; omega
    · simp only [hm, ite_false]; omega

/-- The COM part of the range loop, from `start` on, is the theory's sum over group firsts. -/
theorem rangesHalf_com (bath : Nat) (gs : List Nat) {start : Nat} {rs : List Range} {n : Nat}
    (h : Covers start rs n) :
    (rangesHalfModel bath gs start rs).1
      = sumOver gs (fun g => if start ≤ g ∧ g < n ∧ comOf rs g = bath then 6 else 0) := by
  induction rs generalizing start with
  | nil =>
    simp only [Covers] at h
    simp only [rangesHalfModel]
    exact (sumOver_zero _ _ (fun g _ => by simp; omega)).symm
  | cons r rs ih =>
    simp only [Covers] at h
    have hle := covers_le h.2
    simp only [rangesHalfModel, ih h.2]
    rw [groupsStartingIn_eq]
    have e : (if r.comBath = bath then
        6 * sumOver gs (fun g => if start ≤ g ∧ g ≤ r.lastAtom then (1 : Int) else 0) else 0)
        = sumOver gs (fun g =>
            if start ≤ g ∧ g ≤ r.lastAtom ∧ r.comBath = bath then (6 : Int) else 0) := by
      by_cases hb : r.comBath = bath
      · simp only [hb, ite_true, ← sumOver_mul]
        apply sumOver_congr; intro g _
        by_cases hg : start ≤ g ∧ g ≤ r.lastAtom <;> simp [hg] <;> omega
      · simp only [hb, ite_false, and_false]
        exact (sumOver_zero _ _ (fun g _ => by simp)).symm
    rw [e, ← sumOver_add]
    apply sumOver_congr
    intro g _
    simp only [comOf, bathsOfModel]
    by_cases hg : g ≤ r.lastAtom
    · simp only [hg, ite_true]
      by_cases hs : start ≤ g
      · by_cases hb : r.comBath = bath <;> simp [hs, hb] <;> omega
      · simp [hs]; omega
    · simp only [hg, ite_false]
      by_cases hs : r.lastAtom + 1 ≤ g
      · simp [hs, show start ≤ g by omega]
      · omega

/-- The IR part of the range loop, from `start` on: every atom's 6 halves less each group's 6. -/
theorem rangesHalf_ir (bath : Nat) (gs : List Nat) {start : Nat} {rs : List Range} {n : Nat}
    (h : Covers start rs n) :
    (rangesHalfModel bath gs start rs).2
      = sumBelow n (fun a => if start ≤ a ∧ irOf rs a = bath then 6 else 0)
        - sumOver gs (fun g => if start ≤ g ∧ g < n ∧ irOf rs g = bath then 6 else 0) := by
  induction rs generalizing start with
  | nil =>
    simp only [Covers] at h
    simp only [rangesHalfModel]
    rw [sumOver_zero _ _ (fun g _ => by simp; omega),
      sumBelow_congr n _ (fun _ => 0) (fun a ha => by simp; omega), sumBelow_const]
    simp
  | cons r rs ih =>
    simp only [Covers] at h
    have hle := covers_le h.2
    simp only [rangesHalfModel, ih h.2]
    rw [groupsStartingIn_eq, ← sumBelow_interval n start r.lastAtom (by omega)]
    have e : (if r.irBath = bath then
        6 * sumBelow n (fun a => if start ≤ a ∧ a ≤ r.lastAtom then (1 : Int) else 0)
          - 6 * sumOver gs (fun g => if start ≤ g ∧ g ≤ r.lastAtom then (1 : Int) else 0)
        else 0)
        = sumBelow n (fun a =>
            if start ≤ a ∧ a ≤ r.lastAtom ∧ r.irBath = bath then (6 : Int) else 0)
          - sumOver gs (fun g =>
            if start ≤ g ∧ g ≤ r.lastAtom ∧ r.irBath = bath then (6 : Int) else 0) := by
      by_cases hb : r.irBath = bath
      · simp only [hb, ite_true, ← sumOver_mul, ← sumBelow_mul]
        congr 1
        · apply sumBelow_congr; intro a _
          by_cases ha : start ≤ a ∧ a ≤ r.lastAtom <;> simp [ha] <;> omega
        · apply sumOver_congr; intro g _
          by_cases hg : start ≤ g ∧ g ≤ r.lastAtom <;> simp [hg] <;> omega
      · simp only [hb, ite_false, and_false]
        rw [sumOver_zero _ _ (fun g _ => by simp),
          sumBelow_congr n _ (fun _ => 0) (fun a _ => by simp), sumBelow_const]
        simp
    rw [e]
    have ea : sumBelow n (fun a => if start ≤ a ∧ irOf (r :: rs) a = bath then (6 : Int) else 0)
        = sumBelow n (fun a => if start ≤ a ∧ a ≤ r.lastAtom ∧ r.irBath = bath then 6 else 0)
          + sumBelow n (fun a => if r.lastAtom + 1 ≤ a ∧ irOf rs a = bath then 6 else 0) := by
      rw [← sumBelow_add]
      apply sumBelow_congr; intro a _
      simp only [irOf, bathsOfModel]
      by_cases ha : a ≤ r.lastAtom
      · simp only [ha, ite_true]
        by_cases hs : start ≤ a
        · by_cases hb : r.irBath = bath <;> simp [hs, hb] <;> omega
        · simp [hs]; omega
      · simp only [ha, ite_false]
        by_cases hs : r.lastAtom + 1 ≤ a
        · simp [hs, show start ≤ a by omega]
        · omega
    have eg : sumOver gs (fun g => if start ≤ g ∧ g < n ∧ irOf (r :: rs) g = bath then (6 : Int)
        else 0)
        = sumOver gs (fun g => if start ≤ g ∧ g ≤ r.lastAtom ∧ r.irBath = bath then 6 else 0)
          + sumOver gs (fun g => if r.lastAtom + 1 ≤ g ∧ g < n ∧ irOf rs g = bath then 6 else 0) := by
      rw [← sumOver_add]
      apply sumOver_congr; intro g _
      simp only [irOf, bathsOfModel]
      by_cases hg : g ≤ r.lastAtom
      · simp only [hg, ite_true]
        by_cases hs : start ≤ g
        · by_cases hb : r.irBath = bath <;> simp [hs, hb] <;> omega
        · simp [hs]; omega
      · simp only [hg, ite_false]
        by_cases hs : r.lastAtom + 1 ≤ g
        · simp [hs, show start ≤ g by omega]
        · omega
    rw [ea, eg]
    omega

/-- A constraint end is counted in bath `bath` exactly when its atom is one of the `n` atoms and
its IR bath is `bath`. -/
theorem endIn_eq {rs : List Range} {n : Nat} (h : Covers 0 rs n) (bath a : Nat) :
    endInModel bath rs a = if a < n ∧ irOf rs a = bath then 1 else 0 := by
  have hs := bathsOf_some_iff h a (Nat.zero_le a)
  by_cases ha : a < n
  · have e := bathsOf_eq a (hs.mpr ha)
    simp [endInModel, e, ha]
  · have hn : bathsOfModel rs a = none := by
      cases e : bathsOfModel rs a with
      | none => rfl
      | some p => exact absurd (hs.mp (by simp [e])) ha
    simp [endInModel, hn, ha]

theorem ends_eq {rs : List Range} {n : Nat} (h : Covers 0 rs n) (bath : Nat)
    (cs : List (Nat × Nat)) :
    endsModel bath rs cs = sumOver cs (fun c =>
      (if c.1 < n ∧ irOf rs c.1 = bath then 1 else 0)
        + (if c.2 < n ∧ irOf rs c.2 = bath then 1 else 0)) := by
  induction cs with
  | nil => rfl
  | cons c cs ih => simp only [endsModel, sumOver, ih, endIn_eq h, Int.add_assoc]

/-- **The range-by-range count is the atom-by-atom count.** When the ranges cover the `n` atoms
in order, `bath_half_dof`'s model is the theory's `(comHalf, irHalf)` of the system whose atoms
have the baths the ranges give. -/
theorem bathHalfDofModel_eq {ranges : List Range} {n : Nat} (h : Covers 0 ranges n)
    (bath : Nat) (groupFirsts : List Nat) (constraints : List (Nat × Nat)) (rottrans : Bool) :
    bathHalfDofModel bath ranges groupFirsts constraints rottrans
      = (comHalf (systemOf n ranges groupFirsts constraints rottrans) bath,
         irHalf (systemOf n ranges groupFirsts constraints rottrans) bath) := by
  have hcom := rangesHalf_com bath groupFirsts h
  have hir := rangesHalf_ir bath groupFirsts h
  have hends := ends_eq h bath constraints
  simp only [Nat.zero_le, true_and] at hcom hir
  -- atom 0's baths: `some` exactly when there are atoms
  have h0 := bathsOf_some_iff h 0 (Nat.le_refl 0)
  simp only [bathHalfDofModel, comHalf, irHalf, systemOf, BathSystem.rottransOn, hcom, hir,
    hends]
  by_cases hn : 0 < n
  · rw [bathsOf_eq 0 (h0.mpr hn)]
    cases rottrans
    · simp only [Bool.false_and, Bool.false_eq_true, false_and, ite_false, Int.sub_zero]
      first | rfl | exact ⟨rfl, rfl⟩
    · by_cases hc : comOf ranges 0 = bath <;> by_cases hi : irOf ranges 0 = bath <;>
        simp [hc, hi, hn] <;> first | rfl | exact ⟨rfl, rfl⟩
  · have hnone : bathsOfModel ranges 0 = none := by
      cases e : bathsOfModel ranges 0 with
      | none => rfl
      | some p => exact absurd (h0.mp (by simp [e])) hn
    rw [hnone]
    cases rottrans <;> simp [hn] <;> first | rfl | exact ⟨rfl, rfl⟩

/-- The same under the plan's validation. -/
theorem bathHalfDofModel_eq_of_wellFormed {nBaths n : Nat} {ranges : List Range}
    {groupFirsts : List Nat} (h : RangesWellFormed nBaths n ranges groupFirsts) (bath : Nat)
    (constraints : List (Nat × Nat)) (rottrans : Bool) :
    bathHalfDofModel bath ranges groupFirsts constraints rottrans
      = (comHalf (systemOf n ranges groupFirsts constraints rottrans) bath,
         irHalf (systemOf n ranges groupFirsts constraints rottrans) bath) :=
  bathHalfDofModel_eq h.1 bath groupFirsts constraints rottrans

/-- Every bath of `bath_half_dofs`' model is the theory's, in order. -/
theorem bathHalfDofsModel_eq {nBaths n : Nat} {ranges : List Range} {groupFirsts : List Nat}
    (h : RangesWellFormed nBaths n ranges groupFirsts) (constraints : List (Nat × Nat))
    (rottrans : Bool) :
    bathHalfDofsModel nBaths ranges groupFirsts constraints rottrans
      = (List.range nBaths).map (fun b =>
          (comHalf (systemOf n ranges groupFirsts constraints rottrans) b,
           irHalf (systemOf n ranges groupFirsts constraints rottrans) b)) := by
  simp only [bathHalfDofsModel, bathHalfDofModel_eq_of_wellFormed h]

/-- **The spread is the theory's spread.** -/
theorem ndfminSpreadModel_eq (halfDofs : List Int) (ndfmin : Int) :
    ndfminSpreadModel halfDofs ndfmin
      = (spreadNumerators halfDofs ndfmin, spreadDenominator halfDofs) := by
  simp only [ndfminSpreadModel, spreadNumerators, spreadFactor, spreadDenominator]
  by_cases h : 0 < halfDofs.sum <;> simp [h]


/-- **`empty_baths` on `bath_half_dofs`' totals flags the theory's empty baths**, in order:
called as `newen-run` calls it, on `com + ir` of every bath before the spread, it says which
baths gromosXX decouples. -/
theorem emptyBathsModel_eq {nBaths n : Nat} {ranges : List Range} {groupFirsts : List Nat}
    (h : RangesWellFormed nBaths n ranges groupFirsts) (constraints : List (Nat × Nat))
    (rottrans : Bool) :
    emptyBathsModel ((bathHalfDofsModel nBaths ranges groupFirsts constraints rottrans).map
        fun p => p.1 + p.2)
      = (List.range nBaths).map (fun b =>
          decide (EmptyBath (systemOf n ranges groupFirsts constraints rottrans) b)) := by
  rw [bathHalfDofsModel_eq h]
  simp only [emptyBathsModel, List.map_map]
  rfl

end NewenSpec.Thermo.Bath
