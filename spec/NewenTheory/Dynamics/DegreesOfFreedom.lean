/-!
# Kinetic degrees of freedom

A system of `N` point particles in three dimensions has `3 N` coordinates. Every independent
holonomic constraint (a fixed bond length, a rigid water geometry) removes one of them. Overall
motion that the simulation does not let the particles exchange with a heat bath removes more: a
fixed or removed centre-of-mass translation removes 3, and GROMOS takes the number to remove as
`NDFMIN` (BOUNDCOND block). What is left is the number of kinetic degrees of freedom `N_df`.

This number matters because of equipartition: in equilibrium every quadratic kinetic degree of
freedom carries `½ k_B T` on average, so `⟨K⟩ = ½ N_df k_B T`. The instantaneous temperature
`T = 2 K / (N_df k_B)` and every thermostat that drives `K` towards a target therefore depend on
counting `N_df` correctly. Equipartition itself is not proved here; only the counting is.

A system is described molecule by molecule: each `Molecule` has its atoms and its constraints.
The count is the sum over molecules of `3 · atoms − constraints`, minus `NDFMIN` for the whole
system. It is an integer (`Int`) on purpose: more constraints than coordinates is not a physical
system, and the count then becomes negative rather than meaningless.

Only core Lean is used (`omega` and `Int` lemmas), no Mathlib.
-/

namespace NewenTheory.Dynamics.DegreesOfFreedom

/-- The number of coordinates of `n` point particles in three dimensions. -/
def coordinates (n : Nat) : Nat := 3 * n

/-- A molecule as far as counting goes: its atoms and the independent holonomic constraints
imposed on them. -/
structure Molecule where
  atoms : Nat
  constraints : Nat

/-- The degrees of freedom one molecule keeps: its coordinates less its constraints. -/
def Molecule.dof (m : Molecule) : Int := (coordinates m.atoms : Int) - m.constraints

/-- All atoms of a collection of molecules. -/
def atoms : List Molecule → Nat
  | [] => 0
  | m :: ms => m.atoms + atoms ms

/-- All constraints of a collection of molecules. -/
def constraints : List Molecule → Nat
  | [] => 0
  | m :: ms => m.constraints + constraints ms

/-- The sum over molecules of the degrees of freedom each keeps. -/
def internalDof : List Molecule → Int
  | [] => 0
  | m :: ms => m.dof + internalDof ms

/-- The kinetic degrees of freedom of a system of molecules from which `ndfmin` degrees of
freedom of overall motion are removed. Physically `ndfmin ≥ 0`; it is an `Int` so that the count
is defined for every input. -/
def dof (ms : List Molecule) (ndfmin : Int) : Int := internalDof ms - ndfmin

/-- The molecule-by-molecule sum is the global count: `3 N` minus all constraints. -/
theorem internalDof_eq (ms : List Molecule) :
    internalDof ms = (coordinates (atoms ms) : Int) - constraints ms := by
  induction ms with
  | nil => simp [internalDof, atoms, constraints, coordinates]
  | cons m ms ih =>
    simp only [internalDof, atoms, constraints, Molecule.dof, coordinates, ih]
    omega

/-- **Additivity**: the degrees of freedom of two disjoint subsystems add up. -/
theorem internalDof_append (ms₁ ms₂ : List Molecule) :
    internalDof (ms₁ ++ ms₂) = internalDof ms₁ + internalDof ms₂ := by
  induction ms₁ with
  | nil => simp [internalDof]
  | cons m ms ih =>
    simp only [List.cons_append, internalDof, ih]
    omega

/-- Additivity with the overall motion removed once, from the whole system. -/
theorem dof_append (ms₁ ms₂ : List Molecule) (ndfmin : Int) :
    dof (ms₁ ++ ms₂) ndfmin = internalDof ms₁ + internalDof ms₂ - ndfmin := by
  rw [dof, internalDof_append]

/-- **Identical molecules**: `k` copies of one molecule keep `k` times its degrees of freedom.
This is how a solvent stored as one template molecule repeated `k` times is counted. -/
theorem internalDof_replicate (k : Nat) (m : Molecule) :
    internalDof (List.replicate k m) = (k : Int) * m.dof := by
  induction k with
  | zero => simp [internalDof]
  | succ k ih =>
    rw [List.replicate_succ, internalDof, ih, Int.natCast_succ, Int.add_mul, Int.one_mul,
      Int.add_comm]

/-- The atoms of `k` identical molecules. -/
theorem atoms_replicate (k : Nat) (m : Molecule) :
    atoms (List.replicate k m) = k * m.atoms := by
  induction k with
  | zero => simp [atoms]
  | succ k ih => rw [List.replicate_succ, atoms, ih, Nat.succ_mul, Nat.add_comm]

/-- The constraints of `k` identical molecules. -/
theorem constraints_replicate (k : Nat) (m : Molecule) :
    constraints (List.replicate k m) = k * m.constraints := by
  induction k with
  | zero => simp [constraints]
  | succ k ih => rw [List.replicate_succ, constraints, ih, Nat.succ_mul, Nat.add_comm]

/-- The atoms of two subsystems add up. -/
theorem atoms_append (ms₁ ms₂ : List Molecule) : atoms (ms₁ ++ ms₂) = atoms ms₁ + atoms ms₂ := by
  induction ms₁ with
  | nil => simp [atoms]
  | cons m ms ih => simp only [List.cons_append, atoms, ih]; omega

/-- The constraints of two subsystems add up. -/
theorem constraints_append (ms₁ ms₂ : List Molecule) :
    constraints (ms₁ ++ ms₂) = constraints ms₁ + constraints ms₂ := by
  induction ms₁ with
  | nil => simp [constraints]
  | cons m ms ih => simp only [List.cons_append, constraints, ih]; omega

/-- **Physicality**: the count is non-negative exactly when the constraints and the removed
overall motion together do not exceed the `3 N` coordinates. -/
theorem dof_nonneg_iff (ms : List Molecule) (ndfmin : Int) :
    0 ≤ dof ms ndfmin ↔ (constraints ms : Int) + ndfmin ≤ coordinates (atoms ms) := by
  rw [dof, internalDof_eq]
  omega

/-! ## Temperature baths

A thermostat with several baths (GROMOS `MULTIBATH`) couples every atom to two of them: the
centre-of-mass translation of the atom's temperature group to a *COM bath*, everything else
(internal motion and rotation) to an *IR bath*. Equipartition then holds per bath, so each bath
needs its own `N_df`. Here the bath of every atom is given directly, as two functions; how a
`MULTIBATH` block describes it (ranges of atoms) is the code's business and is bridged to this in
`NewenSpec/Thermo/Bath.lean`.

The counts are in *half* degrees of freedom: a distance constraint between atoms in two baths
removes ½ from each, and in halves that stays an integer. Per bath `b`:

* every atom gives its 3 degrees of freedom (6 halves) to its IR bath;
* every temperature group moves the 3 of its centre-of-mass translation (6 halves) from the IR
  bath of its first atom to the COM bath of that atom; for groups given by distinct first atoms
  this is "a group's first atom gives 6 to its COM bath, every other atom 6 to its IR bath";
* every constraint end takes 1 half from the IR bath of its atom;
* roto-translational constraints take 3 (6 halves) from each bath of atom 0.

`NDFMIN` is then spread over the baths in proportion to their degrees of freedom
(`spreadNumerators`, `spreadDenominator`).
-/

/-- `Σ_{a < n} f a`. -/
def sumBelow : Nat → (Nat → Int) → Int
  | 0, _ => 0
  | n + 1, f => sumBelow n f + f n

/-- `Σ_{x ∈ xs} f x`. -/
def sumOver {α : Type} : List α → (α → Int) → Int
  | [], _ => 0
  | x :: xs, f => f x + sumOver xs f

/-- A system as the thermostat sees it: `atoms` atoms, the COM and IR bath of each, the first
atom of every temperature group, the distance constraints (pairs of atoms), and whether
roto-translational constraints are on. -/
structure BathSystem where
  atoms : Nat
  comBath : Nat → Nat
  irBath : Nat → Nat
  groupFirsts : List Nat
  constraints : List (Nat × Nat)
  rottrans : Bool

/-- Roto-translational constraints act only on a non-empty system. -/
def BathSystem.rottransOn (s : BathSystem) : Bool := s.rottrans && 0 < s.atoms

/-- Half degrees of freedom of the centre-of-mass translation in bath `b`. -/
def comHalf (s : BathSystem) (b : Nat) : Int :=
  sumOver s.groupFirsts (fun g => if g < s.atoms ∧ s.comBath g = b then 6 else 0)
    - (if s.rottransOn ∧ s.comBath 0 = b then 6 else 0)

/-- Half degrees of freedom of the internal and rotational motion in bath `b`. -/
def irHalf (s : BathSystem) (b : Nat) : Int :=
  sumBelow s.atoms (fun a => if s.irBath a = b then 6 else 0)
    - sumOver s.groupFirsts (fun g => if g < s.atoms ∧ s.irBath g = b then 6 else 0)
    - sumOver s.constraints (fun c =>
        (if c.1 < s.atoms ∧ s.irBath c.1 = b then 1 else 0)
          + (if c.2 < s.atoms ∧ s.irBath c.2 = b then 1 else 0))
    - (if s.rottransOn ∧ s.irBath 0 = b then 6 else 0)

theorem sumBelow_add (n : Nat) (f g : Nat → Int) :
    sumBelow n (fun a => f a + g a) = sumBelow n f + sumBelow n g := by
  induction n with
  | zero => simp [sumBelow]
  | succ n ih => simp only [sumBelow, ih]; omega

theorem sumBelow_sub (n : Nat) (f g : Nat → Int) :
    sumBelow n (fun a => f a - g a) = sumBelow n f - sumBelow n g := by
  induction n with
  | zero => simp [sumBelow]
  | succ n ih => simp only [sumBelow, ih]; omega

theorem sumBelow_congr (n : Nat) (f g : Nat → Int) (h : ∀ a, a < n → f a = g a) :
    sumBelow n f = sumBelow n g := by
  induction n with
  | zero => rfl
  | succ n ih =>
    simp only [sumBelow]
    rw [ih (fun a ha => h a (by omega)), h n (by omega)]

theorem sumBelow_const (n : Nat) (c : Int) : sumBelow n (fun _ => c) = n * c := by
  induction n with
  | zero => simp [sumBelow]
  | succ n ih => simp only [sumBelow, ih]; rw [Int.natCast_succ, Int.add_mul, Int.one_mul]

theorem sumOver_add {α : Type} (xs : List α) (f g : α → Int) :
    sumOver xs (fun x => f x + g x) = sumOver xs f + sumOver xs g := by
  induction xs with
  | nil => simp [sumOver]
  | cons x xs ih => simp only [sumOver, ih]; omega

theorem sumOver_congr {α : Type} (xs : List α) (f g : α → Int) (h : ∀ x, x ∈ xs → f x = g x) :
    sumOver xs f = sumOver xs g := by
  induction xs with
  | nil => rfl
  | cons x xs ih =>
    simp only [sumOver]
    rw [h x (List.mem_cons_self ..), ih (fun y hy => h y (List.mem_cons_of_mem _ hy))]

theorem sumOver_const {α : Type} (xs : List α) (c : Int) :
    sumOver xs (fun _ => c) = xs.length * c := by
  induction xs with
  | nil => simp [sumOver]
  | cons x xs ih =>
    simp only [sumOver, ih, List.length_cons]
    rw [Int.natCast_succ, Int.add_mul, Int.one_mul, Int.add_comm]

/-- Sums over baths and sums over atoms commute. -/
theorem sumBelow_sumBelow (n m : Nat) (f : Nat → Nat → Int) :
    sumBelow n (fun b => sumBelow m (fun a => f a b))
      = sumBelow m (fun a => sumBelow n (fun b => f a b)) := by
  induction m with
  | zero => simp only [sumBelow]; exact (sumBelow_const n 0).trans (by simp)
  | succ m ih => simp only [sumBelow]; rw [sumBelow_add, ih]

/-- Sums over baths and sums over a list commute. -/
theorem sumBelow_sumOver {α : Type} (n : Nat) (xs : List α) (f : α → Nat → Int) :
    sumBelow n (fun b => sumOver xs (fun x => f x b))
      = sumOver xs (fun x => sumBelow n (fun b => f x b)) := by
  induction xs with
  | nil => simp only [sumOver]; exact (sumBelow_const n 0).trans (by simp)
  | cons x xs ih => simp only [sumOver]; rw [sumBelow_add, ih]

/-- Summed over all baths `b < n`, "`c` if `x = b`" gives `c` once, if `x` is a bath. -/
theorem sumBelow_ite_eq (n x : Nat) (c : Int) :
    sumBelow n (fun b => if x = b then c else 0) = if x < n then c else 0 := by
  induction n with
  | zero => simp [sumBelow]
  | succ n ih =>
    simp only [sumBelow, ih]
    by_cases h1 : x < n
    · simp [h1, show x ≠ n by omega, show x < n + 1 by omega]
    · by_cases h2 : x = n
      · simp [h2]
      · simp [h1, h2, show ¬ x < n + 1 by omega]

/-- The same, with a side condition `p` under which `x` is a bath. -/
theorem sumBelow_ite_and_eq (n x : Nat) (p : Prop) [Decidable p] (c : Int) (hx : p → x < n) :
    sumBelow n (fun b => if p ∧ x = b then c else 0) = if p then c else 0 := by
  by_cases hp : p
  · simp only [hp, true_and, ite_true]
    rw [sumBelow_ite_eq]
    simp [hx hp]
  · simp only [hp, false_and, ite_false]
    exact (sumBelow_const n 0).trans (by simp)

/-- **The baths share out the whole system.** When every atom's two baths are among `0 .. n−1`
and every constraint joins two atoms of the system, the half degrees of freedom of all baths add
up to twice `3 N − C` (atoms `N`, constraints `C`), less 12 halves for the roto-translational
constraints. -/
theorem bathHalf_partition (s : BathSystem) (n : Nat)
    (hcom : ∀ a, a < s.atoms → s.comBath a < n) (hir : ∀ a, a < s.atoms → s.irBath a < n)
    (hcon : ∀ c, c ∈ s.constraints → c.1 < s.atoms ∧ c.2 < s.atoms) :
    sumBelow n (fun b => comHalf s b + irHalf s b)
      = 2 * (3 * (s.atoms : Int) - s.constraints.length)
          - (if s.rottransOn then 12 else 0) := by
  have hrot : s.rottransOn → 0 < s.atoms := by
    simp [BathSystem.rottransOn]
  simp only [comHalf, irHalf, sumBelow_add, sumBelow_sub, sumBelow_sumOver]
  rw [sumBelow_sumBelow n s.atoms, sumBelow_ite_and_eq n _ _ _ (fun h => hcom 0 (hrot h)),
    sumBelow_ite_and_eq n _ _ _ (fun h => hir 0 (hrot h))]
  have e1 : sumOver s.groupFirsts (fun g => sumBelow n (fun b =>
      if g < s.atoms ∧ s.comBath g = b then (6 : Int) else 0))
      = sumOver s.groupFirsts (fun g => if g < s.atoms then 6 else 0) :=
    sumOver_congr _ _ _ (fun g _ => sumBelow_ite_and_eq n _ _ _ (hcom g))
  have e2 : sumOver s.groupFirsts (fun g => sumBelow n (fun b =>
      if g < s.atoms ∧ s.irBath g = b then (6 : Int) else 0))
      = sumOver s.groupFirsts (fun g => if g < s.atoms then 6 else 0) :=
    sumOver_congr _ _ _ (fun g _ => sumBelow_ite_and_eq n _ _ _ (hir g))
  have e3 : sumBelow s.atoms (fun a => sumBelow n (fun b =>
      if s.irBath a = b then (6 : Int) else 0)) = sumBelow s.atoms (fun _ => 6) :=
    sumBelow_congr _ _ _ (fun a ha => by rw [sumBelow_ite_eq]; simp [hir a ha])
  have e4 : sumOver s.constraints (fun c =>
      sumBelow n (fun b => if c.1 < s.atoms ∧ s.irBath c.1 = b then (1 : Int) else 0)
        + sumBelow n (fun b => if c.2 < s.atoms ∧ s.irBath c.2 = b then 1 else 0))
      = sumOver s.constraints (fun _ => 2) := by
    apply sumOver_congr
    intro c hc
    obtain ⟨h1, h2⟩ := hcon c hc
    rw [sumBelow_ite_and_eq n _ _ _ (fun _ => hir _ h1),
      sumBelow_ite_and_eq n _ _ _ (fun _ => hir _ h2)]
    simp [h1, h2]
  rw [e1, e2, e3, e4, sumBelow_const, sumOver_const]
  by_cases h : s.rottransOn = true
  · simp only [h, ite_true]; omega
  · simp only [h, Bool.false_eq_true, ite_false]; omega

/-- **One bath is the whole system.** With every atom in bath 0 for both parts, the bath holds
twice the system's degrees of freedom of step 1 (`dof ms 0`, no `NDFMIN`) less the 12 halves of
the roto-translational constraints. -/
theorem bathHalf_single (s : BathSystem) (ms : List Molecule)
    (hcom : ∀ a, s.comBath a = 0) (hir : ∀ a, s.irBath a = 0)
    (hcon : ∀ c, c ∈ s.constraints → c.1 < s.atoms ∧ c.2 < s.atoms)
    (hatoms : atoms ms = s.atoms) (hconstraints : constraints ms = s.constraints.length) :
    comHalf s 0 + irHalf s 0 = 2 * dof ms 0 - (if s.rottransOn then 12 else 0) := by
  have h := bathHalf_partition s 1 (fun a _ => by rw [hcom]; omega)
    (fun a _ => by rw [hir]; omega) hcon
  simp only [sumBelow] at h
  rw [dof, internalDof_eq, hatoms, hconstraints, coordinates]
  push_cast
  omega

/-! ### Spreading NDFMIN over the baths

`NDFMIN` degrees of freedom of overall motion are removed from the system as a whole, and gromosXX
shares the removal out over the baths in proportion to each bath's degrees of freedom. With half
counts `h_b` and their total `H = Σ_b h_b > 0`, bath `b` keeps `(h_b / 2) · (1 − 2 · NDFMIN / H)`,
the exact rational `h_b · (H − 2 · NDFMIN) / (2 H)`: numerators over one common denominator. With
`H ≤ 0` there is nothing to spread over, and every bath keeps `h_b / 2`. -/

/-- The factor every half count is multiplied by. -/
def spreadFactor (h : List Int) (ndfmin : Int) : Int :=
  if 0 < h.sum then h.sum - 2 * ndfmin else 1

/-- The numerator of the degrees of freedom bath `b` keeps after the spread. -/
def spreadNumerators (h : List Int) (ndfmin : Int) : List Int :=
  h.map (fun x => x * spreadFactor h ndfmin)

/-- The common denominator. -/
def spreadDenominator (h : List Int) : Int := if 0 < h.sum then 2 * h.sum else 2

theorem sum_map_mul (h : List Int) (k : Int) : (h.map (fun x => x * k)).sum = h.sum * k := by
  induction h with
  | nil => simp
  | cons x xs ih => simp only [List.map_cons, List.sum_cons, ih, Int.add_mul]

/-- **The spread removes exactly `NDFMIN`.** With `H > 0` the numerators add up to
`(H − 2·NDFMIN)·H` over the denominator `2H`, i.e. `Σ_b num_b / den = H/2 − NDFMIN`: the baths
together keep the system's degrees of freedom less `NDFMIN`. -/
theorem spread_sum (h : List Int) (ndfmin : Int) (hpos : 0 < h.sum) :
    (spreadNumerators h ndfmin).sum = (h.sum - 2 * ndfmin) * h.sum
      ∧ spreadDenominator h = 2 * h.sum := by
  refine ⟨?_, ?_⟩
  · simp only [spreadNumerators, sum_map_mul, spreadFactor, hpos, ite_true]
    exact Int.mul_comm _ _
  · simp [spreadDenominator, hpos]

/-- **A bath without degrees of freedom keeps none.** -/
theorem spread_zero (h : List Int) (ndfmin : Int) (i : Nat) (hi : h[i]? = some 0) :
    (spreadNumerators h ndfmin)[i]? = some 0 := by
  simp [spreadNumerators, List.getElem?_map, hi]

/-- **Nothing to spread over.** With `H ≤ 0` every bath keeps its own `h_b / 2`. -/
theorem spread_total_nonpos (h : List Int) (ndfmin : Int) (hnp : h.sum ≤ 0) :
    spreadNumerators h ndfmin = h ∧ spreadDenominator h = 2 := by
  have hn : ¬ 0 < h.sum := by omega
  simp [spreadNumerators, spreadDenominator, spreadFactor, hn]


/-! ### Empty baths

A bath with no degrees of freedom cannot be coupled: its temperature is undefined. gromosXX
decouples such a bath (`tau = -1`, "removing coupling of bath with 0 dof",
`topology.cc:1404-1411`). It decides on the count *before* NDFMIN is spread over the baths
(`multibath.cc:139-153` runs after it); `emptyBath_iff` shows that this is the same as deciding on
the spread numerators, except when NDFMIN takes all the system's degrees of freedom
(`0 < H ≤ 2·NDFMIN`), which newen refuses. -/

/-- Bath `b` has no degrees of freedom before NDFMIN: `comHalf + irHalf ≤ 0` halves. -/
def EmptyBath (s : BathSystem) (b : Nat) : Prop := comHalf s b + irHalf s b ≤ 0

instance (s : BathSystem) (b : Nat) : Decidable (EmptyBath s b) :=
  inferInstanceAs (Decidable (_ ≤ _))

/-- **A bath no atom belongs to is empty**: it has exactly no degrees of freedom. -/
theorem emptyBath_of_unused (s : BathSystem) (b : Nat)
    (h : ∀ a, a < s.atoms → s.comBath a ≠ b ∧ s.irBath a ≠ b) :
    comHalf s b + irHalf s b = 0 ∧ EmptyBath s b := by
  have hz : ∀ (f : Nat → Nat), (∀ a, a < s.atoms → f a ≠ b) → ∀ (c : Int) (xs : List Nat),
      sumOver xs (fun g => if g < s.atoms ∧ f g = b then c else 0) = 0 := by
    intro f hf c xs
    induction xs with
    | nil => rfl
    | cons x xs ih =>
      simp only [sumOver, ih]
      by_cases hx : x < s.atoms
      · simp [hf x hx]
      · simp [hx]
  have hcom := hz s.comBath (fun a ha => (h a ha).1) 6 s.groupFirsts
  have hir := hz s.irBath (fun a ha => (h a ha).2) 6 s.groupFirsts
  have hatoms : sumBelow s.atoms (fun a => if s.irBath a = b then 6 else 0) = 0 := by
    rw [sumBelow_congr s.atoms _ (fun _ => 0) (fun a ha => by simp [(h a ha).2])]
    rw [sumBelow_const]; simp
  have hcs : sumOver s.constraints (fun c =>
      (if c.1 < s.atoms ∧ s.irBath c.1 = b then 1 else 0)
        + (if c.2 < s.atoms ∧ s.irBath c.2 = b then 1 else 0)) = 0 := by
    rw [sumOver_congr s.constraints _ (fun _ => 0)]
    · rw [sumOver_const]; simp
    · intro c _
      have e1 : ¬ (c.1 < s.atoms ∧ s.irBath c.1 = b) := fun ⟨ha, he⟩ => (h c.1 ha).2 he
      have e2 : ¬ (c.2 < s.atoms ∧ s.irBath c.2 = b) := fun ⟨ha, he⟩ => (h c.2 ha).2 he
      simp only [e1, e2, ite_false]; rfl
  have hrt : ∀ f : Nat → Nat, (∀ a, a < s.atoms → f a ≠ b) →
      ¬ (s.rottransOn = true ∧ f 0 = b) := by
    intro f hf ⟨hon, he⟩
    simp only [BathSystem.rottransOn, Bool.and_eq_true, decide_eq_true_eq] at hon
    exact hf 0 hon.2 he
  have e := hrt s.comBath (fun a ha => (h a ha).1)
  have e' := hrt s.irBath (fun a ha => (h a ha).2)
  have hsum : comHalf s b + irHalf s b = 0 := by
    simp only [comHalf, irHalf, hcom, hir, hatoms, hcs, e, e', ite_false]; rfl
  exact ⟨hsum, by unfold EmptyBath; omega⟩

theorem nonpos_iff_mul_nonpos (x : Int) {f : Int} (hf : 0 < f) : x ≤ 0 ↔ x * f ≤ 0 := by
  constructor
  · intro hx
    have := Int.mul_le_mul_of_nonneg_right hx (Int.le_of_lt hf)
    simpa using this
  · intro hx
    apply Classical.byContradiction
    intro hn
    have := Int.mul_pos (Int.lt_of_not_ge hn) hf
    omega

/-- **Deciding before the spread is deciding after it**, whenever NDFMIN leaves some degrees of
freedom (`2·NDFMIN < H`) or there are none to spread over (`H ≤ 0`): bath `i` has no degrees of
freedom before the spread exactly when its numerator after the spread is `≤ 0`. -/
theorem emptyBath_iff (h : List Int) (ndfmin : Int)
    (hv : 2 * ndfmin < h.sum ∨ h.sum ≤ 0) (i : Nat) (hi : i < h.length) :
    h[i] ≤ 0 ↔ (spreadNumerators h ndfmin)[i]'(by simpa [spreadNumerators] using hi) ≤ 0 := by
  simp only [spreadNumerators, List.getElem_map]
  have hf : 0 < spreadFactor h ndfmin := by
    unfold spreadFactor; split <;> omega
  exact nonpos_iff_mul_nonpos h[i] hf

end NewenTheory.Dynamics.DegreesOfFreedom
