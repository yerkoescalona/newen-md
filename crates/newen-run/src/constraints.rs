//! Which constraint algorithm(s) a GROMOS input selects.
//!
//! gromosXX dispatches the *solute* algorithm on `NTCP` (only relevant when `NTC > 1`) and the
//! *solvent* algorithm on `NTCS` (only relevant when solvent molecules exist), independently.
//! Moved here from `pyo3-newen` (where it was one of three copies of this logic).

use newen_integrators::constraints::NtcMode;
use newen_io::imd::ImdParameters;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConstraintSelection {
    pub solute_shake: bool,
    pub solvent_shake: bool,
    pub settle_enabled: bool,
    pub lincs_enabled: bool,
    pub solute_lincs: bool,
    pub solvent_lincs: bool,
}

impl ConstraintSelection {
    /// `has_solvent` must reflect the topology's *actual* solvent atoms
    /// (`num_atoms() > num_solute_atoms()`), not `imd.nsm`: the Python object path can solvate
    /// a topology directly with factory-built parameters that never set `nsm`, and relying on
    /// `nsm` alone silently disabled solvent constraints there (rigid water flying apart).
    pub fn from_imd(imd: &ImdParameters, has_solvent: bool) -> Self {
        Self::from_codes(imd.ntc, imd.ntcp, imd.ntcs, has_solvent)
    }

    /// The same dispatch from the raw GROMOS codes (NTC, NTCP, NTCS) — used by the recipe plan.
    pub fn from_codes(ntc: i32, ntcp: i32, ntcs: i32, has_solvent: bool) -> Self {
        let solute_constrained = ntc > 1;
        let solute_lincs = solute_constrained && ntcp == 2;
        let solute_shake = solute_constrained && !solute_lincs;

        let solvent_constrained = ntcs > 0 && has_solvent;
        let solvent_settle = solvent_constrained && ntcs == 3;
        let solvent_lincs = solvent_constrained && ntcs == 2;
        let solvent_shake = solvent_constrained && !solvent_settle && !solvent_lincs;

        ConstraintSelection {
            solute_shake,
            solvent_shake,
            settle_enabled: solvent_settle,
            lincs_enabled: solute_lincs || solvent_lincs,
            solute_lincs,
            solvent_lincs,
        }
    }

    pub fn shake_enabled(&self) -> bool {
        self.solute_shake || self.solvent_shake
    }

    /// Any constraint algorithm at all (drives constraint-force output and the DOF count).
    pub fn any(&self) -> bool {
        self.shake_enabled() || self.settle_enabled || self.lincs_enabled
    }

    pub fn solute_constrained(&self) -> bool {
        self.solute_shake || self.solute_lincs
    }

    pub fn solvent_constrained(&self) -> bool {
        self.solvent_shake || self.settle_enabled || self.solvent_lincs
    }

    /// The solute constraint mode implied by `NTC`.
    pub fn ntc_mode(imd: &ImdParameters) -> NtcMode {
        Self::ntc_mode_of(imd.ntc)
    }

    pub fn ntc_mode_of(ntc: i32) -> NtcMode {
        match ntc {
            3 => NtcMode::AllBonds,
            2 => NtcMode::HydrogenBonds,
            _ => NtcMode::SolventOnly,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ntc_1_with_settle_water() {
        let s = ConstraintSelection::from_codes(1, 1, 3, true);
        assert!(!s.solute_shake && !s.solute_lincs);
        assert!(s.settle_enabled && !s.solvent_shake && !s.lincs_enabled);
        assert!(s.any() && s.solvent_constrained() && !s.solute_constrained());
    }

    #[test]
    fn ntc_3_shake_everywhere() {
        let s = ConstraintSelection::from_codes(3, 1, 1, true);
        assert!(s.solute_shake && s.solvent_shake && s.shake_enabled());
        assert!(!s.settle_enabled && !s.lincs_enabled);
    }

    #[test]
    fn ntcp_2_selects_lincs_for_the_solute() {
        let s = ConstraintSelection::from_codes(2, 2, 2, true);
        assert!(s.solute_lincs && s.solvent_lincs && s.lincs_enabled);
        assert!(!s.shake_enabled());
    }

    #[test]
    fn no_solvent_means_no_solvent_constraints() {
        let s = ConstraintSelection::from_codes(1, 1, 3, false);
        assert!(!s.any());
    }
}
