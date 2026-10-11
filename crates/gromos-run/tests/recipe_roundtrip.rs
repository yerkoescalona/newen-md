//! `RunRecipe` round trips (3.9 G3/G7, step 2): `.imd` → recipe → `.imd` → recipe is
//! the identity on every reference input; the recipe survives TOML and JSON; the defaults are
//! gromosXX's absent-block behaviour; unmodelled blocks are refused unless allowed.

use std::path::{Path, PathBuf};

use gromos_io::imd::{parse_imd_str, read_imd_file, ImdParameters};
use gromos_run::recipe::{Bath, PassthroughPolicy, Thermostat, ThermostatAlgorithm};
use gromos_run::{RunError, RunRecipe, RECIPE_VERSION};

fn reference_inputs() -> Vec<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references");
    let mut files: Vec<PathBuf> = std::fs::read_dir(root)
        .expect("reference dir")
        .flatten()
        .filter(|e| e.path().is_dir())
        .flat_map(|dir| {
            std::fs::read_dir(dir.path())
                .into_iter()
                .flatten()
                .flatten()
                .map(|e| e.path())
                .filter(|p| p.extension().is_some_and(|x| x == "in"))
        })
        .collect();
    files.sort();
    files
}

#[test]
fn every_reference_input_round_trips_through_the_recipe() {
    let inputs = reference_inputs();
    assert!(inputs.len() >= 40);
    for path in &inputs {
        let imd = read_imd_file(path).unwrap();
        // Individual λ (LAMBDAS) and extended TI (PRECALCLAM) are not modelled by the recipe yet
        // (1.8); their reference inputs are refused, by design, until they are.
        let unmodelled: Vec<&str> = ["LAMBDAS", "PRECALCLAM"]
            .into_iter()
            .filter(|b| imd.raw_blocks.contains_key(*b))
            .collect();
        if !unmodelled.is_empty() {
            println!(
                "{}: skipped, the recipe does not model {} yet",
                path.display(),
                unmodelled.join(", ")
            );
            continue;
        }
        // GAMD and EDS are applied by the binary, not the recipe, and INNERLOOP selects a kernel:
        // the same passthrough policy the binary uses.
        let policy = PassthroughPolicy::allow(["GAMD", "EDS", "INNERLOOP"]);
        let (recipe, diag) = RunRecipe::from_imd_with(&imd, &policy)
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        assert_eq!(recipe.version, RECIPE_VERSION);
        let expected_passthrough = imd.raw_blocks.keys().any(|b| policy.allowed.contains(b));
        assert!(
            expected_passthrough || diag.notes.iter().all(|n| !n.contains("passed through")),
            "{}: unexpected passthrough {:?}",
            path.display(),
            diag.notes
        );

        // recipe → ImdParameters → recipe
        let back = recipe.to_imd();
        let (again, _) = RunRecipe::from_imd_with(&back, &policy).unwrap();
        assert_eq!(
            again,
            recipe,
            "{}: recipe → imd → recipe changed",
            path.display()
        );

        // recipe → .imd text → ImdParameters → recipe
        let text = recipe.to_imd_string(None);
        let reparsed =
            parse_imd_str(&text).unwrap_or_else(|e| panic!("{}: {e}\n{text}", path.display()));
        let (again, _) = RunRecipe::from_imd_with(&reparsed, &policy).unwrap();
        assert_eq!(
            again,
            recipe,
            "{}: recipe → text → recipe changed",
            path.display()
        );

        // TOML and JSON
        let toml_text = recipe.to_toml().unwrap();
        assert_eq!(
            RunRecipe::from_toml(&toml_text).unwrap(),
            recipe,
            "{}",
            path.display()
        );
        let json_text = recipe.to_json().unwrap();
        assert_eq!(
            RunRecipe::from_json(&json_text).unwrap(),
            recipe,
            "{}",
            path.display()
        );
    }
}

#[test]
fn default_recipe_is_what_gromosxx_does_with_every_optional_block_absent() {
    let d = RunRecipe::default();
    assert!(d.ensemble.thermostat.is_none(), "no MULTIBATH → no bath");
    assert!(
        d.ensemble.barostat.is_none(),
        "no PRESSURESCALE → no barostat"
    );
    assert_eq!(
        d.control.com_motion.every, 0,
        "no COMTRANSROT → no COM removal"
    );
    assert_eq!(d.perturbation.enabled, 0);
    assert_eq!(d.minimisation.method, 0);
    assert_eq!(d.forcefield.restraints.position.mode, 0);
    assert_eq!(d.forcefield.restraints.distance.mode, 0);
    // Derived from the single defaults table, not a second one.
    let (from_defaults, _) = RunRecipe::from_imd(&ImdParameters::default()).unwrap();
    assert_eq!(d, from_defaults);
}

#[test]
fn absent_blocks_are_reported_in_diagnostics() {
    let minimal = "SYSTEM\n  1 0\nEND\nSTEP\n  10 0.0 0.002\nEND\n";
    let imd = parse_imd_str(minimal).unwrap();
    let (_, diag) = RunRecipe::from_imd(&imd).unwrap();
    let joined = diag.notes.join("\n");
    assert!(
        joined.contains("MULTIBATH block absent: no temperature coupling"),
        "{joined}"
    );
    assert!(joined.contains("COMTRANSROT block absent"), "{joined}");
    assert!(joined.contains("FORCE block absent"), "{joined}");
}

#[test]
fn unmodelled_blocks_are_refused_unless_allowed() {
    let text = "SYSTEM\n  1 0\nEND\nGAMD\n  1 2 3\nEND\n";
    let imd = parse_imd_str(text).unwrap();
    let err = RunRecipe::from_imd(&imd).unwrap_err();
    assert!(
        matches!(err, RunError::UnknownBlock { ref block } if block == "GAMD"),
        "{err}"
    );

    let (recipe, diag) =
        RunRecipe::from_imd_with(&imd, &PassthroughPolicy::allow(["GAMD"])).unwrap();
    assert_eq!(recipe.passthrough.get("GAMD").map(|v| v.len()), Some(1));
    assert!(diag
        .notes
        .iter()
        .any(|n| n.contains("GAMD block passed through")));
    // …and it comes back out in the .imd text.
    assert!(recipe.to_imd_string(None).contains("GAMD\n1 2 3\nEND"));
}

#[test]
fn factories_build_the_same_runs_as_the_parameter_factories() {
    let nvt = RunRecipe::nvt(0.002, 100, 298.0, "hbonds").unwrap();
    assert_eq!(nvt.control.steps, 100);
    assert_eq!(
        nvt.ensemble.thermostat,
        Some(Thermostat {
            algorithm: ThermostatAlgorithm::Berendsen,
            baths: vec![Bath {
                temperature: 298.0,
                tau: 0.1
            }],
            dof_sets: vec![],
        })
    );
    assert_eq!(nvt.constraints.solute.ntc(), 2);
    let em = RunRecipe::minimize(50);
    assert!(em.is_minimization());
    assert!(RunRecipe::npt(0.002, 10, 300.0, 1.0, "none")
        .unwrap()
        .ensemble
        .barostat
        .is_some());
    assert!(RunRecipe::nve(0.002, 10, "nonsense").is_err());
}

#[test]
fn a_typo_in_a_recipe_field_is_an_error() {
    let mut json: serde_json::Value =
        serde_json::from_str(&RunRecipe::default().to_json().unwrap()).unwrap();
    json["control"]["stepz"] = serde_json::json!(5);
    let err = RunRecipe::from_json(&json.to_string()).unwrap_err();
    assert!(err.to_string().contains("stepz"), "{err}");
}
