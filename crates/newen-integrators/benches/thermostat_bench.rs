use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use newen_core::{
    configuration::Configuration,
    math::Vec3,
    topology::{Atom, Topology},
};
use newen_integrators::thermostats::{
    berendsen_thermostat, nose_hoover_thermostat, BerendsenThermostatParameters,
    NoseHooverThermostatParameters,
};

fn make_system(n_atoms: usize) -> (Topology, Configuration) {
    let mut topo = Topology::new();
    for i in 0..n_atoms {
        topo.moltypes[0].atoms.push(Atom {
            name: format!("C{i}"),
            residue_nr: 1,
            residue_name: "BNZ".to_string(),
            iac: 0,
            mass: 12.0,
            charge: 0.0,
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: false,
        });
    }
    topo.mass = vec![12.0_f64; n_atoms];

    let mut conf = Configuration::new(n_atoms, 1, 1);
    for i in 0..n_atoms {
        let t = i as f64 * 0.1;
        conf.current_mut().vel[i] = Vec3::new(t.sin() * 0.1, t.cos() * 0.1, t * 0.01);
    }

    (topo, conf)
}

fn bench_berendsen(c: &mut Criterion) {
    let params = BerendsenThermostatParameters::default();
    let mut group = c.benchmark_group("berendsen_thermostat");

    for n in [100usize, 1_000, 10_000] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            let (topo, mut conf) = make_system(n);
            b.iter(|| berendsen_thermostat(&topo, black_box(&mut conf), 0.002, &params))
        });
    }

    group.finish();
}

fn bench_nose_hoover(c: &mut Criterion) {
    let mut params = NoseHooverThermostatParameters::default();
    let mut group = c.benchmark_group("nose_hoover_thermostat");

    for n in [100usize, 1_000, 10_000] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            let (topo, mut conf) = make_system(n);
            b.iter(|| {
                nose_hoover_thermostat(&topo, black_box(&mut conf), 0.002, black_box(&mut params))
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_berendsen, bench_nose_hoover);
criterion_main!(benches);
