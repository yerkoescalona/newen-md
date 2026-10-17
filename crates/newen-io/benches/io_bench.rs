use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use newen_core::{configuration::Configuration, math::Vec3};
use newen_io::{
    energy::{EnergyFrame, EnergyWriter},
    force::ForceWriter,
    trajectory::TrajectoryWriter,
};

fn make_conf(n_atoms: usize) -> Configuration {
    let mut conf = Configuration::new(n_atoms, 1, 1);
    for i in 0..n_atoms {
        let t = i as f64 * 0.1;
        conf.current_mut().pos[i] = Vec3::new(t.sin(), t.cos(), t * 0.01);
    }
    conf
}

fn bench_trajectory_write(c: &mut Criterion) {
    let mut group = c.benchmark_group("trajectory_write_frame");

    for n in [10usize, 100, 1_000] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            b.iter_with_setup(
                || {
                    let path = std::env::temp_dir().join(format!("bench_traj_{n}.trc"));
                    let writer = TrajectoryWriter::new(&path, "bench", false, false).unwrap();
                    let conf = make_conf(n);
                    (writer, conf, path)
                },
                |(mut writer, conf, path)| {
                    writer
                        .write_frame(black_box(0), black_box(0.0), &conf)
                        .unwrap();
                    drop(writer);
                    std::fs::remove_file(path).ok();
                },
            )
        });
    }

    group.finish();
}

fn bench_energy_write(c: &mut Criterion) {
    c.bench_function("energy_write_frame", |b| {
        b.iter_with_setup(
            || {
                let path = std::env::temp_dir().join("bench_energy.tre");
                let writer = EnergyWriter::new(&path, "bench").unwrap();
                (writer, path)
            },
            |(mut writer, path)| {
                let frame = EnergyFrame::new(black_box(0.0), 100.0, -200.0, 300.0);
                writer.write_frame(black_box(&frame)).unwrap();
                drop(writer);
                std::fs::remove_file(path).ok();
            },
        )
    });
}

fn bench_force_write(c: &mut Criterion) {
    let mut group = c.benchmark_group("force_write_frame");

    for n in [10usize, 100, 1_000] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            let forces: Vec<Vec3> = (0..n)
                .map(|i| Vec3::new(i as f64, i as f64 * 0.5, i as f64 * 0.25))
                .collect();

            b.iter_with_setup(
                || {
                    let path = std::env::temp_dir().join(format!("bench_force_{n}.trf"));
                    let writer = ForceWriter::new(&path, "bench", false).unwrap();
                    (writer, path)
                },
                |(mut writer, path)| {
                    writer
                        .write_frame(black_box(0), black_box(0.0), &forces, None)
                        .unwrap();
                    drop(writer);
                    std::fs::remove_file(path).ok();
                },
            )
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_trajectory_write,
    bench_energy_write,
    bench_force_write
);
criterion_main!(benches);
