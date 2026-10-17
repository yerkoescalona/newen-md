use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use newen_core::math::Vec3;
use newen_forces::nonbonded::{lj_crf_interaction, CRFParameters, ForceStorage, LJParameters};

fn crf_params() -> CRFParameters {
    let cutoff = 1.4_f64;
    CRFParameters {
        crf_cut: cutoff,
        crf_2cut3i: 1.0 / (2.0 * cutoff.powi(3)),
        crf_cut3i: 1.0 / cutoff,
        cutoff_sq: 1.4_f64.powi(2),
    }
}

fn bench_lj_crf_single(c: &mut Criterion) {
    let crf = crf_params();
    let r = Vec3::new(0.5, 0.0, 0.0);

    c.bench_function("lj_crf_interaction/single_pair", |b| {
        b.iter(|| {
            lj_crf_interaction(
                black_box(r),
                black_box(0.001),
                black_box(0.0001),
                black_box(0.5),
                &crf,
            )
        })
    });
}

fn bench_nonbonded_n_pairs(c: &mut Criterion) {
    let crf = crf_params();
    let lj = LJParameters {
        c6: 0.001,
        c12: 0.0001,
        cs6: 0.001,
        cs12: 0.0001,
    };
    let mut group = c.benchmark_group("nonbonded_n_pairs");

    for n in [100usize, 1_000, 10_000] {
        let pairs: Vec<Vec3> = (0..n)
            .map(|i| {
                let t = i as f64 * 0.05 + 0.3;
                Vec3::new(t, t * 0.5, t * 0.3)
            })
            .collect();

        group.bench_with_input(BenchmarkId::from_parameter(n), &pairs, |b, pairs| {
            b.iter(|| {
                pairs.iter().fold((0.0_f64, 0.0_f64), |acc, r| {
                    let (_, e_lj, e_crf) =
                        lj_crf_interaction(black_box(*r), lj.c6, lj.c12, 0.5, &crf);
                    (acc.0 + e_lj, acc.1 + e_crf)
                })
            })
        });
    }

    group.finish();
}

fn bench_force_storage_clear(c: &mut Criterion) {
    let mut group = c.benchmark_group("force_storage_clear");

    for n in [100usize, 1_000, 10_000] {
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            let mut storage = ForceStorage::new(n);
            b.iter(|| storage.clear())
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_lj_crf_single,
    bench_nonbonded_n_pairs,
    bench_force_storage_clear
);
criterion_main!(benches);
