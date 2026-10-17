use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use newen_core::math::{BoundaryCondition, Rectangular, Vacuum, Vec3};

fn bench_nearest_image(c: &mut Criterion) {
    let mut group = c.benchmark_group("nearest_image");

    let ri = Vec3::new(0.1, 0.2, 0.3);
    let rj = Vec3::new(4.9, 5.1, 0.4);

    group.bench_function("vacuum", |b| {
        let bc = Vacuum;
        b.iter(|| bc.nearest_image(black_box(ri), black_box(rj)))
    });

    group.bench_function("rectangular", |b| {
        let bc = Rectangular::new(Vec3::new(5.0, 5.0, 5.0));
        b.iter(|| bc.nearest_image(black_box(ri), black_box(rj)))
    });

    group.finish();
}

fn bench_put_into_box(c: &mut Criterion) {
    let bc = Rectangular::new(Vec3::new(5.0, 5.0, 5.0));
    let pos = Vec3::new(6.5, -0.5, 5.2);

    c.bench_function("put_into_box/rectangular", |b| {
        b.iter(|| bc.put_into_box(black_box(pos)))
    });
}

fn bench_vec3_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("vec3");

    let a = Vec3::new(1.0, 2.0, 3.0);
    let b = Vec3::new(4.0, 5.0, 6.0);

    group.bench_function("dot", |b_| b_.iter(|| black_box(a).dot(black_box(b))));

    group.bench_function("cross", |b_| b_.iter(|| black_box(a).cross(black_box(b))));

    group.bench_function("length", |b_| b_.iter(|| black_box(a).length()));

    group.bench_function("normalize", |b_| b_.iter(|| black_box(a).normalize()));

    group.finish();
}

fn bench_nearest_image_scaling(c: &mut Criterion) {
    let bc = Rectangular::new(Vec3::new(10.0, 10.0, 10.0));
    let mut group = c.benchmark_group("nearest_image_n_pairs");

    for n in [100usize, 1_000, 10_000] {
        let pairs: Vec<(Vec3, Vec3)> = (0..n)
            .map(|i| {
                let t = i as f64 * 0.1;
                (
                    Vec3::new(t.sin(), t.cos(), t * 0.01),
                    Vec3::new(t.cos(), t.sin(), t * 0.02),
                )
            })
            .collect();

        group.bench_with_input(BenchmarkId::from_parameter(n), &pairs, |b, pairs| {
            b.iter(|| {
                pairs
                    .iter()
                    .map(|(ri, rj)| bc.nearest_image(black_box(*ri), black_box(*rj)))
                    .fold(Vec3::ZERO, |acc, v| acc + v)
            })
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_nearest_image,
    bench_put_into_box,
    bench_vec3_ops,
    bench_nearest_image_scaling
);
criterion_main!(benches);
