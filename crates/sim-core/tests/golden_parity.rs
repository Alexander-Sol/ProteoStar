//! Parity against C# goldens in `parity/golden/`, each a `<case>.simspec.json` plus the
//! `<case>.scans.tsv` C# produced from it.
//!
//! - `golden-v1` — the Phase-0 data-free golden (contract §5): `Simulator.SimulateCentroided`,
//!   i.e. forward model + centroider only. Written by
//!   `mzLib/Test/TopDownSimulator/SimSpecGoldenGenerator.cs`. Rust's unpruned centroider
//!   reproduces it **bit for bit**.
//! - `reduced-v1`, `noise-jitter-v1`, `noise-nojitter-v1` — the same models through
//!   `Simulator.PrepareMs1` (reduction, then noise floor ± jitter). Written by a throwaway console
//!   app against mzLib `TopDownSimClean` @ 733955c0 + the uncommitted Phase-0 `SimulationSpec.cs`;
//!   to regenerate, call `PrepareMs1` on each spec exactly as `Simulator.WriteMzml` would, except
//!   that a disabled jitter is passed as `PeakJitterModel.None` (C#'s `jitter ?? Orbitrap` would
//!   otherwise turn "disabled" into the default jitter). Rust's pruned centroiding matches
//!   `reduced-v1` within the contract tolerance; its native noise streams match the noise goldens
//!   statistically.

use std::path::PathBuf;

use sim_core::spec::SimSpec;
use sim_core::Simulator;

/// Contract §5 tolerances.
const MZ_TOLERANCE_TH: f64 = 1e-6;
const INTENSITY_REL_TOLERANCE: f64 = 1e-4;

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("parity/golden").join(name)
}

fn read_spec(case: &str) -> SimSpec {
    SimSpec::read(golden(&format!("{case}.simspec.json"))).unwrap()
}

type Spectrum = (Vec<f64>, Vec<f64>);

/// Per scan, the (m/z, intensity) arrays of `<case>.scans.tsv`.
fn read_expected(case: &str, spec: &SimSpec) -> Vec<Spectrum> {
    let text = std::fs::read_to_string(golden(&format!("{case}.scans.tsv"))).unwrap();
    let mut scans: Vec<Spectrum> = vec![(Vec::new(), Vec::new()); spec.scan_times.len()];
    for line in text.lines().skip(1).filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let s: usize = f[0].parse().unwrap();
        assert_eq!(f[1].parse::<usize>().unwrap(), s + 1);
        assert_eq!(f[2].parse::<f64>().unwrap(), spec.scan_times[s], "{case} scan {s}: retention time");
        scans[s].0.push(f[3].parse().unwrap());
        scans[s].1.push(f[4].parse().unwrap());
    }
    scans
}

/// Asserts `actual` matches the case's golden peak for peak, within contract tolerance, and
/// reports how much of it is bit-exact.
fn assert_matches_golden(case: &str, spec: &SimSpec, actual: &[Spectrum]) {
    let expected = read_expected(case, spec);
    assert_eq!(actual.len(), expected.len(), "{case}: scan count");

    let (mut peaks, mut bit_exact) = (0usize, 0usize);
    let (mut max_dmz, mut max_rel) = (0.0f64, 0.0f64);
    for (s, ((want_mz, want_y), (mz, y))) in expected.iter().zip(actual).enumerate() {
        assert_eq!(mz.len(), want_mz.len(), "{case} scan {s}: peak count");
        for i in 0..want_mz.len() {
            let dmz = (mz[i] - want_mz[i]).abs();
            let rel = (y[i] - want_y[i]).abs() / want_y[i].abs().max(f64::MIN_POSITIVE);
            max_dmz = max_dmz.max(dmz);
            max_rel = max_rel.max(rel);
            bit_exact += usize::from(mz[i] == want_mz[i] && y[i] == want_y[i]);
            assert!(dmz <= MZ_TOLERANCE_TH, "{case} scan {s} peak {i}: m/z {} vs {}", mz[i], want_mz[i]);
            assert!(
                rel <= INTENSITY_REL_TOLERANCE,
                "{case} scan {s} peak {i} (m/z {}): intensity {} vs {}",
                want_mz[i],
                y[i],
                want_y[i]
            );
        }
        peaks += want_mz.len();
    }

    println!(
        "{case}: {peaks} peaks over {} scans; {bit_exact} bit-exact; max |Δm/z| {max_dmz:e} Th, \
         max rel Δintensity {max_rel:e}",
        actual.len()
    );
}

/// `simulate_run`, after checking that the live path (`simulate_scan`, on a fresh instance so it
/// has to find any run-level state itself) produces the same scans.
fn simulate_both_paths(case: &str, spec: &SimSpec) -> Vec<Spectrum> {
    let run = Simulator::new(spec).unwrap().simulate_run();
    let live = Simulator::new(spec).unwrap();
    for (s, from_run) in run.iter().enumerate() {
        let scan = live.simulate_scan(s).unwrap();
        assert_eq!(scan.one_based_scan_number, s as i32 + 1);
        assert_eq!(scan.retention_time, spec.scan_times[s]);
        assert_eq!(
            (&scan.mz, &scan.intensity),
            (&from_run.mz, &from_run.intensity),
            "{case}: simulate_scan({s}) != simulate_run"
        );
    }
    run.into_iter().map(|s| (s.mz, s.intensity)).collect()
}

#[test]
fn golden_v1_centroids_match_csharp_bit_for_bit() {
    let spec = read_spec("golden-v1");
    let actual: Vec<Spectrum> = Simulator::new(&spec)
        .unwrap()
        .centroid_run()
        .into_iter()
        .map(|c| (c.mz, c.intensity))
        .collect();
    assert_matches_golden("golden-v1", &spec, &actual);

    let expected = read_expected("golden-v1", &spec);
    assert_eq!(actual, expected, "unpruned centroiding should be bit-exact");
}

#[test]
fn reduced_matches_csharp() {
    let spec = read_spec("reduced-v1");
    let actual = simulate_both_paths("reduced-v1", &spec);
    assert_matches_golden("reduced-v1", &spec, &actual);
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    sorted[((sorted.len() - 1) as f64 * q).round() as usize]
}

/// Native noise streams can't match C# draw for draw, so compare what a feature finder sees:
/// how many peaks there are, and the distributions of their intensities and positions.
fn noise_matches_csharp_statistically(case: &str) {
    let spec = read_spec(case);
    let rust = simulate_both_paths(case, &spec);
    let csharp = read_expected(case, &spec);

    let pooled = |scans: &[Spectrum], pick: fn(&Spectrum) -> &Vec<f64>| {
        let mut v: Vec<f64> = scans.iter().flat_map(|s| pick(s).iter().copied()).collect();
        v.sort_by(f64::total_cmp);
        v
    };
    let (rust_mz, csharp_mz) = (pooled(&rust, |s| &s.0), pooled(&csharp, |s| &s.0));
    let (rust_y, csharp_y) = (pooled(&rust, |s| &s.1), pooled(&csharp, |s| &s.1));

    let count_ratio = rust_mz.len() as f64 / csharp_mz.len() as f64;
    println!("{case}: {} peaks (C# {}), ratio {count_ratio:.4}", rust_mz.len(), csharp_mz.len());
    assert!((count_ratio - 1.0).abs() < 0.03, "{case}: peak count ratio {count_ratio}");

    for q in [0.1, 0.25, 0.5, 0.75, 0.9] {
        let y_ratio = quantile(&rust_y, q) / quantile(&csharp_y, q);
        let mz_ratio = quantile(&rust_mz, q) / quantile(&csharp_mz, q);
        println!("  q{q}: intensity ratio {y_ratio:.4}, m/z ratio {mz_ratio:.5}");
        assert!((y_ratio - 1.0).abs() < 0.05, "{case}: intensity q{q} ratio {y_ratio}");
        assert!((mz_ratio - 1.0).abs() < 0.01, "{case}: m/z q{q} ratio {mz_ratio}");
    }
}

#[test]
fn noise_with_jitter_matches_csharp_statistically() {
    noise_matches_csharp_statistically("noise-jitter-v1");
}

#[test]
fn noise_without_jitter_matches_csharp_statistically() {
    noise_matches_csharp_statistically("noise-nojitter-v1");
}

#[test]
fn same_seed_reproduces_and_different_seed_differs() {
    let mut spec = read_spec("noise-jitter-v1");
    let a = Simulator::new(&spec).unwrap().simulate_scan(10).unwrap();
    let b = Simulator::new(&spec).unwrap().simulate_scan(10).unwrap();
    assert_eq!((&a.mz, &a.intensity), (&b.mz, &b.intensity));

    spec.noise.as_mut().unwrap().seed += 1;
    let c = Simulator::new(&spec).unwrap().simulate_scan(10).unwrap();
    assert_ne!(a.mz, c.mz);
}

#[test]
fn spec_round_trips_through_json() {
    for case in ["golden-v1", "noise-jitter-v1"] {
        let spec = read_spec(case);
        let again = SimSpec::from_json(&spec.to_json().unwrap()).unwrap();
        assert_eq!(spec, again);
    }
}

#[test]
fn unknown_schema_version_is_refused() {
    let text = std::fs::read_to_string(golden("golden-v1.simspec.json")).unwrap();
    let bumped = text.replacen("\"schemaVersion\": 1", "\"schemaVersion\": 2", 1);
    assert!(matches!(
        SimSpec::from_json(&bumped),
        Err(sim_core::SimError::UnsupportedSchema(2))
    ));
}

#[test]
fn out_of_range_scan_is_an_error() {
    let spec = read_spec("golden-v1");
    assert!(matches!(
        Simulator::new(&spec).unwrap().simulate_scan(spec.scan_times.len()),
        Err(sim_core::SimError::ScanOutOfRange { .. })
    ));
}
