//! `TicModel` against the full pipeline it stands in for, and its updates against a fresh build.

use std::time::Instant;

use sim_core::spec::*;
use sim_core::tic::TicModel;
use sim_core::Simulator;

/// The ProteoStar workbench fixture: a SimSpec fitted to a real top-down run.
fn fixture_spec() -> SimSpec {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../apps/desktop/e2e/data/02-18-20_jurkat_td_rep2_fract7.rt44-45.simspec.json"
    );
    SimSpec::read(path).expect("fixture spec")
}

/// Every scan's summed centroids, as the full pipeline produces them.
fn full_tic(spec: &SimSpec) -> Vec<f64> {
    Simulator::new(spec).unwrap().simulate_run().iter().map(|s| s.intensity.iter().sum()).collect()
}

fn model_tic(spec: &SimSpec) -> Vec<f64> {
    let t = TicModel::new(spec).unwrap().trace(None);
    t.signal.iter().zip(&t.noise).map(|(s, n)| s + n).collect()
}

/// Largest per-scan relative error, over scans holding at least 1% of the run's peak TIC.
fn worst_relative_error(model: &[f64], full: &[f64]) -> f64 {
    let peak = full.iter().copied().fold(0.0f64, f64::max);
    model
        .iter()
        .zip(full)
        .filter(|(_, &f)| f >= 0.01 * peak)
        .map(|(&m, &f)| (m / f - 1.0).abs())
        .fold(0.0, f64::max)
}

#[test]
fn matches_the_full_pipeline_on_a_fitted_spec() {
    let with_noise = fixture_spec();
    assert!(with_noise.noise.is_some());
    let mut clean = with_noise.clone();
    clean.noise = None;
    let mut no_jitter = with_noise.clone();
    no_jitter.noise.as_mut().unwrap().jitter.enabled = false;

    for (name, spec, tolerance) in [
        // Measured within ±2.3% on every scan. Without jitter the only approximations are overlap
        // between proteoforms (in a crowded scan neighbours' tails lift each other's maxima, which
        // the model does not count) and, without noise, the floor's brightest centroid.
        ("noise, no jitter", &no_jitter, 0.03),
        ("noise + jitter", &with_noise, 0.03),
        ("clean", &clean, 0.03),
    ] {
        let started = Instant::now();
        let model = model_tic(spec);
        let model_ms = started.elapsed().as_secs_f64() * 1e3;
        let full = full_tic(spec);
        let worst = worst_relative_error(&model, &full);
        eprintln!("{name}: worst scan off by {:.3}% (model built in {model_ms:.0} ms)", worst * 100.0);
        assert!(worst < tolerance, "{name}: {worst}");
    }
}

#[test]
fn an_edit_costs_little_and_matches_a_fresh_build() {
    let base = fixture_spec();
    let mut model = TicModel::new(&base).unwrap();

    let mut spec = base.clone();
    let edits: Vec<(&str, Box<dyn Fn(&mut SimSpec)>)> = vec![
        ("abundance", Box::new(|s| s.proteoforms[3].abundance *= 4.0)),
        ("rt", Box::new(|s| s.proteoforms[5].rt.mu += 0.3)),
        ("charge", Box::new(|s| {
            let SimChargeDistribution::Gaussian { sigma_z, .. } = s.proteoforms[7].charge;
            s.proteoforms[7].charge = SimChargeDistribution::Gaussian { mu_z: 20.5, sigma_z };
        })),
        ("remove", Box::new(|s| {
            s.proteoforms.remove(0);
        })),
        ("duplicate", Box::new(|s| {
            let p = s.proteoforms[4].clone();
            s.proteoforms.insert(5, p);
        })),
        ("noise level", Box::new(|s| s.noise.as_mut().unwrap().floor.noise_level_at_reference_mz *= 2.0)),
        ("noise off", Box::new(|s| s.noise = None)),
        ("peak width", Box::new(|s| s.peak_width_model = SimPeakWidthModel::Orbitrap { k: 2e-7 })),
    ];
    for (name, edit) in edits {
        edit(&mut spec);
        let started = Instant::now();
        model.update(&spec).unwrap();
        let updated = model.trace(Some(2));
        let ms = started.elapsed().as_secs_f64() * 1e3;
        eprintln!("{name}: update + trace {ms:.1} ms");
        assert_eq!(updated, TicModel::new(&spec).unwrap().trace(Some(2)), "after {name}");
    }
}

#[test]
fn the_selected_xic_is_part_of_the_tic() {
    let spec = fixture_spec();
    let t = TicModel::new(&spec).unwrap().trace(Some(0));
    let xic = t.selected.expect("proteoform 0 exists");
    assert!(!xic.target_mzs.is_empty());
    assert!(xic.target_mzs.windows(2).all(|w| w[0] < w[1]));
    assert!(xic.intensity.iter().zip(&t.signal).all(|(x, s)| *x <= *s));
    assert!(xic.intensity.iter().any(|&x| x > 0.0));
    // Out of range, or a background species: no XIC.
    assert!(TicModel::new(&spec).unwrap().trace(Some(spec.proteoforms.len())).selected.is_none());
}
