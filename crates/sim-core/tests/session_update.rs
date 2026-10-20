//! `Simulator::update` must be invisible in the output: after any sequence of edits, a
//! long-lived simulator produces exactly — bit for bit — what a fresh `Simulator::new` on the
//! final spec produces. Its caches (envelopes, per-scan maxima, the viewed scan's signal) may
//! only skip work.

use sim_core::spec::*;
use sim_core::Simulator;

fn base_spec(noise: Option<SimNoise>) -> SimSpec {
    // Spread over a long gradient so each proteoform's EMG covers only some of the scans.
    let frac = |i: usize, k: u64| ((i as u64 * k) % 1000) as f64 / 1000.0;
    let proteoform = |i: usize| SimProteoform {
        monoisotopic_mass: 5000.0 + 30_000.0 * frac(i, 7919),
        abundance: 1e8 * (1.0 + 20.0 * frac(i, 104_729)),
        rt: SimRtProfile { mu: 12.0 + 36.0 * frac(i, 1543), sigma: 0.2, tau: 0.1 },
        charge: SimChargeDistribution::Gaussian { mu_z: 8.0 + 25.0 * frac(i, 389), sigma_z: 2.5 },
        identifier: Some(format!("P{i}")),
    };
    SimSpec {
        schema_version: CURRENT_SCHEMA_VERSION,
        generator: None,
        source_run: None,
        instrument: SimInstrument::default(),
        charge_range: SimChargeRange { min: 3, max: 45 },
        peak_width_model: SimPeakWidthModel::Orbitrap { k: 1.2e-7 },
        scan_times: (0..60).map(|s| 10.0 + 40.0 * s as f64 / 60.0).collect(),
        reduction: None,
        noise,
        proteoforms: (0..60).map(proteoform).collect(),
        background: (60..70).map(proteoform).collect(),
    }
}

type Edit = (&'static str, fn(&mut SimSpec));

fn charge_mu(spec: &mut SimSpec, i: usize, mu: f64) {
    let SimChargeDistribution::Gaussian { sigma_z, .. } = spec.proteoforms[i].charge;
    spec.proteoforms[i].charge = SimChargeDistribution::Gaussian { mu_z: mu, sigma_z };
}

/// Applied in order, each on top of the last.
const EDITS: &[Edit] = &[
    ("abundance", |s| s.proteoforms[5].abundance *= 3.0),
    ("rt mu", |s| s.proteoforms[7].rt.mu += 1.5),
    ("rt sigma + tau", |s| {
        s.proteoforms[8].rt.sigma = 0.35;
        s.proteoforms[8].rt.tau = 0.0;
    }),
    ("charge mu", |s| charge_mu(s, 3, 20.0)),
    ("mass", |s| s.proteoforms[9].monoisotopic_mass += 1234.5),
    ("identifier only", |s| s.proteoforms[10].identifier = Some("renamed".into())),
    ("two scattered edits", |s| {
        s.proteoforms[2].abundance *= 0.5;
        s.proteoforms[40].rt.mu -= 2.0;
    }),
    ("append", |s| {
        let mut p = s.proteoforms[11].clone();
        p.monoisotopic_mass += 77.0;
        p.rt.mu = 30.0;
        s.proteoforms.push(p);
    }),
    ("insert in the middle", |s| {
        let mut p = s.proteoforms[12].clone();
        p.abundance *= 4.0;
        s.proteoforms.insert(20, p);
    }),
    ("remove first", |s| {
        s.proteoforms.remove(0);
    }),
    ("background", |s| s.background[3].abundance *= 10.0),
    ("dim the brightest", |s| {
        let brightest = (0..s.proteoforms.len())
            .max_by(|&a, &b| s.proteoforms[a].abundance.total_cmp(&s.proteoforms[b].abundance))
            .unwrap();
        s.proteoforms[brightest].abundance *= 1e-3;
    }),
    ("reduction", |s| {
        s.reduction = Some(SimReduction { relative_intensity_threshold: 3e-4, absolute_intensity_floor: 50.0 })
    }),
    ("noise on", |s| s.noise = Some(SimNoise::orbitrap_defaults(11))),
    ("noise level", |s| s.noise.as_mut().unwrap().floor.noise_level_at_reference_mz *= 1.1),
    ("noise level, across a power of two", |s| {
        s.noise.as_mut().unwrap().floor.noise_level_at_reference_mz *= 2.5
    }),
    ("jitter off", |s| s.noise.as_mut().unwrap().jitter.enabled = false),
    ("abundance under noise", |s| s.proteoforms[15].abundance *= 2.0),
    ("noise off", |s| s.noise = None),
    ("peak width", |s| s.peak_width_model = SimPeakWidthModel::Orbitrap { k: 1.5e-7 }),
    ("charge range", |s| s.charge_range = SimChargeRange { min: 5, max: 40 }),
    ("scan times", |s| {
        s.scan_times.pop();
        s.scan_times[0] -= 0.1;
    }),
    ("replace everything", |s| {
        let mut other = base_spec(None);
        for p in &mut other.proteoforms {
            p.monoisotopic_mass += 500.0;
        }
        s.proteoforms = other.proteoforms;
    }),
];

fn assert_same_scans(what: &str, a: &[flashlfq_core::peak_indexing::Scan], b: &[flashlfq_core::peak_indexing::Scan]) {
    assert_eq!(a.len(), b.len(), "{what}: scan count");
    for (s, (x, y)) in a.iter().zip(b).enumerate() {
        let bits = |v: &[f64]| v.iter().map(|f| f.to_bits()).collect::<Vec<_>>();
        assert_eq!(x.mz.len(), y.mz.len(), "{what}: scan {s} peak count");
        assert!(bits(&x.mz) == bits(&y.mz), "{what}: scan {s} m/z differs");
        assert!(bits(&x.intensity) == bits(&y.intensity), "{what}: scan {s} intensity differs");
        assert_eq!(x.retention_time.to_bits(), y.retention_time.to_bits(), "{what}: scan {s} RT");
    }
}

fn run_edits(start: SimSpec) {
    let mut spec = start;
    let mut live = Simulator::new(&spec).unwrap();
    // Warm every cache the edits could wrongly keep.
    live.global_max();
    live.simulate_scan(30).unwrap();

    for &(name, edit) in EDITS {
        edit(&mut spec);
        live.update(&spec).unwrap();
        let fresh = Simulator::new(&spec).unwrap();

        let viewed = 30.min(spec.scan_times.len() - 1);
        let picks = [0, viewed, spec.scan_times.len() - 1];
        let live_scans: Vec<_> = picks.iter().map(|&i| live.simulate_scan(i).unwrap()).collect();
        let fresh_scans: Vec<_> = picks.iter().map(|&i| fresh.simulate_scan(i).unwrap()).collect();
        assert_same_scans(&format!("{name} (simulate_scan)"), &live_scans, &fresh_scans);

        // The viewed scan again, now from the viewed-signal cache.
        assert_same_scans(
            &format!("{name} (cached viewed scan)"),
            &[live.simulate_scan(viewed).unwrap()],
            &[fresh.simulate_scan(viewed).unwrap()],
        );
        assert_same_scans(&format!("{name} (simulate_run)"), &live.simulate_run(), &fresh.simulate_run());
        assert_eq!(live.global_max().to_bits(), fresh.global_max().to_bits(), "{name}: global max");
    }
}

#[test]
fn update_matches_new_starting_without_noise() {
    run_edits(base_spec(None));
}

#[test]
fn update_matches_new_starting_with_noise() {
    run_edits(base_spec(Some(SimNoise::orbitrap_defaults(3))));
}

#[test]
fn a_local_edit_keeps_most_scan_maxima() {
    let mut spec = base_spec(None);
    let mut live = Simulator::new(&spec).unwrap();
    live.global_max();
    let n = live.scan_count();
    assert_eq!(live.cached_scan_maxima(), n);

    spec.proteoforms[5].abundance *= 3.0;
    live.update(&spec).unwrap();
    let kept = live.cached_scan_maxima();
    assert!(kept > n / 2 && kept < n, "kept {kept} of {n} scan maxima after one abundance edit");

    // Noise and reduction edits keep all of them.
    live.global_max();
    spec.noise = Some(SimNoise::orbitrap_defaults(0));
    spec.reduction = Some(SimReduction { relative_intensity_threshold: 1e-3, absolute_intensity_floor: 0.0 });
    live.update(&spec).unwrap();
    spec.noise = None;
    live.update(&spec).unwrap();
    assert_eq!(live.cached_scan_maxima(), n);

    // A layout edit keeps none.
    spec.peak_width_model = SimPeakWidthModel::Orbitrap { k: 1.3e-7 };
    live.update(&spec).unwrap();
    assert_eq!(live.cached_scan_maxima(), 0);
}

#[test]
fn a_failed_update_changes_nothing() {
    let spec = base_spec(None);
    let mut live = Simulator::new(&spec).unwrap();
    let before = live.simulate_run();

    let mut bad = spec.clone();
    bad.proteoforms[0].rt.sigma = -1.0;
    assert!(live.update(&bad).is_err());
    assert_eq!(live.spec(), &spec);
    assert_same_scans("after a failed update", &live.simulate_run(), &before);
}
