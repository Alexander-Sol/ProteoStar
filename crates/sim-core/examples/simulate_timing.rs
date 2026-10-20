//! Rough timings for the two workbench paths on a run-sized synthetic spec.
//!
//! `cargo run -p sim-core --release --example simulate_timing [proteoforms] [scans] [run-only]`
//! (defaults 1000 × 3000, full-density noise). `run-only` times just `simulate_run` in a fresh
//! process, unaffected by the gigabytes the exact `centroid_run` allocates and frees.

use std::time::Instant;

use sim_core::spec::*;
use sim_core::Simulator;

fn main() {
    let run_only = std::env::args().any(|a| a == "run-only");
    let args: Vec<usize> = std::env::args().skip(1).filter_map(|a| a.parse().ok()).collect();
    let n_proteoforms = args.first().copied().unwrap_or(1000);
    let n_scans = args.get(1).copied().unwrap_or(3000);

    // Deterministic spread of masses / RTs / charges; no RNG dependency needed.
    let frac = |i: usize, k: u64| ((i as u64 * k) % 1000) as f64 / 1000.0;
    let proteoforms = (0..n_proteoforms)
        .map(|i| SimProteoform {
            monoisotopic_mass: 5000.0 + 45_000.0 * frac(i, 7919),
            abundance: 1e8 * (1.0 + 20.0 * frac(i, 104_729)),
            rt: SimRtProfile { mu: 10.0 + 50.0 * frac(i, 1543), sigma: 0.2, tau: 0.1 },
            charge: SimChargeDistribution::Gaussian {
                mu_z: 8.0 + 30.0 * frac(i, 389),
                sigma_z: 2.5,
            },
            identifier: None,
        })
        .collect();

    let spec = SimSpec {
        schema_version: CURRENT_SCHEMA_VERSION,
        generator: None,
        source_run: None,
        instrument: SimInstrument::default(),
        charge_range: SimChargeRange { min: 3, max: 60 },
        peak_width_model: SimPeakWidthModel::Orbitrap { k: 1.2e-7 },
        scan_times: (0..n_scans).map(|s| 5.0 + 60.0 * s as f64 / n_scans as f64).collect(),
        reduction: None,
        noise: Some(SimNoise::orbitrap_defaults(0)),
        proteoforms,
        background: Vec::new(),
    };

    let t = Instant::now();
    drop(Simulator::new(&spec).unwrap());
    println!("Simulator::new        {:>8.1} ms  (first: includes the averagine parse)", t.elapsed().as_secs_f64() * 1e3);
    let t = Instant::now();
    let sim = Simulator::new(&spec).unwrap();
    println!("Simulator::new        {:>8.1} ms  (again: what each live edit pays)", t.elapsed().as_secs_f64() * 1e3);

    if run_only {
        let t = Instant::now();
        let run = sim.simulate_run();
        println!("simulate_run          {:>8.1} ms  ({} scans)", t.elapsed().as_secs_f64() * 1e3, run.len());
        return;
    }

    let mid = n_scans / 2;
    let t = Instant::now();
    let centroids = sim.centroid_scan(mid).unwrap();
    println!(
        "centroid_scan({mid})  {:>8.1} ms  ({} unreduced signal peaks)",
        t.elapsed().as_secs_f64() * 1e3,
        centroids.mz.len()
    );

    let t = Instant::now();
    let scan = sim.simulate_scan(mid).unwrap();
    println!(
        "simulate_scan({mid})  {:>8.1} ms  ({} peaks)",
        t.elapsed().as_secs_f64() * 1e3,
        scan.mz.len()
    );

    let t = Instant::now();
    let signal_run = sim.centroid_run();
    let signal_peaks: usize = signal_run.iter().map(|c| c.mz.len()).sum();
    println!(
        "centroid_run          {:>8.1} ms  ({signal_peaks} unreduced signal peaks)",
        t.elapsed().as_secs_f64() * 1e3
    );
    drop(signal_run);

    let t = Instant::now();
    let run = sim.simulate_run();
    let peaks: usize = run.iter().map(|s| s.mz.len()).sum();
    println!(
        "simulate_run          {:>8.1} ms  ({} scans, {peaks} peaks)",
        t.elapsed().as_secs_f64() * 1e3,
        run.len()
    );
    drop(run);

    session_timings(spec, mid);
}

/// What one live edit costs on a long-lived simulator: `update` to the edited spec, then the
/// viewed scan.
fn session_timings(mut spec: SimSpec, viewed: usize) {
    println!("\nsession: update + simulate_scan({viewed}) per edit");
    let mut sim = Simulator::new(&spec).unwrap();
    sim.simulate_scan(viewed).unwrap();

    let edit = |name: &str, sim: &mut Simulator, spec: &mut SimSpec, f: &dyn Fn(&mut SimSpec)| {
        f(spec);
        let t = Instant::now();
        sim.update(spec).unwrap();
        let updated = t.elapsed().as_secs_f64() * 1e3;
        let kept = sim.cached_scan_maxima();
        sim.simulate_scan(viewed).unwrap();
        println!(
            "  {name:<30} {:>8.1} ms  (update {updated:.1} ms; scan maxima kept {kept}/{})",
            t.elapsed().as_secs_f64() * 1e3,
            sim.scan_count()
        );
    };

    edit("abundance (noise on)", &mut sim, &mut spec, &|s| s.proteoforms[0].abundance *= 2.0);
    edit("RT mu (noise on)", &mut sim, &mut spec, &|s| s.proteoforms[1].rt.mu += 0.5);
    edit("mass (noise on)", &mut sim, &mut spec, &|s| s.proteoforms[2].monoisotopic_mass += 100.0);
    edit("noise level", &mut sim, &mut spec, &|s| {
        s.noise.as_mut().unwrap().floor.noise_level_at_reference_mz *= 1.05
    });
    edit("jitter off", &mut sim, &mut spec, &|s| s.noise.as_mut().unwrap().jitter.enabled = false);
    edit("noise off (cold scan maxima)", &mut sim, &mut spec, &|s| s.noise = None);
    edit("abundance (noise off)", &mut sim, &mut spec, &|s| s.proteoforms[0].abundance *= 2.0);
    edit("RT mu (noise off)", &mut sim, &mut spec, &|s| s.proteoforms[1].rt.mu += 0.5);
    edit("noise on", &mut sim, &mut spec, &|s| s.noise = Some(SimNoise::orbitrap_defaults(0)));
    edit("noise off again (warm)", &mut sim, &mut spec, &|s| s.noise = None);
    edit("peak width k (noise off)", &mut sim, &mut spec, &|s| {
        s.peak_width_model = SimPeakWidthModel::Orbitrap { k: 1.3e-7 }
    });
}
