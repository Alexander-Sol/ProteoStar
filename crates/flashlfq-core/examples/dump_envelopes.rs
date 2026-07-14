//! Dump the isotope-envelope comb weights of the target averagine and each decoy model at a given
//! most-abundant (mode) mass, for a visual comparison plot. Uses `EnvelopeModel::comb_weights`, which is
//! keyed by the **most-intense** mass and normalized to max 1.0 — so every model's mode sits at the same
//! mass and the same height, and only the envelope *shape* differs.
//!
//! Usage: dump_envelopes <mode_mass> <out.tsv>
//! Emits `model<TAB>k<TAB>weight` where k is the ¹³C index from the model's monoisotope (the mode index
//! is argmax of the weights; the plotter aligns modes and lays teeth at the model's own spacing).

use std::io::Write;
use flashlfq_core::deconvolution::EnvelopeModel;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: dump_envelopes <mode_mass> <out.tsv>");
        std::process::exit(2);
    }
    let mode_mass: f64 = args[1].parse().expect("mode_mass must be a float");
    let out_path = &args[2];

    // (label, model). "shifted" is not a distinct envelope model — it is the real averagine weights laid
    // on a 0.94-Da lattice — so the plotter derives it from the averagine row; it is not dumped here.
    let models: &[(&str, EnvelopeModel)] = &[
        ("averagine", EnvelopeModel::Averagine),
        ("weird", EnvelopeModel::CustomDecoy),
        ("shuffled", EnvelopeModel::ShuffledDecoy),
    ];

    let mut out = std::fs::File::create(out_path).expect("create out");
    writeln!(out, "model\tk\tweight").unwrap();
    for (label, model) in models {
        // Long template so the full 12 kDa envelope is captured (mono-to-mode teeth are all retained).
        let w = model.comb_weights(mode_mass, 1e-4, 200);
        for (k, wk) in w.iter().enumerate() {
            writeln!(out, "{label}\t{k}\t{wk:.6}").unwrap();
        }
        eprintln!("{label}: {} teeth, mode at index {}", w.len(),
                  w.iter().enumerate().fold((0usize, 0.0f64), |(bi, bv), (i, &v)| if v > bv { (i, v) } else { (bi, bv) }).0);
    }
    eprintln!("wrote {out_path}");
}
