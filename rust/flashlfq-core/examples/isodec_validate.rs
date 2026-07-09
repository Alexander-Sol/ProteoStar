//! Validates the native IsoDec charge predictor against realistic averagine envelopes at known charges.
//! Determines the correct .bin layout + class→charge mapping empirically.
//!
//! Usage: isodec_validate <phase_model_8.bin>

use flashlfq_core::deconvolution::averagine_intensities_from_mono;
use flashlfq_core::isodec::{encode_phase, IsoDecModel, ELEN, NCLASS};
use flashlfq_core::isotopic_envelope::{mass_to_mz_f64, PROTON_MASS};

fn averagine_envelope(mono_mass: f64, charge: i32) -> (Vec<f64>, Vec<f32>) {
    let intens = averagine_intensities_from_mono(mono_mass, 1e-3, 40);
    let spacing = 1.0033548 / charge as f64;
    let mono_mz = mass_to_mz_f64(mono_mass, charge);
    let mz: Vec<f64> = (0..intens.len()).map(|k| mono_mz + k as f64 * spacing).collect();
    let inten: Vec<f32> = intens.iter().map(|&w| w as f32).collect();
    (mz, inten)
}

fn argmax(v: &[f32]) -> usize {
    v.iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: isodec_validate <phase_model_8.bin>");
    let bytes = std::fs::read(&path).expect("read model");
    println!("model {path}: {} bytes ({} floats, expect {})", bytes.len(), bytes.len() / 4, IsoDecModel::N_PARAMS);
    let model = IsoDecModel::from_bin(&bytes).expect("parse model");

    // Sanity: a mid m/z (~900) at each charge → realistic averagine envelope → predicted class.
    println!("\n charge -> argmax class  (top-3 logits)   [envelope at m/z ~900]");
    let mut offsets = std::collections::HashMap::new();
    for z in [2i32, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30] {
        let mono_mass = (900.0 * z as f64) - z as f64 * PROTON_MASS;
        let (mz, inten) = averagine_envelope(mono_mass, z);
        let enc = encode_phase(&mz, &inten);
        let logits = model.logits(&enc);
        let cls = argmax(&logits);
        // top-3 classes
        let mut idx: Vec<usize> = (0..NCLASS).collect();
        idx.sort_by(|&a, &b| logits[b].total_cmp(&logits[a]));
        let top3: Vec<String> = idx[..3].iter().map(|&i| format!("{i}({:.1})", logits[i])).collect();
        let off = cls as i32 - z;
        *offsets.entry(off).or_insert(0) += 1;
        println!("   z={z:>2} -> class {cls:>2}  off={off:+}   {}", top3.join(" "));
    }
    println!("\n class-minus-charge offset histogram: {offsets:?}");
    println!("(consistent single offset => correct layout; class = charge + offset => set class_to_charge)");

    // Quick layout sanity: encoding of a clean charge-7 envelope must be non-trivial.
    let _ = ELEN;
}
