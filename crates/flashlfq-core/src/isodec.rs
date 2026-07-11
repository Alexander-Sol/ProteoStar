//! Native Rust port of **IsoDec** charge-state assignment (Marty lab, JACS 2025) — the neural-network
//! deconvoluter that produced one of the top-down ground-truth sets. See
//! `agent_info/IsoDec-Port-Investigation.md`.
//!
//! IsoDec turns a centroid peak cluster into a fixed **phase encoding** and classifies its charge with a
//! tiny single-hidden-layer MLP. This module ports the deterministic pieces (encoding + forward pass);
//! the trained weights are loaded from the shipped `phase_model_8.bin` (see [`IsoDecModel::from_bin`]).
//!
//! Reference encoding (`unidec/IsoDec/encoding.py::encode_phase`, phaseres=8):
//! ```text
//! rescale[j] = mz[j] / mass_diff_c
//! for z_index in 0..maxz:               # charge = z_index + 1
//!     for each peak j:
//!         phase = (rescale[j] * (z_index+1)) mod 1
//!         bin   = floor(phase * phaseres)          # 0..phaseres-1
//!         phases[z_index][bin] += intensity[j]
//! phases /= max(phases)                            # normalise to peak 1.0
//! ```
//! Network (`models.py::Fast8PhaseNeuralNetwork`): `Linear(400→400) + ReLU + Linear(400→50)`, argmax
//! over the 50 logits → charge. Input is the row-major flattened `(maxz=50, phaseres=8)` matrix.

/// ¹³C−¹²C mass difference used for the phase (IsoDec `mass_diff_c`).
pub const MASS_DIFF_C: f64 = 1.0033548;
/// Maximum charge considered by the default phase model.
pub const MAXZ: usize = 50;
/// Phase resolution (bins) of the default (Fast8) model.
pub const PHASERES: usize = 8;
/// Flattened encoding length / MLP input width.
pub const ELEN: usize = MAXZ * PHASERES; // 400
/// Hidden width of the Fast8 MLP.
pub const HIDDEN: usize = 400;
/// Number of output classes (charge logits).
pub const NCLASS: usize = 50;

/// Encodes a centroid cluster into the row-major `(MAXZ × PHASERES)` phase histogram, normalised so its
/// maximum bin is 1.0. `mz`/`intensity` are parallel; peaks may be in any order. Ports
/// `encode_phase(..., phaseres=8)` verbatim.
pub fn encode_phase(mz: &[f64], intensity: &[f32]) -> [f32; ELEN] {
    let mut phases = [0f32; ELEN];
    if mz.is_empty() {
        return phases;
    }
    for (j, &m) in mz.iter().enumerate() {
        let rescale = m / MASS_DIFF_C;
        let inten = intensity[j];
        for zi in 0..MAXZ {
            let phase = (rescale * (zi as f64 + 1.0)).rem_euclid(1.0);
            let mut bin = (phase * PHASERES as f64).floor() as usize;
            if bin >= PHASERES {
                bin = PHASERES - 1; // guard the phase==~1.0 edge
            }
            phases[zi * PHASERES + bin] += inten;
        }
    }
    let maxv = phases.iter().cloned().fold(0f32, f32::max);
    if maxv > 0.0 {
        for p in phases.iter_mut() {
            *p /= maxv;
        }
    }
    phases
}

/// The Fast8 phase MLP weights: `Linear(ELEN→HIDDEN) + ReLU + Linear(HIDDEN→NCLASS)`.
/// PyTorch `Linear` computes `y = x·Wᵀ + b`, so `w1` is `[HIDDEN][ELEN]` and `w2` is `[NCLASS][HIDDEN]`.
pub struct IsoDecModel {
    pub w1: Vec<f32>, // HIDDEN * ELEN, row-major [out][in]
    pub b1: Vec<f32>, // HIDDEN
    pub w2: Vec<f32>, // NCLASS * HIDDEN
    pub b2: Vec<f32>, // NCLASS
}

impl IsoDecModel {
    /// Total number of f32 parameters in the Fast8 model.
    pub const N_PARAMS: usize = HIDDEN * ELEN + HIDDEN + NCLASS * HIDDEN + NCLASS; // 180450

    /// Loads the model from a raw little-endian f32 blob laid out `[w1, b1, w2, b2]`. A leading header
    /// (any bytes before the last `N_PARAMS` floats) is skipped. Returns `None` if the file is too small
    /// or its float count doesn't end on the expected layout.
    pub fn from_bin(bytes: &[u8]) -> Option<IsoDecModel> {
        let need = Self::N_PARAMS * 4;
        if bytes.len() < need {
            return None;
        }
        // Take the LAST N_PARAMS floats (skip any header prefix).
        let start = bytes.len() - need;
        let floats: Vec<f32> = bytes[start..]
            .chunks_exact(4)
            .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();
        let mut o = 0;
        let take = |o: &mut usize, n: usize| {
            let s = floats[*o..*o + n].to_vec();
            *o += n;
            s
        };
        let w1 = take(&mut o, HIDDEN * ELEN);
        let b1 = take(&mut o, HIDDEN);
        let w2 = take(&mut o, NCLASS * HIDDEN);
        let b2 = take(&mut o, NCLASS);
        Some(IsoDecModel { w1, b1, w2, b2 })
    }

    /// Runs the forward pass on a flattened phase encoding and returns the 50 output logits.
    pub fn logits(&self, x: &[f32; ELEN]) -> [f32; NCLASS] {
        // Hidden = ReLU(W1 x + b1)
        let mut h = [0f32; HIDDEN];
        for (o, hv) in h.iter_mut().enumerate() {
            let mut acc = self.b1[o];
            let row = &self.w1[o * ELEN..(o + 1) * ELEN];
            for (wi, xv) in row.iter().zip(x.iter()) {
                acc += wi * xv;
            }
            *hv = acc.max(0.0);
        }
        // Out = W2 h + b2
        let mut out = [0f32; NCLASS];
        for (o, ov) in out.iter_mut().enumerate() {
            let mut acc = self.b2[o];
            let row = &self.w2[o * HIDDEN..(o + 1) * HIDDEN];
            for (wi, hv) in row.iter().zip(h.iter()) {
                acc += wi * hv;
            }
            *ov = acc;
        }
        out
    }

    /// Predicts the charge of a centroid cluster: encode → forward → argmax. The mapping from output
    /// class index to charge is validated empirically (see tests); exposed as [`class_to_charge`].
    pub fn predict_charge(&self, mz: &[f64], intensity: &[f32]) -> i32 {
        let enc = encode_phase(mz, intensity);
        let logits = self.logits(&enc);
        let cls = argmax(&logits);
        class_to_charge(cls)
    }
}

/// Maps an output class index (0..NCLASS-1) to a charge. **Validated as the identity** against
/// realistic averagine envelopes at charges 2..30 (every one predicted its own class exactly, offset 0;
/// runner-up logits are the z↔2z harmonics). Class 0 is the no-call/charge-0 sentinel.
#[inline]
pub fn class_to_charge(cls: usize) -> i32 {
    cls as i32
}

/// The Fast8 phase model, embedded from `models/phase_model_8.bin` (raw `[w1,b1,w2,b2]` f32 blob shipped
/// in the mzLib/UniDec distribution). Built once; use for charge prediction without an external file.
pub fn default_model() -> IsoDecModel {
    static BYTES: &[u8] = include_bytes!("../models/phase_model_8.bin");
    IsoDecModel::from_bin(BYTES).expect("embedded phase_model_8.bin is the expected Fast8 layout")
}

#[inline]
fn argmax(v: &[f32]) -> usize {
    let mut bi = 0;
    let mut bv = f32::NEG_INFINITY;
    for (i, &x) in v.iter().enumerate() {
        if x > bv {
            bv = x;
            bi = i;
        }
    }
    bi
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A clean synthetic isotope envelope at a given charge encodes so that charge z's phase row
    /// concentrates all intensity in bin 0 (peaks sit exactly on the 1.0033/z grid → phase ≈ 0).
    fn synth(mono_mz: f64, charge: i32, n: usize) -> (Vec<f64>, Vec<f32>) {
        let spacing = MASS_DIFF_C / charge as f64;
        let mz: Vec<f64> = (0..n).map(|k| mono_mz + k as f64 * spacing).collect();
        let inten: Vec<f32> = (0..n).map(|k| (1.0 / (1.0 + k as f32))).collect();
        (mz, inten)
    }

    #[test]
    fn encoding_concentrates_true_charge_in_one_bin() {
        let z = 7usize;
        let (mz, inten) = synth(900.0, z as i32, 6);
        let enc = encode_phase(&mz, &inten);
        // A clean charge-z envelope's peaks all share ONE phase bin (phase is constant mod 1 across the
        // 1.0033/z-spaced teeth) — so the true-charge row's max bin holds essentially the whole row.
        let row: Vec<f32> = (0..PHASERES).map(|b| enc[(z - 1) * PHASERES + b]).collect();
        let rmax = row.iter().cloned().fold(0f32, f32::max);
        let rsum: f32 = row.iter().sum();
        assert!(rmax >= 0.99 * rsum, "true-charge row should sit in one bin: max {rmax} of sum {rsum}");
    }

    #[test]
    fn encoding_is_normalised() {
        let (mz, inten) = synth(900.0, 5, 5);
        let enc = encode_phase(&mz, &inten);
        let maxv = enc.iter().cloned().fold(0f32, f32::max);
        assert!((maxv - 1.0).abs() < 1e-6);
    }

    #[test]
    fn embedded_model_predicts_known_charges() {
        use crate::deconvolution::averagine_intensities_from_mono;
        use crate::isotopic_envelope::{mass_to_mz_f64, PROTON_MASS};
        let model = default_model();
        // Realistic averagine envelope at a fixed m/z (~900) for a range of charges → predicts each.
        for z in [2i32, 4, 6, 8, 10, 15, 20, 30] {
            let mono_mass = 900.0 * z as f64 - z as f64 * PROTON_MASS;
            let intens = averagine_intensities_from_mono(mono_mass, 1e-3, 40);
            let spacing = MASS_DIFF_C / z as f64;
            let mono_mz = mass_to_mz_f64(mono_mass, z);
            let mz: Vec<f64> = (0..intens.len()).map(|k| mono_mz + k as f64 * spacing).collect();
            let inten: Vec<f32> = intens.iter().map(|&w| w as f32).collect();
            assert_eq!(model.predict_charge(&mz, &inten), z, "IsoDec should recover charge {z}");
        }
    }
}
