//! Dump MS2 precursor / isolation-window metadata for MS2-linkability analysis.
//!
//! Reads a spectra file (mzML or Thermo `.raw`) and, for every MS2 spectrum, writes one TSV row:
//!   rt_min      retention time of the MS2 scan (minutes)
//!   iso_target  isolation-window target m/z (the recorded center)
//!   iso_lo      lower m/z bound of the isolation window
//!   iso_hi      upper m/z bound of the isolation window
//!   sel_mz      selected-ion (precursor) m/z
//!   charge      precursor charge if the instrument assigned one, else 0
//!
//! Only precursor metadata is decoded (DetailLevel::MetadataOnly) — the MS2 peak arrays are never
//! read, so this is fast even on the big files. `.raw` reading needs a .NET 8 runtime.
//!
//! Usage:
//!   cargo run --release --example dump_ms2_precursors -- <spectra_file> <out.tsv>

use std::fs::File;
use std::io::{BufWriter, Write};

use mzdata::io::{DetailLevel, MZReader, ThermoRawReader};
use mzdata::prelude::*;
use mzdata::spectrum::Spectrum;

fn write_ms2<I: Iterator<Item = Spectrum>>(reader: I, out: &mut BufWriter<File>) -> usize {
    let mut n = 0usize;
    for spectrum in reader {
        if spectrum.ms_level() != 2 {
            continue;
        }
        let rt = spectrum.start_time();
        let Some(prec) = spectrum.precursor() else {
            continue;
        };
        let iso = &prec.isolation_window;
        let (sel_mz, charge) = prec
            .ions
            .first()
            .map(|ion| (ion.mz, ion.charge.unwrap_or(0)))
            .unwrap_or((iso.target as f64, 0));
        writeln!(
            out,
            "{:.5}\t{:.5}\t{:.5}\t{:.5}\t{:.5}\t{}",
            rt, iso.target, iso.lower_bound, iso.upper_bound, sel_mz, charge
        )
        .expect("write row");
        n += 1;
    }
    n
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: dump_ms2_precursors <spectra_file> <out.tsv>");
        std::process::exit(1);
    }
    let path = &args[1];
    let out_path = &args[2];

    let mut out = BufWriter::new(File::create(out_path).expect("create out tsv"));
    writeln!(out, "rt_min\tiso_target\tiso_lo\tiso_hi\tsel_mz\tcharge").expect("header");

    eprintln!("reading MS2 precursors from {path} ...");
    let is_thermo_raw = std::path::Path::new(path)
        .extension()
        .map(|e| e.eq_ignore_ascii_case("raw"))
        .unwrap_or(false);

    let n = if is_thermo_raw {
        let reader = ThermoRawReader::new_with_detail_level_and_centroiding(
            path,
            DetailLevel::MetadataOnly,
            false,
        )
        .expect("open thermo raw");
        write_ms2(reader, &mut out)
    } else {
        let mut reader = MZReader::open_path(path).expect("open spectra file");
        reader.set_detail_level(DetailLevel::MetadataOnly);
        write_ms2(reader, &mut out)
    };

    out.flush().ok();
    eprintln!("wrote {n} MS2 precursor rows to {out_path}");
}
