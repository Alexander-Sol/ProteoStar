//! Runs the full Phase-1 MS2 engine + Phase-3 match-between-runs transfer over the real K562
//! corpus and writes the resulting MBR **feature table** to the Parquet path given as the first
//! CLI arg. Used by the P3.2e Python acceptance check (read it back with pyarrow). See PLAN.md
//! P3.2e.
//!
//! Usage: `write_mbr_feature_table <output.parquet>` — resolves the corpus relative to the
//! crate manifest (`mzLib/Test/FlashLFQ/TestData`), so no data paths need to be passed.

use std::collections::HashMap;
use std::path::PathBuf;

use flashlfq_core::engine::run_msms;
use flashlfq_core::mbr_search::run_mbr;
use flashlfq_core::parquet_output::write_feature_table_parquet;
use flashlfq_core::psm_tsv::read_identifications;

const FILE_3: &str = "20100614_Velos1_TaGe_SA_K562_3";
const FILE_4: &str = "20100614_Velos1_TaGe_SA_K562_4";

fn test_data(relative: &str) -> PathBuf {
    flashlfq_core::mzlib_test_data(relative)
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: write_mbr_feature_table <output.parquet>");

    let ids = read_identifications(test_data("AllPSMs.psmtsv"))
        .expect("AllPSMs.psmtsv parses into identifications");
    let mut file_to_mzml = HashMap::new();
    file_to_mzml.insert(FILE_3.to_string(), test_data(&format!("{FILE_3}.mzML")));
    file_to_mzml.insert(FILE_4.to_string(), test_data(&format!("{FILE_4}.mzML")));

    let result = run_msms(ids, &file_to_mzml).expect("MS2 engine runs over the corpus");
    let mbr = run_mbr(
        &result.peaks_by_file,
        &result.engines_by_file,
        &result.peptide_sequences_to_quantify,
    );

    let targets = mbr.feature_rows.iter().filter(|r| !r.random_rt).count();
    let decoys = mbr.feature_rows.iter().filter(|r| r.random_rt).count();

    let rows = write_feature_table_parquet(&mbr.feature_rows, &path).expect("write parquet");
    println!("wrote {rows} feature rows ({targets} target / {decoys} decoy) to {path}");
}
