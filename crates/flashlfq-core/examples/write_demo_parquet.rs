//! Writes a small demo peptide-results Parquet file to the path given as the first CLI arg.
//! Used by the P1.16 Python acceptance check (read it back with pyarrow). See PLAN.md P1.16.

use std::collections::HashMap;

use flashlfq_core::detection_type::DetectionType;
use flashlfq_core::parquet_output::write_peptide_results_parquet;
use flashlfq_core::results::{PeptideQuant, PeptideResults};

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: write_demo_parquet <output.parquet>");

    let mut quant: HashMap<String, HashMap<String, PeptideQuant>> = HashMap::new();

    let mut a: HashMap<String, PeptideQuant> = HashMap::new();
    a.insert(
        "fileA".to_string(),
        PeptideQuant {
            intensity: 1.0e6,
            retention_time: 12.5,
            detection_type: DetectionType::MSMS,
        },
    );
    a.insert(
        "fileB".to_string(),
        PeptideQuant {
            intensity: 0.0,
            retention_time: 0.0,
            detection_type: DetectionType::NotDetected,
        },
    );
    quant.insert("PEP[ox]TIDE".to_string(), a);

    let mut b: HashMap<String, PeptideQuant> = HashMap::new();
    b.insert(
        "fileA".to_string(),
        PeptideQuant {
            intensity: 2.5e5,
            retention_time: 30.1,
            detection_type: DetectionType::MSMSAmbiguousPeakfinding,
        },
    );
    quant.insert("PEPTIDEK".to_string(), b);

    let results = PeptideResults { quant };
    let rows = write_peptide_results_parquet(&results, &path).expect("write parquet");
    println!("wrote {rows} rows to {path}");
}
