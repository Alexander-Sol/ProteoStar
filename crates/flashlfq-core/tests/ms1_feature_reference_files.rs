//! Parses the three real-world `_ms1.feature` dialects against files produced by the actual
//! tools, not hand-written approximations of them.
//!
//! Fixtures in `data/ms1_feature/` are verbatim copies of mzLib's reader-test corpus
//! (`Test/FileReadingTests/ExternalFileTypes/`), so anything mzLib is known to read is
//! covered here. The unit tests in `src/ms1_feature.rs` use inlined excerpts for speed; this
//! is the guard against those excerpts drifting from the genuine article — a real file's
//! trailing tabs, scientific notation, and column ordering all have to survive.

use flashlfq_core::ms1_feature::{read_ms1_feature_file, Ms1FeatureDialect, RtUnit};

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("data/ms1_feature")
        .join(name)
}

#[test]
fn reads_real_topfd_v162_file() {
    let set = read_ms1_feature_file(fixture("topfd_v1.6.2.ms1.feature")).unwrap();
    assert_eq!(set.dialect, Ms1FeatureDialect::TopFdV1_6);
    assert_eq!(set.rt_unit, RtUnit::Seconds);
    assert_eq!(set.records.len(), 4);

    let first = &set.records[0];
    assert!((first.mass - 10835.85272090354).abs() < 1e-9);
    assert!((first.intensity - 10849947123.04).abs() < 1e-2);
    // File holds 2390.8 s; records are always minutes.
    assert!((first.rt_apex - 2390.8 / 60.0).abs() < 1e-9);
    assert!((first.apex_intensity.unwrap() - 912795138.8).abs() < 1e-1);
    assert_eq!((first.charge_min, first.charge_max), (7, 17));

    // Every row: a sane mass, an ordered RT window, and an apex inside it.
    for r in &set.records {
        assert!(r.mass > 0.0);
        assert!(r.rt_begin <= r.rt_end, "id {}", r.id);
        assert!(r.rt_apex >= r.rt_begin && r.rt_apex <= r.rt_end, "id {}", r.id);
        assert!(r.charge_min <= r.charge_max);
        assert!(r.apex_intensity.is_some(), "v1.6 always carries Apex_intensity");
    }
}

#[test]
fn reads_real_topfd_v170_file() {
    let set = read_ms1_feature_file(fixture("topfd_v1.7.0.ms1.feature")).unwrap();
    assert_eq!(set.dialect, Ms1FeatureDialect::TopFdV1_7);
    assert_eq!(set.rt_unit, RtUnit::Minutes, "v1.7 writes minutes, unlike v1.6");
    assert_eq!(set.records.len(), 4);

    let first = &set.records[0];
    assert!((first.mass - 6983.005647734753).abs() < 1e-9);
    assert!((first.rt_begin - 44.45833333333334).abs() < 1e-12);
    assert_eq!(
        first.file_name.as_deref(),
        Some("E:/Projects/LVS_TD_Yeast/05-26-17_B7A_yeast_td_fract7_rep1.mzML")
    );
    assert_eq!((first.scan_min, first.scan_max), (Some(1539), Some(1593)));
    assert_eq!(first.apex_scan, Some(1557));
    assert_eq!(first.rep_charge, Some(10));
    assert_eq!(first.envelope_num, Some(129));
    assert!((first.ec_score.unwrap() - 0.9979024819216775).abs() < 1e-12);

    // The wide-schema fields are what distinguishes this dialect: all rows must have them,
    // and the trailing tab on every line must not shift a column.
    for r in &set.records {
        assert!(r.scan_min.is_some() && r.apex_scan.is_some() && r.ec_score.is_some());
        assert!(r.charge_min <= r.charge_max && r.charge_min > 0);
        assert!(r.rep_charge.unwrap() >= r.charge_min && r.rep_charge.unwrap() <= r.charge_max);
        // Scan and time bounds must agree on ordering — a shifted column would break this.
        assert!(r.scan_min.unwrap() <= r.scan_max.unwrap());
        assert!(r.rt_begin <= r.rt_end);
    }
}

#[test]
fn reads_real_flashdeconv_file() {
    let set = read_ms1_feature_file(fixture("flashdeconv_openms3.0.0.ms1.feature")).unwrap();
    assert_eq!(set.dialect, Ms1FeatureDialect::FlashDeconv);
    assert_eq!(set.rt_unit, RtUnit::Seconds);
    assert_eq!(set.records.len(), 7);

    let first = &set.records[0];
    // Intensities here are in scientific notation (1.21E+10).
    assert!((first.intensity - 1.21e10).abs() < 1.0);
    assert!((first.rt_apex - 2390.8 / 60.0).abs() < 1e-9);
    assert_eq!((first.charge_min, first.charge_max), (7, 18));

    // The defining absence: FLASHDeconv has no apex-intensity column, and we must not
    // invent a zero for it.
    assert!(set.records.iter().all(|r| r.apex_intensity.is_none()));
}

#[test]
fn every_real_file_round_trips_in_its_own_dialect() {
    for name in [
        "topfd_v1.6.2.ms1.feature",
        "topfd_v1.7.0.ms1.feature",
        "flashdeconv_openms3.0.0.ms1.feature",
    ] {
        let set = read_ms1_feature_file(fixture(name)).unwrap();
        let mut buf: Vec<u8> = Vec::new();
        set.write(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let back = flashlfq_core::ms1_feature::read_ms1_feature(&text).unwrap();

        assert_eq!(back.dialect, set.dialect, "{name}: dialect changed on write");
        assert_eq!(back.rt_unit, set.rt_unit, "{name}: rt unit changed on write");
        assert_eq!(back.records.len(), set.records.len(), "{name}: row count");
        for (a, b) in set.records.iter().zip(&back.records) {
            assert!((a.mass - b.mass).abs() < 1e-9, "{name}: mass drifted");
            assert!((a.intensity - b.intensity).abs() < 1e-3, "{name}: intensity drifted");
            assert!((a.rt_begin - b.rt_begin).abs() < 1e-9, "{name}: rt_begin drifted");
            assert!((a.rt_end - b.rt_end).abs() < 1e-9, "{name}: rt_end drifted");
            assert!((a.rt_apex - b.rt_apex).abs() < 1e-9, "{name}: rt_apex drifted");
            assert_eq!(a.charge_min, b.charge_min, "{name}: charge_min");
            assert_eq!(a.charge_max, b.charge_max, "{name}: charge_max");
            assert_eq!(a.apex_intensity, b.apex_intensity, "{name}: apex intensity");
            assert_eq!(a.rep_charge, b.rep_charge, "{name}: rep_charge");
            assert_eq!(a.scan_min, b.scan_min, "{name}: scan_min");
        }
    }
}
