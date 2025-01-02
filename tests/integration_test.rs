use std::{fs, path::PathBuf, process::Command};

#[test]
fn test_methylation_pattern() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let data_dir = PathBuf::from(manifest_dir).join("tests/data");

    let pileup = data_dir.join("geobacillus-plasmids.pileup.bed");
    let assembly = data_dir.join("geobacillus-plasmids.assembly.fasta");
    let expected_out = data_dir.join("expected_out.tsv");

    let out_file = PathBuf::from(manifest_dir)
        .join("target")
        .join("test_out.tsv");

    let status = Command::new("cargo")
        .args(&[
            "run",
            "--quiet",
            "--",
            "methylation-pattern",
            "-p",
            pileup.to_str().unwrap(),
            "-a",
            assembly.to_str().unwrap(),
            "-m",
            "GATC_a_1",
            "GATC_m_3",
            "RGATCY_a_2",
            "-o",
            out_file.to_str().unwrap(),
            "--batch-size",
            "2",
            "--min-valid-read-coverage",
            "3",
        ])
        .status()
        .expect("Failed to execute cargo run");

    assert!(
        status.success(),
        "Process ended with non-success status: {:?}",
        status
    );

    let actual = fs::read_to_string(&out_file).expect("Could not read output file");
    let expected = fs::read_to_string(&expected_out).expect("Could not read expected output file");

    assert_eq!(
        actual.trim(),
        expected.trim(),
        "Output did not match expected"
    );
}
