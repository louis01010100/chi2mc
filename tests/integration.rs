use std::io::Write;
use std::path::Path;
use std::process::Command;

/// Create a temporary TSV input file with test data and return its path.
fn write_test_input(dir: &Path) -> std::path::PathBuf {
    let path = dir.join("test_input.tsv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "probeset_id\tbatch_name\tn_AA\tn_AB\tn_BB\tn_NC").unwrap();
    // PS001: small n_NC → likely needs MC
    writeln!(f, "PS001\tbatch1\t45\t52\t38\t2").unwrap();
    writeln!(f, "PS001\tbatch2\t48\t50\t40\t1").unwrap();
    writeln!(f, "PS001\tbatch3\t42\t55\t37\t3").unwrap();
    // PS002: large counts → standard test
    writeln!(f, "PS002\tbatch1\t120\t85\t30\t15").unwrap();
    writeln!(f, "PS002\tbatch2\t115\t90\t32\t18").unwrap();
    writeln!(f, "PS002\tbatch3\t118\t88\t31\t17").unwrap();
    // PS003: very small counts → MC
    writeln!(f, "PS003\tbatch1\t2\t3\t1\t1").unwrap();
    writeln!(f, "PS003\tbatch2\t1\t4\t2\t1").unwrap();
    writeln!(f, "PS003\tbatch3\t3\t2\t1\t1").unwrap();
    // PS004: zero NC column
    writeln!(f, "PS004\tbatch1\t50\t60\t40\t0").unwrap();
    writeln!(f, "PS004\tbatch2\t55\t58\t42\t0").unwrap();
    writeln!(f, "PS004\tbatch3\t52\t62\t38\t0").unwrap();
    path
}

fn binary_path() -> std::path::PathBuf {
    // Look for debug binary first, then release
    let debug = Path::new(env!("CARGO_BIN_EXE_chi2mc"));
    debug.to_path_buf()
}

#[test]
fn test_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_test_input(dir.path());
    let output = dir.path().join("output.tsv");

    let status = Command::new(binary_path())
        .arg(input.to_str().unwrap())
        .arg(output.to_str().unwrap())
        .arg("--random-seed")
        .arg("42")
        .arg("--n-simulations")
        .arg("1000")
        .status()
        .expect("failed to run chi2mc");

    assert!(status.success(), "chi2mc exited with non-zero status");
    assert!(output.exists(), "output TSV not created");

    // Read and validate output
    let content = std::fs::read_to_string(&output).unwrap();
    let lines: Vec<&str> = content.lines().collect();

    // Header + 4 probesets
    assert_eq!(lines.len(), 5, "expected 5 lines (header + 4 probesets)");
    assert!(lines[0].starts_with("probeset_id\t"));

    // Check all 4 probesets are present
    let body = &lines[1..];
    let ids: Vec<&str> = body.iter().map(|l| l.split('\t').next().unwrap()).collect();
    assert!(ids.contains(&"PS001"));
    assert!(ids.contains(&"PS002"));
    assert!(ids.contains(&"PS003"));
    assert!(ids.contains(&"PS004"));

    // PS002 should use standard test
    let ps002_line = body.iter().find(|l| l.starts_with("PS002")).unwrap();
    let fields: Vec<&str> = ps002_line.split('\t').collect();
    assert_eq!(fields[4], "standard", "PS002 should use standard test");

    // PS003 should use monte_carlo
    let ps003_line = body.iter().find(|l| l.starts_with("PS003")).unwrap();
    let fields: Vec<&str> = ps003_line.split('\t').collect();
    assert_eq!(fields[4], "monte_carlo", "PS003 should use monte_carlo test");

    // PS004 should have excluded NC genotype
    let ps004_line = body.iter().find(|l| l.starts_with("PS004")).unwrap();
    let fields: Vec<&str> = ps004_line.split('\t').collect();
    assert!(
        fields[9].contains("NC"),
        "PS004 should have NC in excluded_genotypes, got: {}",
        fields[9]
    );

    // Summary report should exist
    let summary = dir.path().join("output_summary.txt");
    assert!(summary.exists(), "summary report not created");
}

#[test]
fn test_reproducibility() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_test_input(dir.path());
    let output1 = dir.path().join("out1.tsv");
    let output2 = dir.path().join("out2.tsv");

    for out in [output1.as_path(), output2.as_path()] {
        let status = Command::new(binary_path())
            .arg(input.to_str().unwrap())
            .arg(out.to_str().unwrap())
            .arg("--random-seed")
            .arg("42")
            .arg("--n-simulations")
            .arg("1000")
            .arg("--n-workers")
            .arg("1")
            .status()
            .expect("failed to run chi2mc");
        assert!(status.success());
    }

    let c1 = std::fs::read_to_string(&output1).unwrap();
    let c2 = std::fs::read_to_string(&output2).unwrap();
    assert_eq!(c1, c2, "outputs should be identical with same seed and single worker");
}

#[test]
fn test_show_progress_flag() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_test_input(dir.path());
    let output = dir.path().join("output.tsv");

    let status = Command::new(binary_path())
        .arg(input.to_str().unwrap())
        .arg(output.to_str().unwrap())
        .arg("--show-progress")
        .arg("--n-simulations")
        .arg("100")
        .status()
        .expect("failed to run chi2mc");

    assert!(status.success());
}
