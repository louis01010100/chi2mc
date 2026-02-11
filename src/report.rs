use std::error::Error;
use std::io::Write;
use std::path::Path;

use chrono::Local;

use crate::probeset::ProbesetResult;

/// Generate a detailed summary report text file.
pub fn write_summary_report(
    path: &Path,
    input_file: &str,
    output_file: &str,
    elapsed_secs: f64,
    results: &[ProbesetResult],
) -> Result<(), Box<dyn Error>> {
    let mut f = std::fs::File::create(path)?;

    writeln!(f, "{}", "=".repeat(80))?;
    writeln!(f, "CHI-SQUARED HOMOGENEITY TEST - DETAILED SUMMARY REPORT")?;
    writeln!(f, "{}\n", "=".repeat(80))?;

    writeln!(f, "Input file: {input_file}")?;
    writeln!(f, "Output file: {output_file}")?;
    writeln!(
        f,
        "Analysis date: {}",
        Local::now().format("%Y-%m-%d %H:%M:%S")
    )?;
    writeln!(f, "Processing time: {elapsed_secs:.2} seconds\n")?;

    let n_total = results.len();
    let n_standard = results.iter().filter(|r| r.method == "standard").count();
    let n_mc = results.iter().filter(|r| r.method == "monte_carlo").count();
    let skipped: Vec<_> = results.iter().filter(|r| r.method == "skipped").collect();
    let valid: Vec<_> = results.iter().filter(|r| r.method != "skipped").collect();

    writeln!(f, "{}", "-".repeat(80))?;
    writeln!(f, "OVERALL STATISTICS")?;
    writeln!(f, "{}", "-".repeat(80))?;
    writeln!(f, "Total probesets tested: {n_total}")?;
    writeln!(f, "Standard chi-squared tests: {n_standard}")?;
    writeln!(f, "Monte Carlo tests: {n_mc}")?;
    writeln!(f, "Skipped probesets (errors): {}\n", skipped.len())?;

    if !valid.is_empty() {
        let n_significant = valid
            .iter()
            .filter(|r| r.pvalue.map_or(false, |p| p < 0.05))
            .count();
        let pct = 100.0 * n_significant as f64 / valid.len() as f64;
        writeln!(f, "Significant results (p < 0.05): {n_significant} ({pct:.1}%)\n")?;
    }

    // Excluded genotypes section
    let excluded: Vec<_> = results
        .iter()
        .filter(|r| r.excluded_genotypes.is_some())
        .collect();
    if !excluded.is_empty() {
        writeln!(f, "{}", "-".repeat(80))?;
        writeln!(
            f,
            "PROBESETS WITH EXCLUDED GENOTYPES ({} total)",
            excluded.len()
        )?;
        writeln!(f, "{}", "-".repeat(80))?;
        writeln!(
            f,
            "The following probesets had one or more genotype columns with all zeros."
        )?;
        writeln!(
            f,
            "These genotypes were automatically excluded from the chi-squared test.\n"
        )?;
        for r in &excluded {
            writeln!(
                f,
                "  - {}: excluded {} (all zeros)",
                r.probeset_id,
                r.excluded_genotypes.as_deref().unwrap()
            )?;
        }
        writeln!(f)?;
    }

    // Skipped probesets section
    if !skipped.is_empty() {
        writeln!(f, "{}", "-".repeat(80))?;
        writeln!(f, "SKIPPED PROBESETS ({} total)", skipped.len())?;
        writeln!(f, "{}", "-".repeat(80))?;
        writeln!(f, "The following probesets could not be tested:\n")?;
        for r in &skipped {
            writeln!(
                f,
                "  - {}: {}",
                r.probeset_id,
                r.error.as_deref().unwrap_or("Unknown error")
            )?;
        }
        writeln!(f)?;
    }

    // Significant results section
    if !valid.is_empty() {
        let mut significant: Vec<_> = valid
            .iter()
            .filter(|r| r.pvalue.map_or(false, |p| p < 0.05))
            .collect();
        if !significant.is_empty() {
            significant.sort_by(|a, b| {
                a.pvalue
                    .unwrap_or(f64::INFINITY)
                    .partial_cmp(&b.pvalue.unwrap_or(f64::INFINITY))
                    .unwrap()
            });
            writeln!(f, "{}", "-".repeat(80))?;
            writeln!(
                f,
                "SIGNIFICANT RESULTS (p < 0.05) - {} total",
                significant.len()
            )?;
            writeln!(f, "{}", "-".repeat(80))?;
            writeln!(
                f,
                "Probesets showing significant heterogeneity across batches:\n"
            )?;
            for r in &significant {
                let excl_note = match &r.excluded_genotypes {
                    Some(e) => format!(" (excluded: {e})"),
                    None => String::new(),
                };
                writeln!(
                    f,
                    "  - {}: p={:.6}, chi2={:.2}, method={}{}",
                    r.probeset_id,
                    r.pvalue.unwrap_or(f64::NAN),
                    r.chi2_statistic.unwrap_or(f64::NAN),
                    r.method,
                    excl_note,
                )?;
            }
            writeln!(f)?;
        }
    }

    writeln!(f, "{}", "=".repeat(80))?;
    writeln!(f, "END OF REPORT")?;
    writeln!(f, "{}", "=".repeat(80))?;

    Ok(())
}
