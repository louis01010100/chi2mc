use std::collections::BTreeMap;
use std::error::Error;
use std::io::Write;
use std::path::Path;

use serde::Deserialize;

use crate::probeset::ProbesetResult;

/// A single row from the input TSV file.
#[derive(Debug, Clone, Deserialize)]
#[allow(non_snake_case)]
pub struct InputRecord {
    pub probeset_id: String,
    pub batch_name: String,
    pub n_AA: i64,
    pub n_AB: i64,
    pub n_BB: i64,
    pub n_NC: i64,
}

/// Read a TSV input file and group records by probeset_id.
///
/// Returns a BTreeMap so probesets are iterated in sorted order.
/// Each group is sorted by batch_name for reproducibility.
pub fn read_input(path: &Path) -> Result<BTreeMap<String, Vec<InputRecord>>, Box<dyn Error>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)?;

    let mut groups: BTreeMap<String, Vec<InputRecord>> = BTreeMap::new();

    for result in rdr.deserialize() {
        let record: InputRecord = result?;
        groups
            .entry(record.probeset_id.clone())
            .or_default()
            .push(record);
    }

    // Sort each group by batch_name
    for records in groups.values_mut() {
        records.sort_by(|a, b| a.batch_name.cmp(&b.batch_name));
    }

    Ok(groups)
}

/// Format a float for output, using "NaN" for NaN values and ensuring .0 suffix for integers.
fn fmt_float(v: f64) -> String {
    if v.is_nan() {
        "NaN".to_string()
    } else if v.fract() == 0.0 && v.is_finite() {
        format!("{v:.1}")
    } else {
        format!("{v}")
    }
}

/// Format an optional float: NaN when None.
fn fmt_opt(v: Option<f64>) -> String {
    match v {
        Some(val) => fmt_float(val),
        None => "NaN".to_string(),
    }
}

/// Write results to a TSV file, matching the Python output format.
pub fn write_output(path: &Path, results: &[ProbesetResult]) -> Result<(), Box<dyn Error>> {
    let mut f = std::fs::File::create(path)?;

    // Header
    writeln!(
        f,
        "probeset_id\tchi2_statistic\tpvalue\tdof\tmethod\tn_batches\ttotal_samples\tn_simulations\terror\texcluded_genotypes"
    )?;

    for r in results {
        let chi2 = fmt_opt(r.chi2_statistic);
        let pvalue = fmt_opt(r.pvalue);
        let dof = fmt_opt(r.dof);
        let n_batches = fmt_float(r.n_batches as f64);
        let total_samples = fmt_float(r.total_samples as f64);
        let n_simulations = fmt_opt(r.n_simulations);
        let error = r.error.as_deref().unwrap_or("");
        let excluded = r.excluded_genotypes.as_deref().unwrap_or("");

        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.probeset_id,
            chi2,
            pvalue,
            dof,
            r.method,
            n_batches,
            total_samples,
            n_simulations,
            error,
            excluded,
        )?;
    }

    Ok(())
}
