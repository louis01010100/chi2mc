use crate::io::InputRecord;
use crate::monte_carlo::monte_carlo_chi2;
use crate::stats::{check_minimum_expected, chi2_cdf_pvalue, chi2_statistic, compute_expected};

/// Result of processing a single probeset.
#[derive(Debug, Clone)]
pub struct ProbesetResult {
    pub probeset_id: String,
    pub chi2_statistic: Option<f64>,
    pub pvalue: Option<f64>,
    pub dof: Option<f64>,
    pub method: String,
    pub n_batches: usize,
    pub total_samples: i64,
    pub n_simulations: Option<f64>,
    pub error: Option<String>,
    pub excluded_genotypes: Option<String>,
}

/// Process a single probeset: build contingency table, remove zero columns,
/// decide standard vs Monte Carlo test, return result.
pub fn process_probeset(
    probeset_id: &str,
    records: &[InputRecord],
    min_expected_count: usize,
    n_simulations: usize,
    random_seed: u64,
) -> ProbesetResult {
    let n_batches = records.len();
    let genotype_names = ["AA", "AB", "BB", "NC"];

    // Build full contingency table (n_batches x 4) in row-major order
    let mut full_table: Vec<f64> = Vec::with_capacity(n_batches * 4);
    for rec in records {
        full_table.push(rec.n_AA as f64);
        full_table.push(rec.n_AB as f64);
        full_table.push(rec.n_BB as f64);
        full_table.push(rec.n_NC as f64);
    }

    let total_samples: i64 = records
        .iter()
        .map(|r| r.n_AA + r.n_AB + r.n_BB + r.n_NC)
        .sum();

    // Compute column sums to find zero columns
    let mut col_sums = [0.0f64; 4];
    for i in 0..n_batches {
        for j in 0..4 {
            col_sums[j] += full_table[i * 4 + j];
        }
    }

    let non_zero_cols: Vec<usize> = (0..4).filter(|&j| col_sums[j] > 0.0).collect();
    let zero_cols: Vec<usize> = (0..4).filter(|&j| col_sums[j] == 0.0).collect();

    let excluded_genotypes = if zero_cols.is_empty() {
        None
    } else {
        Some(
            zero_cols
                .iter()
                .map(|&j| genotype_names[j])
                .collect::<Vec<_>>()
                .join(", "),
        )
    };

    // If all columns zero, skip
    if non_zero_cols.is_empty() {
        return ProbesetResult {
            probeset_id: probeset_id.to_string(),
            chi2_statistic: None,
            pvalue: None,
            dof: None,
            method: "skipped".to_string(),
            n_batches,
            total_samples,
            n_simulations: None,
            error: Some("All genotype columns are zero".to_string()),
            excluded_genotypes: None,
        };
    }

    let n_cols = non_zero_cols.len();

    // Build reduced table (only non-zero columns)
    let table: Vec<f64> = if n_cols == 4 {
        full_table
    } else {
        let mut reduced = Vec::with_capacity(n_batches * n_cols);
        for i in 0..n_batches {
            for &j in &non_zero_cols {
                reduced.push(full_table[i * 4 + j]);
            }
        }
        reduced
    };

    // Degrees of freedom
    let dof = (n_batches.saturating_sub(1)) * (n_cols.saturating_sub(1));
    if dof == 0 {
        return ProbesetResult {
            probeset_id: probeset_id.to_string(),
            chi2_statistic: None,
            pvalue: None,
            dof: None,
            method: "skipped".to_string(),
            n_batches,
            total_samples,
            n_simulations: None,
            error: Some("Degrees of freedom is zero (single row or column after removal)".to_string()),
            excluded_genotypes,
        };
    }

    // Decide test method
    let use_standard = check_minimum_expected(&table, n_batches, n_cols, min_expected_count as f64);

    if use_standard {
        // Standard chi-squared test via CDF
        let expected = compute_expected(&table, n_batches, n_cols);
        let stat = chi2_statistic(&table, &expected);
        let pval = chi2_cdf_pvalue(stat, dof);

        ProbesetResult {
            probeset_id: probeset_id.to_string(),
            chi2_statistic: Some(stat),
            pvalue: Some(pval),
            dof: Some(dof as f64),
            method: "standard".to_string(),
            n_batches,
            total_samples,
            n_simulations: None,
            error: None,
            excluded_genotypes,
        }
    } else {
        // Monte Carlo test
        let mc = monte_carlo_chi2(&table, n_batches, n_cols, n_simulations, random_seed);

        ProbesetResult {
            probeset_id: probeset_id.to_string(),
            chi2_statistic: Some(mc.statistic),
            pvalue: Some(mc.pvalue),
            dof: Some(mc.dof as f64),
            method: "monte_carlo".to_string(),
            n_batches,
            total_samples,
            n_simulations: Some(mc.n_simulations as f64),
            error: None,
            excluded_genotypes,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_records(data: &[(i64, i64, i64, i64)]) -> Vec<InputRecord> {
        data.iter()
            .enumerate()
            .map(|(i, &(aa, ab, bb, nc))| InputRecord {
                probeset_id: "TEST".to_string(),
                batch_name: format!("batch{}", i + 1),
                n_AA: aa,
                n_AB: ab,
                n_BB: bb,
                n_NC: nc,
            })
            .collect()
    }

    #[test]
    fn test_standard_test() {
        let records = make_records(&[(120, 85, 30, 15), (115, 90, 32, 18), (118, 88, 31, 17)]);
        let r = process_probeset("PS002", &records, 5, 10000, 42);
        assert_eq!(r.method, "standard");
        assert!(r.chi2_statistic.is_some());
        assert!(r.pvalue.is_some());
        assert!(r.n_simulations.is_none());
    }

    #[test]
    fn test_monte_carlo_test() {
        let records = make_records(&[(2, 3, 1, 1), (1, 4, 2, 1), (3, 2, 1, 1)]);
        let r = process_probeset("PS003", &records, 5, 1000, 42);
        assert_eq!(r.method, "monte_carlo");
        assert!(r.chi2_statistic.is_some());
        assert!(r.pvalue.is_some());
        assert!(r.n_simulations.is_some());
    }

    #[test]
    fn test_zero_column() {
        let records = make_records(&[(45, 52, 38, 0), (48, 50, 40, 0), (42, 55, 37, 0)]);
        let r = process_probeset("PS_ZERO", &records, 5, 10000, 42);
        assert!(r.excluded_genotypes.as_deref() == Some("NC"));
        assert_ne!(r.method, "skipped");
    }

    #[test]
    fn test_all_zero() {
        let records = make_records(&[(0, 0, 0, 0), (0, 0, 0, 0)]);
        let r = process_probeset("PS_ALLZERO", &records, 5, 10000, 42);
        assert_eq!(r.method, "skipped");
        assert!(r.error.is_some());
    }
}
