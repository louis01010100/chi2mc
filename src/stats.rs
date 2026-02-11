use statrs::distribution::{ChiSquared, ContinuousCDF};

/// Compute expected frequencies for a contingency table stored as a flat Vec
/// in row-major order with dimensions (n_rows, n_cols).
///
/// Returns a flat Vec of the same dimensions.
pub fn compute_expected(table: &[f64], n_rows: usize, n_cols: usize) -> Vec<f64> {
    let grand_total: f64 = table.iter().sum();

    let row_totals: Vec<f64> = (0..n_rows)
        .map(|i| {
            let start = i * n_cols;
            table[start..start + n_cols].iter().sum()
        })
        .collect();

    let col_totals: Vec<f64> = (0..n_cols)
        .map(|j| (0..n_rows).map(|i| table[i * n_cols + j]).sum())
        .collect();

    let mut expected = vec![0.0; n_rows * n_cols];
    for i in 0..n_rows {
        for j in 0..n_cols {
            expected[i * n_cols + j] = row_totals[i] * col_totals[j] / grand_total;
        }
    }
    expected
}

/// Compute chi-squared statistic: sum((observed - expected)^2 / expected).
pub fn chi2_statistic(observed: &[f64], expected: &[f64]) -> f64 {
    let mut chi2 = 0.0;
    for (o, e) in observed.iter().zip(expected.iter()) {
        if *e > 0.0 {
            let diff = o - e;
            chi2 += diff * diff / e;
        }
    }
    chi2
}

/// Check if all expected frequencies meet the minimum threshold.
pub fn check_minimum_expected(table: &[f64], n_rows: usize, n_cols: usize, min_count: f64) -> bool {
    let expected = compute_expected(table, n_rows, n_cols);
    expected.iter().all(|&e| e >= min_count)
}

/// Compute p-value from chi-squared statistic and degrees of freedom
/// using the chi-squared CDF: p = 1 - CDF(chi2_stat).
pub fn chi2_cdf_pvalue(chi2_stat: f64, dof: usize) -> f64 {
    if dof == 0 {
        return f64::NAN;
    }
    let dist = ChiSquared::new(dof as f64).unwrap();
    1.0 - dist.cdf(chi2_stat)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_expected() {
        // 2x2 table: [[10, 20], [30, 40]]
        let table = vec![10.0, 20.0, 30.0, 40.0];
        let expected = compute_expected(&table, 2, 2);
        let grand = 100.0;
        assert!((expected[0] - 30.0 * 40.0 / grand).abs() < 1e-10);
        assert!((expected[1] - 30.0 * 60.0 / grand).abs() < 1e-10);
        assert!((expected[2] - 70.0 * 40.0 / grand).abs() < 1e-10);
        assert!((expected[3] - 70.0 * 60.0 / grand).abs() < 1e-10);
    }

    #[test]
    fn test_chi2_statistic_perfect_match() {
        let obs = vec![10.0, 20.0, 30.0, 40.0];
        let exp = vec![10.0, 20.0, 30.0, 40.0];
        assert!((chi2_statistic(&obs, &exp)).abs() < 1e-10);
    }

    #[test]
    fn test_chi2_cdf_pvalue() {
        // chi2 = 0, dof > 0 -> p = 1.0
        let p = chi2_cdf_pvalue(0.0, 1);
        assert!((p - 1.0).abs() < 1e-10);

        // Large chi2 -> p near 0
        let p = chi2_cdf_pvalue(100.0, 1);
        assert!(p < 1e-10);
    }

    #[test]
    fn test_check_minimum_expected() {
        // Large table: all expected >= 5
        let table = vec![120.0, 85.0, 30.0, 15.0, 115.0, 90.0, 32.0, 18.0];
        assert!(check_minimum_expected(&table, 2, 4, 5.0));

        // Small table: some expected < 5
        let table = vec![2.0, 3.0, 1.0, 1.0, 1.0, 4.0, 2.0, 1.0];
        assert!(!check_minimum_expected(&table, 2, 4, 5.0));
    }
}
