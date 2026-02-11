use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use crate::stats::{chi2_statistic, compute_expected};

/// Result of a Monte Carlo chi-squared permutation test.
pub struct MonteCarloResult {
    pub statistic: f64,
    pub pvalue: f64,
    pub dof: usize,
    pub n_simulations: usize,
}

/// Run the Monte Carlo chi-squared permutation test.
///
/// `table` is a flat row-major contingency table with dimensions (n_rows, n_cols).
/// Same algorithm as the original Rust PyO3 extension: ChaCha8Rng + Fisher-Yates shuffle.
pub fn monte_carlo_chi2(
    table: &[f64],
    n_rows: usize,
    n_cols: usize,
    n_simulations: usize,
    seed: u64,
) -> MonteCarloResult {
    let expected = compute_expected(table, n_rows, n_cols);
    let obs_chi2 = chi2_statistic(table, &expected);
    let dof = (n_rows - 1) * (n_cols - 1);

    // Build batch and genotype assignment arrays
    let total_samples: usize = table.iter().map(|&v| v as usize).sum();
    let mut batch_assignments = Vec::with_capacity(total_samples);
    let mut genotype_assignments = Vec::with_capacity(total_samples);

    for i in 0..n_rows {
        for j in 0..n_cols {
            let count = table[i * n_cols + j] as usize;
            for _ in 0..count {
                batch_assignments.push(i);
                genotype_assignments.push(j);
            }
        }
    }

    let mut rng = ChaCha8Rng::seed_from_u64(seed);

    // Monte Carlo loop — zero per-iteration allocations
    let mut simulated_table = vec![0.0; n_rows * n_cols];
    let mut count_ge = 0usize;

    for _ in 0..n_simulations {
        // Fisher-Yates shuffle of batch_assignments in-place
        batch_assignments.shuffle(&mut rng);

        // Reconstruct contingency table in-place
        simulated_table.fill(0.0);
        for k in 0..total_samples {
            simulated_table[batch_assignments[k] * n_cols + genotype_assignments[k]] += 1.0;
        }

        let sim_chi2 = chi2_statistic(&simulated_table, &expected);
        if sim_chi2 >= obs_chi2 {
            count_ge += 1;
        }
    }

    let pvalue = count_ge as f64 / n_simulations as f64;

    MonteCarloResult {
        statistic: obs_chi2,
        pvalue,
        dof,
        n_simulations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monte_carlo_reproducibility() {
        let table = vec![
            2.0, 3.0, 1.0, 1.0, 4.0, 2.0, 3.0, 1.0, 3.0, 2.0, 1.0, 1.0,
        ];
        let r1 = monte_carlo_chi2(&table, 3, 4, 1000, 42);
        let r2 = monte_carlo_chi2(&table, 3, 4, 1000, 42);
        assert_eq!(r1.pvalue, r2.pvalue);
        assert_eq!(r1.statistic, r2.statistic);
    }

    #[test]
    fn test_monte_carlo_perfect_table() {
        // Identical rows → chi2 ≈ 0 → p ≈ 1.0
        let table = vec![100.0, 100.0, 100.0, 100.0, 100.0, 100.0];
        let r = monte_carlo_chi2(&table, 2, 3, 1000, 42);
        assert!(r.pvalue > 0.9, "p-value should be near 1.0, got {}", r.pvalue);
    }
}
