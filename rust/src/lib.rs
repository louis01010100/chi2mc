use numpy::ndarray::Array2;
use numpy::PyReadonlyArray2;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

/// Compute expected frequencies from a contingency table.
/// Returns (expected, grand_total).
fn compute_expected(table: &Array2<f64>) -> Array2<f64> {
    let (n_rows, n_cols) = table.dim();
    let grand_total: f64 = table.sum();
    let mut expected = Array2::<f64>::zeros((n_rows, n_cols));

    // row totals
    let row_totals: Vec<f64> = (0..n_rows).map(|i| table.row(i).sum()).collect();
    // col totals
    let col_totals: Vec<f64> = (0..n_cols).map(|j| table.column(j).sum()).collect();

    for i in 0..n_rows {
        for j in 0..n_cols {
            expected[[i, j]] = row_totals[i] * col_totals[j] / grand_total;
        }
    }
    expected
}

/// Compute chi-squared statistic: sum((observed - expected)^2 / expected)
fn chi2_statistic(observed: &Array2<f64>, expected: &Array2<f64>) -> f64 {
    let (n_rows, n_cols) = observed.dim();
    let mut chi2 = 0.0;
    for i in 0..n_rows {
        for j in 0..n_cols {
            let e = expected[[i, j]];
            if e > 0.0 {
                let diff = observed[[i, j]] - e;
                chi2 += diff * diff / e;
            }
        }
    }
    chi2
}

/// Reconstruct a contingency table from batch and genotype assignments.
/// Writes into `table` in-place (must be pre-zeroed).
#[inline]
fn reconstruct_table(
    table: &mut Array2<f64>,
    batch_assignments: &[usize],
    genotype_assignments: &[usize],
) {
    // zero out
    table.fill(0.0);
    for k in 0..batch_assignments.len() {
        table[[batch_assignments[k], genotype_assignments[k]]] += 1.0;
    }
}

/// Run the full Monte Carlo chi-squared permutation test.
///
/// Arguments:
///   contingency_table: int64 numpy array of shape (n_batches, n_genotypes)
///   n_simulations: number of MC iterations
///   seed: optional RNG seed for reproducibility
///
/// Returns a dict with keys: statistic, pvalue, dof, method, n_simulations
#[pyfunction]
#[pyo3(signature = (contingency_table, n_simulations, seed=None))]
fn monte_carlo_chi2<'py>(
    py: Python<'py>,
    contingency_table: PyReadonlyArray2<'py, i64>,
    n_simulations: usize,
    seed: Option<u64>,
) -> PyResult<Bound<'py, PyDict>> {
    let table_i64 = contingency_table.as_array();
    let (n_batches, n_genotypes) = table_i64.dim();

    // Convert to f64
    let observed = table_i64.mapv(|v| v as f64);

    // Compute expected frequencies
    let expected = compute_expected(&observed);

    // Observed chi-squared statistic
    let obs_chi2 = chi2_statistic(&observed, &expected);

    // Degrees of freedom
    let dof = (n_batches - 1) * (n_genotypes - 1);

    // Build batch_assignments and genotype_assignments arrays
    let total_samples: i64 = table_i64.sum();
    let n = total_samples as usize;
    let mut batch_assignments = Vec::with_capacity(n);
    let mut genotype_assignments = Vec::with_capacity(n);

    for i in 0..n_batches {
        for j in 0..n_genotypes {
            let count = table_i64[[i, j]] as usize;
            for _ in 0..count {
                batch_assignments.push(i);
                genotype_assignments.push(j);
            }
        }
    }

    // RNG
    let mut rng = match seed {
        Some(s) => ChaCha8Rng::seed_from_u64(s),
        None => ChaCha8Rng::from_entropy(),
    };

    // Monte Carlo loop — entirely in Rust, zero allocations per iteration
    let mut simulated_table = Array2::<f64>::zeros((n_batches, n_genotypes));
    let mut count_ge = 0usize;

    for _ in 0..n_simulations {
        // Fisher-Yates shuffle of batch_assignments in-place
        batch_assignments.shuffle(&mut rng);

        // Reconstruct contingency table
        reconstruct_table(
            &mut simulated_table,
            &batch_assignments,
            &genotype_assignments,
        );

        // Compute chi-squared for simulated table using same expected frequencies
        let sim_chi2 = chi2_statistic(&simulated_table, &expected);

        if sim_chi2 >= obs_chi2 {
            count_ge += 1;
        }
    }

    let pvalue = count_ge as f64 / n_simulations as f64;

    // Build result dict
    let dict = PyDict::new(py);
    dict.set_item("statistic", obs_chi2)?;
    dict.set_item("pvalue", pvalue)?;
    dict.set_item("dof", dof)?;
    dict.set_item("method", "monte_carlo_rust")?;
    dict.set_item("n_simulations", n_simulations)?;

    Ok(dict)
}

/// Python module definition
#[pymodule]
fn chi2mc_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(monte_carlo_chi2, m)?)?;
    Ok(())
}
