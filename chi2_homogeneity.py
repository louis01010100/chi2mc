"""
Chi-squared homogeneity test for genotype frequencies across batches.

This module implements a two-tier approach:
1. Standard chi-squared test for contingency tables with all expected frequencies >= 5
2. Monte Carlo chi-squared test for tables with any expected frequency < 5
"""
import polars as pl
import numpy as np
from scipy import stats
from typing import Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing


def load_data(file_path: str) -> pl.DataFrame:
    """
    Load genotype data from TSV file.

    Parameters:
    -----------
    file_path : str
        Path to the TSV file containing genotype data.

    Returns:
    --------
    pl.DataFrame
        DataFrame with columns: probeset_id, batch_name, n_AA, n_AB, n_BB, n_NC
    """
    df = pl.read_csv(file_path, separator='\t')

    # Validate required columns
    required_cols = ['probeset_id', 'batch_name', 'n_AA', 'n_AB', 'n_BB', 'n_NC']
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    return df


def create_contingency_table(probeset_data: pl.DataFrame) -> np.ndarray:
    """
    Create a contingency table for a single probeset.

    Parameters:
    -----------
    probeset_data : pl.DataFrame
        DataFrame containing data for a single probeset across multiple batches.

    Returns:
    --------
    np.ndarray
        Contingency table with shape (n_batches, 4) containing counts for
        AA, AB, BB, NC genotypes across batches.
    """
    # Sort by batch_name for consistency
    probeset_data = probeset_data.sort('batch_name')

    # Extract genotype counts as numpy array
    table = probeset_data.select(['n_AA', 'n_AB', 'n_BB', 'n_NC']).to_numpy()

    return table


def check_minimum_expected_frequency(
    contingency_table: np.ndarray,
    min_count: int = 5
) -> bool:
    """
    Check if all expected frequencies in the contingency table meet minimum threshold.

    Parameters:
    -----------
    contingency_table : np.ndarray
        Contingency table with shape (n_rows, n_cols).
    min_count : int, optional
        Minimum expected count threshold (default: 5).

    Returns:
    --------
    bool
        True if all expected frequencies >= min_count, False otherwise.
    """
    # Calculate expected frequencies
    row_totals = contingency_table.sum(axis=1)
    col_totals = contingency_table.sum(axis=0)
    grand_total = contingency_table.sum()

    # Compute expected frequencies for each cell
    expected = np.outer(row_totals, col_totals) / grand_total

    # Check if all expected frequencies meet minimum threshold
    return np.all(expected >= min_count)


def chi2_standard_test(contingency_table: np.ndarray) -> Dict:
    """
    Perform standard chi-squared test of homogeneity.

    Parameters:
    -----------
    contingency_table : np.ndarray
        Contingency table with shape (n_batches, n_genotypes).

    Returns:
    --------
    dict
        Dictionary containing:
        - statistic: Chi-squared test statistic
        - pvalue: P-value
        - dof: Degrees of freedom
        - method: 'standard'
    """
    # Perform chi-squared test
    chi2_stat, pvalue, dof, expected = stats.chi2_contingency(contingency_table)

    return {
        'statistic': chi2_stat,
        'pvalue': pvalue,
        'dof': dof,
        'method': 'standard'
    }


def chi2_monte_carlo_test(
    contingency_table: np.ndarray,
    n_simulations: int = 10000,
    random_seed: int = None
) -> Dict:
    """
    Perform chi-squared test with Monte Carlo simulation.

    This approach simulates the null distribution by randomly permuting
    the batch assignments while maintaining marginal totals.

    Parameters:
    -----------
    contingency_table : np.ndarray
        Contingency table with shape (n_batches, n_genotypes).
    n_simulations : int, optional
        Number of Monte Carlo simulations (default: 10000).
    random_seed : int, optional
        Random seed for reproducibility (default: None).

    Returns:
    --------
    dict
        Dictionary containing:
        - statistic: Chi-squared test statistic
        - pvalue: Monte Carlo p-value
        - dof: Degrees of freedom
        - method: 'monte_carlo'
        - n_simulations: Number of simulations performed
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    # Calculate observed chi-squared statistic
    chi2_stat, _, dof, expected = stats.chi2_contingency(contingency_table)

    # Prepare data for permutation
    # Create a 1D array where each element represents a genotype for a sample
    n_batches, n_genotypes = contingency_table.shape
    batch_assignments = []
    genotype_assignments = []

    for batch_idx in range(n_batches):
        for genotype_idx in range(n_genotypes):
            count = int(contingency_table[batch_idx, genotype_idx])
            batch_assignments.extend([batch_idx] * count)
            genotype_assignments.extend([genotype_idx] * count)

    batch_assignments = np.array(batch_assignments)
    genotype_assignments = np.array(genotype_assignments)

    # Monte Carlo simulation
    simulated_chi2_stats = []

    for _ in range(n_simulations):
        # Permute batch assignments
        permuted_batches = np.random.permutation(batch_assignments)

        # Reconstruct contingency table
        simulated_table = np.zeros_like(contingency_table)
        for batch_idx, genotype_idx in zip(permuted_batches, genotype_assignments):
            simulated_table[batch_idx, genotype_idx] += 1

        # Calculate chi-squared statistic for simulated table
        try:
            sim_chi2_stat, _, _, _ = stats.chi2_contingency(simulated_table)
            simulated_chi2_stats.append(sim_chi2_stat)
        except ValueError:
            # Skip simulations that result in invalid tables
            continue

    simulated_chi2_stats = np.array(simulated_chi2_stats)

    # Calculate p-value
    pvalue = np.mean(simulated_chi2_stats >= chi2_stat)

    return {
        'statistic': chi2_stat,
        'pvalue': pvalue,
        'dof': dof,
        'method': 'monte_carlo',
        'n_simulations': len(simulated_chi2_stats)
    }


def _process_probeset(
    file_path: str,
    probeset_id: str,
    min_expected_count: int,
    n_simulations: int,
    random_seed: Optional[int]
) -> Dict:
    """
    Worker function to process a single probeset.

    This function is designed to be called by parallel workers.

    Parameters:
    -----------
    file_path : str
        Path to the TSV file containing genotype data.
    probeset_id : str
        Probeset identifier to process.
    min_expected_count : int
        Minimum expected count threshold for using standard test.
    n_simulations : int
        Number of Monte Carlo simulations for small samples.
    random_seed : int, optional
        Random seed for reproducibility.

    Returns:
    --------
    dict
        Result dictionary for this probeset. If the test cannot be computed
        (e.g., due to zero columns), returns a dict with error information.
    """
    # Load data
    df = load_data(file_path)
    probeset_data = df.filter(pl.col('probeset_id') == probeset_id)

    # Create contingency table
    contingency_table = create_contingency_table(probeset_data)

    # Check for zero columns and remove them
    col_sums = contingency_table.sum(axis=0)
    zero_cols = np.where(col_sums == 0)[0]

    genotype_names = ['AA', 'AB', 'BB', 'NC']
    excluded_genotypes = None

    if len(zero_cols) > 0:
        # Find which columns are zero
        excluded_genotypes = [genotype_names[i] for i in zero_cols]

        # Remove zero columns
        non_zero_cols = np.where(col_sums > 0)[0]

        # If all columns are zero, we cannot perform the test
        if len(non_zero_cols) == 0:
            return {
                'probeset_id': probeset_id,
                'chi2_statistic': np.nan,
                'pvalue': np.nan,
                'dof': np.nan,
                'method': 'skipped',
                'n_batches': len(probeset_data),
                'total_samples': int(contingency_table.sum()),
                'n_simulations': np.nan,
                'error': 'All genotype columns are zero',
                'excluded_genotypes': None
            }

        # Keep only non-zero columns
        contingency_table = contingency_table[:, non_zero_cols]

    try:
        # Check if we can use standard test
        use_standard = check_minimum_expected_frequency(
            contingency_table,
            min_count=min_expected_count
        )

        # Perform appropriate test
        if use_standard:
            test_result = chi2_standard_test(contingency_table)
        else:
            test_result = chi2_monte_carlo_test(
                contingency_table,
                n_simulations=n_simulations,
                random_seed=random_seed
            )

        # Store results
        original_total = int(create_contingency_table(probeset_data).sum())
        result_dict = {
            'probeset_id': probeset_id,
            'chi2_statistic': test_result['statistic'],
            'pvalue': test_result['pvalue'],
            'dof': test_result['dof'],
            'method': test_result['method'],
            'n_batches': len(probeset_data),
            'total_samples': original_total,
            'error': None,
            'excluded_genotypes': ', '.join(excluded_genotypes) if excluded_genotypes else None
        }

        if test_result['method'] == 'monte_carlo':
            result_dict['n_simulations'] = float(test_result['n_simulations'])
        else:
            result_dict['n_simulations'] = np.nan

        return result_dict

    except ValueError as e:
        # Handle any other ValueError (e.g., from chi2_contingency)
        original_total = int(create_contingency_table(probeset_data).sum())
        return {
            'probeset_id': probeset_id,
            'chi2_statistic': np.nan,
            'pvalue': np.nan,
            'dof': np.nan,
            'method': 'skipped',
            'n_batches': len(probeset_data),
            'total_samples': original_total,
            'n_simulations': np.nan,
            'error': str(e),
            'excluded_genotypes': None
        }


def run_two_tier_chi2_test(
    file_path: str,
    min_expected_count: int = 5,
    n_simulations: int = 10000,
    random_seed: int = None,
    n_workers: Optional[int] = None
) -> pl.DataFrame:
    """
    Run two-tier chi-squared homogeneity test on genotype data.

    For each probeset, tests whether genotype frequencies are homogeneous
    across batches using:
    - Standard chi-squared test if all expected frequencies >= min_expected_count
    - Monte Carlo chi-squared test otherwise

    Probesets are processed in parallel for improved performance.

    Parameters:
    -----------
    file_path : str
        Path to the TSV file containing genotype data.
    min_expected_count : int, optional
        Minimum expected count threshold for using standard test (default: 5).
    n_simulations : int, optional
        Number of Monte Carlo simulations for small samples (default: 10000).
    random_seed : int, optional
        Random seed for reproducibility (default: None).
    n_workers : int, optional
        Number of parallel workers to use. If None, uses all available CPU cores.
        Set to 1 to disable parallel processing (default: None).

    Returns:
    --------
    pl.DataFrame
        Results DataFrame with columns:
        - probeset_id: Probeset identifier
        - chi2_statistic: Chi-squared test statistic
        - pvalue: P-value
        - dof: Degrees of freedom
        - method: Test method used ('standard' or 'monte_carlo')
        - n_batches: Number of batches
        - total_samples: Total number of samples
    """
    # Load data to get list of probesets
    df = load_data(file_path)
    probeset_ids = df['probeset_id'].unique().to_list()

    # Determine number of workers
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()

    # Process probesets in parallel
    results = []

    if n_workers == 1:
        # Sequential processing (useful for debugging)
        for probeset_id in probeset_ids:
            result = _process_probeset(
                file_path,
                probeset_id,
                min_expected_count,
                n_simulations,
                random_seed
            )
            results.append(result)
    else:
        # Parallel processing using spawn context to prevent deadlocks
        # The 'spawn' method prevents deadlocks that can occur with 'fork' (default on Linux)
        # when running under pytest or in threaded environments
        ctx = multiprocessing.get_context('spawn')
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            # Submit all tasks
            future_to_probeset = {
                executor.submit(
                    _process_probeset,
                    file_path,
                    probeset_id,
                    min_expected_count,
                    n_simulations,
                    random_seed
                ): probeset_id
                for probeset_id in probeset_ids
            }

            # Collect results as they complete
            for future in as_completed(future_to_probeset):
                probeset_id = future_to_probeset[future]
                try:
                    result = future.result()
                    results.append(result)
                except Exception as exc:
                    # Log unexpected exceptions (shouldn't happen with error handling in _process_probeset)
                    print(f'WARNING: Probeset {probeset_id} generated an unexpected exception: {exc}')
                    # Add a skipped result for this probeset
                    results.append({
                        'probeset_id': probeset_id,
                        'chi2_statistic': np.nan,
                        'pvalue': np.nan,
                        'dof': np.nan,
                        'method': 'skipped',
                        'n_batches': np.nan,
                        'total_samples': np.nan,
                        'n_simulations': np.nan,
                        'error': str(exc),
                        'excluded_genotypes': None
                    })

    # Create results DataFrame
    results_df = pl.DataFrame(results)

    # Sort by probeset_id for consistent output
    results_df = results_df.sort('probeset_id')

    return results_df


def main():
    """
    Main function to run the two-tier chi-squared test on the data.
    """
    import time
    import argparse
    from datetime import datetime

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Run two-tier chi-squared homogeneity test on genotype data.'
    )
    parser.add_argument(
        'input_file',
        type=str,
        help='Path to the input TSV file containing genotype data'
    )
    parser.add_argument(
        'output_file',
        type=str,
        help='Path to the output TSV file for results'
    )
    parser.add_argument(
        '--min-expected-count',
        type=int,
        default=5,
        help='Minimum expected count threshold for using standard test (default: 5)'
    )
    parser.add_argument(
        '--n-simulations',
        type=int,
        default=10000,
        help='Number of Monte Carlo simulations for small samples (default: 10000)'
    )
    parser.add_argument(
        '--random-seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    parser.add_argument(
        '--n-workers',
        type=int,
        default=None,
        help='Number of parallel workers to use (default: all CPU cores)'
    )

    args = parser.parse_args()

    print("Two-Tier Chi-Squared Homogeneity Test")
    print("=" * 80)
    print(f"Input file: {args.input_file}")
    print(f"Output file: {args.output_file}")

    # Get number of available CPUs
    n_cpus = multiprocessing.cpu_count()
    print(f"Available CPU cores: {n_cpus}")
    n_workers_display = args.n_workers if args.n_workers else n_cpus
    print(f"Processing probesets in parallel using {n_workers_display} worker(s)...")
    print()

    # Run the test with parallel processing
    start_time = time.time()
    results = run_two_tier_chi2_test(
        args.input_file,
        min_expected_count=args.min_expected_count,
        n_simulations=args.n_simulations,
        random_seed=args.random_seed,
        n_workers=args.n_workers
    )
    elapsed_time = time.time() - start_time

    # Display results
    print("Results:")
    print("=" * 80)
    print(results)
    print("\n")

    # Save results
    results.write_csv(args.output_file, separator='\t')
    print(f"Results saved to: {args.output_file}")

    # Summary statistics
    print("\n" + "=" * 80)
    print("Summary:")
    print(f"Total probesets tested: {len(results)}")
    print(f"Standard chi-squared tests: {results.filter(pl.col('method') == 'standard').height}")
    print(f"Monte Carlo tests: {results.filter(pl.col('method') == 'monte_carlo').height}")

    # Check for probesets with excluded genotypes
    excluded = results.filter(pl.col('excluded_genotypes').is_not_null())
    if len(excluded) > 0:
        print(f"Probesets with excluded genotypes: {len(excluded)}")
        # Only show first 10 in console output
        for i, row in enumerate(excluded.iter_rows(named=True)):
            if i < 10:
                print(f"  - {row['probeset_id']}: excluded {row['excluded_genotypes']} (all zeros)")
            elif i == 10:
                print(f"  ... and {len(excluded) - 10} more (see summary report)")
                break

    # Check for skipped probesets
    skipped = results.filter(pl.col('method') == 'skipped')
    if len(skipped) > 0:
        print(f"Skipped probesets (errors): {len(skipped)}")
        # Show details of skipped probesets
        for row in skipped.iter_rows(named=True):
            print(f"  - {row['probeset_id']}: {row['error']}")

    # Count significant results (excluding skipped)
    valid_results = results.filter(pl.col('method') != 'skipped')
    if len(valid_results) > 0:
        n_significant = valid_results.filter(pl.col('pvalue') < 0.05).height
        print(f"Significant results (p < 0.05): {n_significant}")

    print(f"Processing time: {elapsed_time:.2f} seconds")

    # Generate detailed summary report
    summary_file = args.output_file.replace('.tsv', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CHI-SQUARED HOMOGENEITY TEST - DETAILED SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Input file: {args.input_file}\n")
        f.write(f"Output file: {args.output_file}\n")
        f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Processing time: {elapsed_time:.2f} seconds\n\n")

        f.write("-" * 80 + "\n")
        f.write("OVERALL STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total probesets tested: {len(results)}\n")
        f.write(f"Standard chi-squared tests: {results.filter(pl.col('method') == 'standard').height}\n")
        f.write(f"Monte Carlo tests: {results.filter(pl.col('method') == 'monte_carlo').height}\n")
        f.write(f"Skipped probesets (errors): {len(skipped)}\n\n")

        if len(valid_results) > 0:
            n_significant = valid_results.filter(pl.col('pvalue') < 0.05).height
            f.write(f"Significant results (p < 0.05): {n_significant} ({100*n_significant/len(valid_results):.1f}%)\n\n")

        # Excluded genotypes section
        if len(excluded) > 0:
            f.write("-" * 80 + "\n")
            f.write(f"PROBESETS WITH EXCLUDED GENOTYPES ({len(excluded)} total)\n")
            f.write("-" * 80 + "\n")
            f.write("The following probesets had one or more genotype columns with all zeros.\n")
            f.write("These genotypes were automatically excluded from the chi-squared test.\n\n")

            for row in excluded.iter_rows(named=True):
                f.write(f"  - {row['probeset_id']}: excluded {row['excluded_genotypes']} (all zeros)\n")
            f.write("\n")

        # Skipped probesets section
        if len(skipped) > 0:
            f.write("-" * 80 + "\n")
            f.write(f"SKIPPED PROBESETS ({len(skipped)} total)\n")
            f.write("-" * 80 + "\n")
            f.write("The following probesets could not be tested:\n\n")

            for row in skipped.iter_rows(named=True):
                f.write(f"  - {row['probeset_id']}: {row['error']}\n")
            f.write("\n")

        # Significant results section
        if len(valid_results) > 0:
            significant = valid_results.filter(pl.col('pvalue') < 0.05)
            if len(significant) > 0:
                f.write("-" * 80 + "\n")
                f.write(f"SIGNIFICANT RESULTS (p < 0.05) - {len(significant)} total\n")
                f.write("-" * 80 + "\n")
                f.write("Probesets showing significant heterogeneity across batches:\n\n")

                for row in significant.sort('pvalue').iter_rows(named=True):
                    excluded_note = f" (excluded: {row['excluded_genotypes']})" if row['excluded_genotypes'] else ""
                    f.write(f"  - {row['probeset_id']}: p={row['pvalue']:.6f}, chi2={row['chi2_statistic']:.2f}, method={row['method']}{excluded_note}\n")
                f.write("\n")

        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    print(f"\nDetailed summary report saved to: {summary_file}")


if __name__ == '__main__':
    main()
