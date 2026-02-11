"""
Chi-squared homogeneity test for genotype frequencies across batches.

This module implements a two-tier approach:
1. Standard chi-squared test for contingency tables with all expected frequencies >= 5
2. Monte Carlo chi-squared test for tables with any expected frequency < 5
"""
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Optional

import numpy as np
import polars as pl
from scipy import stats

# Try to import Rust extension for maximum MC performance
try:
    import chi2mc_rust
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False

# Try to import numba for JIT compilation (optional but highly recommended for speed)
try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Define a no-op decorator if numba is not available
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return decorator

# Try to import tqdm for progress bars (optional)
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False


@njit
def _reconstruct_table_numba(permuted_batches, genotype_assignments, n_batches, n_genotypes):
    """
    JIT-compiled function to reconstruct contingency table from permuted batches.

    This is significantly faster than the Python loop version.
    """
    table = np.zeros((n_batches, n_genotypes), dtype=np.int64)
    for i in range(len(permuted_batches)):
        table[permuted_batches[i], genotype_assignments[i]] += 1
    return table


@njit
def _compute_chi2_statistic_numba(observed, expected):
    """
    JIT-compiled chi-squared statistic calculation.

    Computes: sum((observed - expected)^2 / expected) for non-zero expected values.
    """
    chi2 = 0.0
    for i in range(observed.shape[0]):
        for j in range(observed.shape[1]):
            if expected[i, j] > 0:
                chi2 += (observed[i, j] - expected[i, j]) ** 2 / expected[i, j]
    return chi2


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
    random_seed: int = None,
    use_numba: bool = None
) -> Dict:
    """
    Perform chi-squared test with Monte Carlo simulation.

    This approach simulates the null distribution by randomly permuting
    the batch assignments while maintaining marginal totals.

    OPTIMIZED: Now uses vectorized operations and optional numba JIT compilation.

    Parameters:
    -----------
    contingency_table : np.ndarray
        Contingency table with shape (n_batches, n_genotypes).
    n_simulations : int, optional
        Number of Monte Carlo simulations (default: 10000).
    random_seed : int, optional
        Random seed for reproducibility (default: None).
    use_numba : bool, optional
        Whether to use numba JIT compilation if available (default: None = auto-detect).
        Set to False to disable numba even if installed.

    Returns:
    --------
    dict
        Dictionary containing:
        - statistic: Chi-squared test statistic
        - pvalue: Monte Carlo p-value
        - dof: Degrees of freedom
        - method: 'monte_carlo' (or 'monte_carlo_numba' if numba used)
        - n_simulations: Number of simulations performed
    """
    # Prefer Rust extension when available — entire MC loop runs in native code
    if RUST_AVAILABLE:
        seed_val = int(random_seed) if random_seed is not None else None
        table_i64 = contingency_table.astype(np.int64, copy=False)
        result = chi2mc_rust.monte_carlo_chi2(table_i64, n_simulations, seed=seed_val)
        return dict(result)

    if random_seed is not None:
        np.random.seed(random_seed)

    # Auto-detect numba usage if not specified
    if use_numba is None:
        use_numba = NUMBA_AVAILABLE

    # Calculate observed chi-squared statistic and expected frequencies
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

    batch_assignments = np.array(batch_assignments, dtype=np.int64)
    genotype_assignments = np.array(genotype_assignments, dtype=np.int64)

    # Monte Carlo simulation - OPTIMIZED with vectorization and optional numba
    simulated_chi2_stats = []

    if use_numba and NUMBA_AVAILABLE:
        # Use numba-optimized version for maximum speed
        for _ in range(n_simulations):
            permuted_batches = np.random.permutation(batch_assignments)
            simulated_table = _reconstruct_table_numba(
                permuted_batches, genotype_assignments, n_batches, n_genotypes
            )

            # Calculate chi-squared statistic using numba
            sim_chi2_stat = _compute_chi2_statistic_numba(
                simulated_table.astype(np.float64), expected
            )
            simulated_chi2_stats.append(sim_chi2_stat)

        method_name = 'monte_carlo_numba'
    else:
        # Use vectorized numpy version (still much faster than original)
        for _ in range(n_simulations):
            permuted_batches = np.random.permutation(batch_assignments)

            # Vectorized contingency table reconstruction
            combined_indices = permuted_batches * n_genotypes + genotype_assignments
            flat_counts = np.bincount(combined_indices, minlength=n_batches * n_genotypes)
            simulated_table = flat_counts.reshape(n_batches, n_genotypes)

            # Calculate chi-squared statistic for simulated table
            try:
                sim_chi2_stat, _, _, _ = stats.chi2_contingency(simulated_table)
                simulated_chi2_stats.append(sim_chi2_stat)
            except ValueError:
                # Skip simulations that result in invalid tables
                continue

        method_name = 'monte_carlo'

    simulated_chi2_stats = np.array(simulated_chi2_stats)

    # Calculate p-value
    pvalue = np.mean(simulated_chi2_stats >= chi2_stat)

    return {
        'statistic': chi2_stat,
        'pvalue': pvalue,
        'dof': dof,
        'method': method_name,
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
    DEPRECATED: Use _process_probeset_with_data for better performance.

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
    return _process_probeset_core(probeset_data, probeset_id, min_expected_count, n_simulations, random_seed)


def _process_probeset_with_data(
    df: pl.DataFrame,
    probeset_id: str,
    min_expected_count: int,
    n_simulations: int,
    random_seed: Optional[int]
) -> Dict:
    """
    OPTIMIZED worker function to process a single probeset with pre-loaded data.

    This version receives the full DataFrame instead of loading it from disk,
    which is much faster when processing multiple probesets.

    Parameters:
    -----------
    df : pl.DataFrame
        Pre-loaded DataFrame containing genotype data for all probesets.
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
        Result dictionary for this probeset.
    """
    probeset_data = df.filter(pl.col('probeset_id') == probeset_id)
    return _process_probeset_core(probeset_data, probeset_id, min_expected_count, n_simulations, random_seed)


def _process_probeset_batch(
    df: pl.DataFrame,
    probeset_ids: list,
    min_expected_count: int,
    n_simulations: int,
    random_seed: Optional[int]
) -> list:
    """
    OPTIMIZED worker function to process a batch of probesets.

    This version processes multiple probesets in a single subprocess,
    which dramatically reduces multiprocessing overhead for large datasets.

    Parameters:
    -----------
    df : pl.DataFrame
        Pre-loaded DataFrame containing genotype data for all probesets.
    probeset_ids : list
        List of probeset identifiers to process in this batch.
    min_expected_count : int
        Minimum expected count threshold for using standard test.
    n_simulations : int
        Number of Monte Carlo simulations for small samples.
    random_seed : int, optional
        Random seed for reproducibility.

    Returns:
    --------
    list
        List of result dictionaries, one for each probeset in the batch.
    """
    results = []
    for probeset_id in probeset_ids:
        result = _process_probeset_with_data(
            df,
            probeset_id,
            min_expected_count,
            n_simulations,
            random_seed
        )
        results.append(result)
    return results


def _process_probeset_core(
    probeset_data: pl.DataFrame,
    probeset_id: str,
    min_expected_count: int,
    n_simulations: int,
    random_seed: Optional[int]
) -> Dict:
    """
    Core processing logic for a single probeset.

    Parameters:
    -----------
    probeset_data : pl.DataFrame
        DataFrame containing data for this specific probeset.
    probeset_id : str
        Probeset identifier.
    min_expected_count : int
        Minimum expected count threshold for using standard test.
    n_simulations : int
        Number of Monte Carlo simulations for small samples.
    random_seed : int, optional
        Random seed for reproducibility.

    Returns:
    --------
    dict
        Result dictionary for this probeset.
    """

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

        if 'monte_carlo' in test_result['method']:  # Handles both 'monte_carlo' and 'monte_carlo_numba'
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


def _format_progress_message(result: Dict, n_completed: int, n_total: int) -> str:
    """
    Format a progress message for a completed probeset.

    Parameters:
    -----------
    result : dict
        Result dictionary for the probeset.
    n_completed : int
        Number of probesets completed so far.
    n_total : int
        Total number of probesets to process.

    Returns:
    --------
    str
        Formatted progress message.
    """
    probeset_id = result['probeset_id']
    method = result['method']

    if method == 'skipped':
        return f"[{n_completed}/{n_total}] {probeset_id}: SKIPPED - {result.get('error', 'Unknown error')}"

    chi2 = result['chi2_statistic']
    pval = result['pvalue']

    # Format p-value
    if pval < 0.0001:
        pval_str = f"{pval:.2e}"
    else:
        pval_str = f"{pval:.6f}"

    # Format chi2
    chi2_str = f"{chi2:.4f}"

    # Indicate significance
    sig_marker = " ***" if pval < 0.001 else " **" if pval < 0.01 else " *" if pval < 0.05 else ""

    return f"[{n_completed}/{n_total}] {probeset_id}: \u03c7\u00b2={chi2_str}, p={pval_str}{sig_marker} ({method})"


def _print_progress(result: Dict, n_completed: int, n_total: int, verbose: bool = False):
    """
    Print progress for a completed probeset.

    Parameters:
    -----------
    result : dict
        Result dictionary for the probeset.
    n_completed : int
        Number of probesets completed so far.
    n_total : int
        Total number of probesets to process.
    verbose : bool
        Whether to show detailed information.
    """
    if verbose:
        print(_format_progress_message(result, n_completed, n_total))
    else:
        # Simple progress counter
        percent = 100 * n_completed / n_total
        print(f"Progress: {n_completed}/{n_total} ({percent:.1f}%)", end='\r')
        if n_completed == n_total:
            print()  # New line when complete


def run_two_tier_chi2_test(
    file_path: str,
    min_expected_count: int = 5,
    n_simulations: int = 10000,
    random_seed: int = None,
    n_workers: Optional[int] = None,
    use_optimized: bool = True,
    show_progress: bool = False,
    verbose_progress: bool = False,
    batch_size: int = 1000
) -> pl.DataFrame:
    """
    Run two-tier chi-squared homogeneity test on genotype data.

    For each probeset, tests whether genotype frequencies are homogeneous
    across batches using:
    - Standard chi-squared test if all expected frequencies >= min_expected_count
    - Monte Carlo chi-squared test otherwise

    Probesets are processed in parallel for improved performance.

    OPTIMIZED: Now uses vectorized operations, optional numba JIT compilation,
    and optimized data loading.

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
    use_optimized : bool, optional
        Whether to use optimized data loading (pass data to workers instead of file path).
        Default: True (recommended for better performance).
    show_progress : bool, optional
        Whether to show progress bar/messages during processing (default: False).
        If tqdm is installed, uses a progress bar. Otherwise, prints progress messages.
    verbose_progress : bool, optional
        Whether to show detailed progress with chi-squared and p-values (default: False).
        Only used when show_progress=True. Shows each result as it completes.
    batch_size : int, optional
        Number of probesets to process in each subprocess (default: 1000).
        Larger batches reduce multiprocessing overhead but may reduce parallelism.
        For datasets with < 10,000 probesets, 1000 is optimal.
        For datasets with > 100,000 probesets, consider 5000-10000.

    Returns:
    --------
    pl.DataFrame
        Results DataFrame with columns:
        - probeset_id: Probeset identifier
        - chi2_statistic: Chi-squared test statistic
        - pvalue: P-value
        - dof: Degrees of freedom
        - method: Test method used ('standard', 'monte_carlo', or 'monte_carlo_numba')
        - n_batches: Number of batches
        - total_samples: Total number of samples
    """
    # Load data to get list of probesets
    df = load_data(file_path)
    probeset_ids = df['probeset_id'].unique().to_list()
    n_total = len(probeset_ids)

    # Determine number of workers
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()

    # Create batches of probesets for more efficient processing
    # Each worker will process a batch of probesets to reduce overhead
    batches = []
    for i in range(0, n_total, batch_size):
        batch = probeset_ids[i:i + batch_size]
        batches.append(batch)

    # Process probesets in parallel
    results = []

    # Setup progress tracking
    # IMPORTANT: Disable progress bar in multiprocessing to prevent I/O contention
    # Use position and leave parameters for better multiprocessing compatibility
    if show_progress and TQDM_AVAILABLE and n_workers > 1:
        # For parallel processing, use special tqdm settings to avoid output contention
        progress_bar = tqdm(
            total=n_total,
            desc="Processing probesets",
            unit="probeset",
            position=0,
            leave=True,
            disable=False  # Keep enabled but with safe settings
        )
    elif show_progress and TQDM_AVAILABLE:
        # For single worker, use standard progress bar
        progress_bar = tqdm(total=n_total, desc="Processing probesets", unit="probeset")
    else:
        progress_bar = None

    if n_workers == 1:
        # Sequential processing (useful for debugging)
        # Process in batches even with single worker for consistency
        n_completed = 0
        for batch_idx, batch in enumerate(batches):
            if use_optimized:
                # Process batch
                batch_results = _process_probeset_batch(
                    df,
                    batch,
                    min_expected_count,
                    n_simulations,
                    random_seed
                )
            else:
                # Process individually (legacy mode)
                batch_results = []
                for probeset_id in batch:
                    result = _process_probeset(
                        file_path,
                        probeset_id,
                        min_expected_count,
                        n_simulations,
                        random_seed
                    )
                    batch_results.append(result)

            # Add batch results
            results.extend(batch_results)

            # Show progress for each probeset in batch
            if show_progress:
                for result in batch_results:
                    n_completed += 1
                    if progress_bar:
                        progress_bar.update(1)
                        if verbose_progress:
                            progress_bar.write(_format_progress_message(result, n_completed, n_total))
                    else:
                        _print_progress(result, n_completed, n_total, verbose_progress)
    else:
        # Parallel processing using spawn context to prevent deadlocks
        # The 'spawn' method prevents deadlocks that can occur with 'fork' (default on Linux)
        # when running under pytest or in threaded environments
        ctx = multiprocessing.get_context('spawn')
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
            # Submit all batch tasks
            if use_optimized:
                # OPTIMIZED: Process batches - each worker handles multiple probesets
                future_to_batch = {
                    executor.submit(
                        _process_probeset_batch,
                        df,
                        batch,
                        min_expected_count,
                        n_simulations,
                        random_seed
                    ): batch
                    for batch in batches
                }
            else:
                # Legacy: Process individually (not recommended for large datasets)
                future_to_batch = {}
                for batch in batches:
                    batch_futures = {
                        executor.submit(
                            _process_probeset,
                            file_path,
                            probeset_id,
                            min_expected_count,
                            n_simulations,
                            random_seed
                        ): [probeset_id]  # Wrap in list for consistency
                        for probeset_id in batch
                    }
                    future_to_batch.update(batch_futures)

            # Collect results as they complete
            n_completed = 0
            for future in as_completed(future_to_batch):
                batch = future_to_batch[future]
                try:
                    batch_results = future.result()

                    # Handle both batch and single results
                    if not isinstance(batch_results, list):
                        batch_results = [batch_results]

                    # Add batch results
                    results.extend(batch_results)

                    # Show progress for each probeset in batch
                    for result in batch_results:
                        n_completed += 1

                        if show_progress:
                            if progress_bar:
                                progress_bar.update(1)
                                # Only show verbose progress every N results in parallel mode to reduce I/O contention
                                if verbose_progress:
                                    # In parallel mode, only show every 10th result to reduce contention
                                    if n_workers == 1 or n_completed % max(1, min(10, n_workers)) == 0:
                                        progress_bar.write(_format_progress_message(result, n_completed, n_total))
                            else:
                                _print_progress(result, n_completed, n_total, verbose_progress)

                except Exception as exc:
                    # Log unexpected exceptions (shouldn't happen with error handling in _process_probeset)
                    print(f'WARNING: Batch processing generated an unexpected exception: {exc}')
                    # Add skipped results for entire batch
                    for probeset_id in batch:
                        result = {
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
                        }
                        results.append(result)
                        n_completed += 1

                        # Show progress even for failed probesets
                        if show_progress:
                            if progress_bar:
                                progress_bar.update(1)
                            else:
                                _print_progress(result, n_completed, n_total, verbose_progress)

    # Close progress bar if used
    if progress_bar:
        progress_bar.close()

    # Create results DataFrame
    results_df = pl.DataFrame(results)

    # Sort by probeset_id for consistent output
    results_df = results_df.sort('probeset_id')

    return results_df
