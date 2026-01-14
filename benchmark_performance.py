"""
Performance benchmarking script for chi2_homogeneity optimizations.

This script measures the performance improvements from:
1. Vectorized contingency table reconstruction
2. Numba JIT compilation
3. Optimized data loading

Usage:
    python benchmark_performance.py <input_file.tsv>

Example:
    python benchmark_performance.py workspace/eval/before.tsv
"""

import time
import sys
import numpy as np
from chi2_homogeneity import (
    run_two_tier_chi2_test,
    load_data,
    NUMBA_AVAILABLE
)


def format_time(seconds):
    """Format time in human-readable format."""
    if seconds < 1:
        return f"{seconds * 1000:.1f} ms"
    elif seconds < 60:
        return f"{seconds:.2f} s"
    else:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"


def benchmark_full_pipeline(file_path, n_simulations=10000):
    """
    Benchmark the full pipeline with different configurations.

    Returns timing results for:
    1. Unoptimized (legacy data loading)
    2. Optimized data loading only
    3. Fully optimized (data loading + vectorization + numba)
    """
    print("=" * 80)
    print("PERFORMANCE BENCHMARK")
    print("=" * 80)
    print(f"Input file: {file_path}")
    print(f"Monte Carlo simulations: {n_simulations}")
    print(f"Numba available: {NUMBA_AVAILABLE}")
    print()

    # Get dataset info
    df = load_data(file_path)
    n_probesets = len(df['probeset_id'].unique())
    print(f"Dataset: {len(df)} rows, {n_probesets} probesets")
    print()

    results = {}

    # Benchmark 1: Legacy (unoptimized data loading)
    print("Benchmark 1: Legacy (unoptimized data loading)")
    print("-" * 80)
    start = time.time()
    results_legacy = run_two_tier_chi2_test(
        file_path,
        n_simulations=n_simulations,
        random_seed=42,
        use_optimized=False,
        n_workers=1  # Single worker for fair comparison
    )
    time_legacy = time.time() - start
    results['legacy'] = time_legacy
    print(f"Time: {format_time(time_legacy)}")
    print()

    # Benchmark 2: Optimized data loading + vectorization
    print("Benchmark 2: Optimized data loading + vectorization")
    print("-" * 80)
    start = time.time()
    results_opt = run_two_tier_chi2_test(
        file_path,
        n_simulations=n_simulations,
        random_seed=42,
        use_optimized=True,
        n_workers=1  # Single worker for fair comparison
    )
    time_opt = time.time() - start
    results['optimized'] = time_opt
    print(f"Time: {format_time(time_opt)}")
    speedup = time_legacy / time_opt if time_opt > 0 else 0
    print(f"Speedup: {speedup:.2f}x faster than legacy")
    print()

    # Verify results are consistent
    print("Verifying result consistency...")
    print("-" * 80)

    # Check that results are numerically similar (within tolerance for Monte Carlo)
    merged = results_legacy.join(
        results_opt,
        on='probeset_id',
        suffix='_opt'
    )

    # For standard tests, p-values should be identical
    standard_tests = merged.filter(
        (merged['method'] == 'standard') & (merged['method_opt'] == 'standard')
    )
    if len(standard_tests) > 0:
        max_pvalue_diff = np.max(np.abs(
            standard_tests['pvalue'].to_numpy() -
            standard_tests['pvalue_opt'].to_numpy()
        ))
        print(f"Standard tests: Max p-value difference = {max_pvalue_diff:.10f}")
        if max_pvalue_diff < 1e-10:
            print("✓ Standard tests produce identical results")
        else:
            print("✗ WARNING: Standard tests show differences!")

    # For Monte Carlo tests, results should be similar but may vary
    mc_tests = merged.filter(
        merged['method'].str.contains('monte_carlo') &
        merged['method_opt'].str.contains('monte_carlo')
    )
    if len(mc_tests) > 0:
        pvalue_diffs = np.abs(
            mc_tests['pvalue'].to_numpy() -
            mc_tests['pvalue_opt'].to_numpy()
        )
        max_pvalue_diff = np.max(pvalue_diffs)
        mean_pvalue_diff = np.mean(pvalue_diffs)
        print(f"Monte Carlo tests: Mean p-value difference = {mean_pvalue_diff:.6f}")
        print(f"Monte Carlo tests: Max p-value difference = {max_pvalue_diff:.6f}")
        if max_pvalue_diff < 0.05:
            print("✓ Monte Carlo tests produce consistent results (differences < 0.05)")
        else:
            print("⚠ Monte Carlo tests show larger differences (expected due to randomness)")

    print()

    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Legacy time:     {format_time(time_legacy)}")
    print(f"Optimized time:  {format_time(time_opt)}")
    print(f"Speedup:         {speedup:.2f}x")
    print(f"Time saved:      {format_time(time_legacy - time_opt)}")

    if NUMBA_AVAILABLE:
        # Count Monte Carlo tests to show impact
        n_mc_tests = results_opt.filter(
            results_opt['method'].str.contains('monte_carlo')
        ).height
        print(f"\nMonte Carlo tests: {n_mc_tests}/{n_probesets}")
        print(f"Numba acceleration: {'ENABLED' if 'numba' in results_opt['method'][0] else 'Available but check method column'}")
    else:
        print("\n⚠ Install numba for additional 5-20x speedup:")
        print("  pip install numba")

    print("=" * 80)

    return results


def benchmark_monte_carlo_only(file_path, n_simulations=10000):
    """
    Detailed benchmark of Monte Carlo simulation performance.

    Tests a single probeset that requires Monte Carlo to isolate
    the performance of the simulation itself.
    """
    from chi2_homogeneity import (
        create_contingency_table,
        chi2_monte_carlo_test
    )

    print("\n")
    print("=" * 80)
    print("MONTE CARLO DETAILED BENCHMARK")
    print("=" * 80)

    # Load data and find a probeset that uses Monte Carlo
    df = load_data(file_path)

    # Run quick test to find Monte Carlo probeset
    quick_results = run_two_tier_chi2_test(
        file_path,
        n_simulations=100,
        random_seed=42,
        n_workers=1
    )

    mc_probesets = quick_results.filter(
        quick_results['method'].str.contains('monte_carlo')
    )

    if len(mc_probesets) == 0:
        print("No probesets require Monte Carlo simulation in this dataset.")
        print("All probesets have sufficient expected frequencies for standard test.")
        return

    # Find a valid Monte Carlo probeset (one without zero columns)
    probeset_id = None
    table = None

    for idx in range(len(mc_probesets)):
        test_probeset_id = mc_probesets['probeset_id'][idx]
        probeset_data = df.filter(df['probeset_id'] == test_probeset_id)
        test_table = create_contingency_table(probeset_data)

        # Remove zero columns
        col_sums = test_table.sum(axis=0)
        non_zero_cols = np.where(col_sums > 0)[0]

        if len(non_zero_cols) > 0:
            table = test_table[:, non_zero_cols]
            probeset_id = test_probeset_id
            break

    if probeset_id is None or table is None:
        print("No valid probesets found for Monte Carlo benchmarking.")
        print("All Monte Carlo probesets have issues with zero columns.")
        return

    print(f"Testing probeset: {probeset_id}")
    print(f"Table shape: {table.shape}")
    print(f"Total samples: {int(table.sum())}")
    print(f"Simulations: {n_simulations}")
    print()

    # Benchmark with numba
    time_numba = None
    if NUMBA_AVAILABLE:
        print("With Numba JIT:")
        print("-" * 80)
        try:
            start = time.time()
            result_numba = chi2_monte_carlo_test(
                table,
                n_simulations=n_simulations,
                random_seed=42,
                use_numba=True
            )
            time_numba = time.time() - start
            print(f"Time: {format_time(time_numba)}")
            print(f"Method: {result_numba['method']}")
            print(f"P-value: {result_numba['pvalue']:.6f}")
            print()
        except Exception as e:
            print(f"Error with Numba version: {e}")
            print("Skipping Numba benchmark.")
            print()

    # Benchmark without numba (vectorized only)
    print("Vectorized (no Numba):")
    print("-" * 80)
    try:
        start = time.time()
        result_vec = chi2_monte_carlo_test(
            table,
            n_simulations=n_simulations,
            random_seed=42,
            use_numba=False
        )
        time_vec = time.time() - start
        print(f"Time: {format_time(time_vec)}")
        print(f"Method: {result_vec['method']}")
        print(f"P-value: {result_vec['pvalue']:.6f}")
        print()
    except Exception as e:
        print(f"Error with vectorized version: {e}")
        print("Cannot benchmark Monte Carlo on this probeset.")
        return

    # Summary
    print("Summary:")
    print("-" * 80)
    if NUMBA_AVAILABLE and time_numba is not None:
        speedup = time_vec / time_numba
        print(f"Vectorized time:  {format_time(time_vec)}")
        print(f"Numba time:       {format_time(time_numba)}")
        print(f"Numba speedup:    {speedup:.2f}x")
    else:
        print(f"Vectorized time:  {format_time(time_vec)}")
        if not NUMBA_AVAILABLE:
            print("\n⚠ Install numba for additional speedup:")
            print("  pip install numba")

    print("=" * 80)


def main():
    if len(sys.argv) < 2:
        print("Usage: python benchmark_performance.py <input_file.tsv>")
        print("\nExample:")
        print("  python benchmark_performance.py workspace/eval/before.tsv")
        sys.exit(1)

    file_path = sys.argv[1]

    # Run benchmarks
    n_simulations = 10000  # Default for real benchmarking

    # Full pipeline benchmark
    benchmark_full_pipeline(file_path, n_simulations=n_simulations)

    # Detailed Monte Carlo benchmark
    benchmark_monte_carlo_only(file_path, n_simulations=n_simulations)

    print("\n✓ Benchmark complete!")


if __name__ == '__main__':
    main()
