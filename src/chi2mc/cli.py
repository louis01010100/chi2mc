"""
CLI entry point for chi2mc.

Usage:
    chi2mc input.tsv output.tsv [options]
"""

import argparse
import multiprocessing
import time
from datetime import datetime

import polars as pl

from chi2mc.chi2_homogeneity import (
    NUMBA_AVAILABLE,
    TQDM_AVAILABLE,
    run_two_tier_chi2_test,
)


def main():
    """
    Main function to run the two-tier chi-squared test on the data.
    """
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
    parser.add_argument(
        '--show-progress',
        action='store_true',
        help='Show progress during processing (uses tqdm if available)'
    )
    parser.add_argument(
        '--verbose-progress',
        action='store_true',
        help='Show detailed progress with chi-squared and p-values for each probeset'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=1000,
        help='Number of probesets to process in each subprocess (default: 1000)'
    )

    args = parser.parse_args()

    print("Two-Tier Chi-Squared Homogeneity Test (OPTIMIZED)")
    print("=" * 80)
    print(f"Input file: {args.input_file}")
    print(f"Output file: {args.output_file}")

    # Get number of available CPUs
    n_cpus = multiprocessing.cpu_count()
    print(f"Available CPU cores: {n_cpus}")
    n_workers_display = args.n_workers if args.n_workers else n_cpus
    print(f"Processing probesets in parallel using {n_workers_display} worker(s)...")

    # Show optimization status
    if NUMBA_AVAILABLE:
        print("Numba JIT compilation: ENABLED (faster Monte Carlo)")
    else:
        print("Numba JIT compilation: Not installed (still using vectorized operations)")
        print("  Install numba for 5-20x additional speedup: pip install numba")

    # Show progress tracking status
    if args.show_progress or args.verbose_progress:
        if TQDM_AVAILABLE:
            print("Progress tracking: ENABLED (using tqdm progress bar)")
        else:
            print("Progress tracking: ENABLED (text-based)")
            if not args.verbose_progress:
                print("  Install tqdm for better progress bars: pip install tqdm")

    # Show actual number of workers and batch configuration
    actual_workers = args.n_workers if args.n_workers else n_cpus
    print(f"Note: Using {actual_workers} parallel workers for processing")
    print(f"Note: Batch size = {args.batch_size} probesets per worker")
    print()

    # Run the test with parallel processing
    start_time = time.time()
    results = run_two_tier_chi2_test(
        args.input_file,
        min_expected_count=args.min_expected_count,
        n_simulations=args.n_simulations,
        random_seed=args.random_seed,
        n_workers=args.n_workers,
        show_progress=args.show_progress or args.verbose_progress,
        verbose_progress=args.verbose_progress,
        batch_size=args.batch_size
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

    # Count Monte Carlo tests (both regular and numba versions)
    mc_tests = results.filter(pl.col('method').str.contains('monte_carlo'))
    mc_numba_tests = results.filter(pl.col('method') == 'monte_carlo_numba')
    print(f"Monte Carlo tests: {len(mc_tests)}")
    if len(mc_numba_tests) > 0:
        print(f"  - With Numba acceleration: {len(mc_numba_tests)}")
        print(f"  - Without Numba: {len(mc_tests) - len(mc_numba_tests)}")

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
            pct = 100 * n_significant / len(valid_results)
            f.write(f"Significant results (p < 0.05): {n_significant} ({pct:.1f}%)\n\n")

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
                    excl = row['excluded_genotypes']
                    excluded_note = f" (excluded: {excl})" if excl else ""
                    chi2 = row['chi2_statistic']
                    f.write(
                        f"  - {row['probeset_id']}: p={row['pvalue']:.6f},"
                        f" chi2={chi2:.2f}, method={row['method']}{excluded_note}\n"
                    )
                f.write("\n")

        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    print(f"\nDetailed summary report saved to: {summary_file}")
