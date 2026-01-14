"""
Example usage of the two-tier chi-squared homogeneity test.

This script demonstrates how to use the chi2_homogeneity module
with your own data.
"""

from chi2_homogeneity import (
    run_two_tier_chi2_test,
    load_data,
    create_contingency_table,
    check_minimum_expected_frequency,
    chi2_standard_test,
    chi2_monte_carlo_test
)
import polars as pl


def example_1_basic_usage():
    """Example 1: Basic usage with the provided test data."""
    print("=" * 80)
    print("Example 1: Basic Usage")
    print("=" * 80)

    # Run the test on the sample data
    results = run_two_tier_chi2_test(
        'workspace/eval/before.tsv',
        min_expected_count=5,
        n_simulations=10000,
        random_seed=42  # For reproducibility
    )

    # Display results
    print("\nResults:")
    print(results)

    # Filter for significant results
    significant = results.filter(pl.col('pvalue') < 0.05)
    print(f"\n\nSignificant results (p < 0.05): {len(significant)}")
    if len(significant) > 0:
        print(significant.select(['probeset_id', 'chi2_statistic', 'pvalue', 'method']))


def example_2_manual_analysis():
    """Example 2: Manual step-by-step analysis of a single probeset."""
    print("\n\n" + "=" * 80)
    print("Example 2: Manual Step-by-Step Analysis")
    print("=" * 80)

    # Load data
    df = load_data('workspace/eval/before.tsv')
    print(f"\nLoaded data with {len(df)} rows and {len(df['probeset_id'].unique())} unique probesets")

    # Analyze a specific probeset
    probeset_id = 'PS004'
    probeset_data = df.filter(pl.col('probeset_id') == probeset_id)

    print(f"\nAnalyzing probeset: {probeset_id}")
    print("\nData:")
    print(probeset_data)

    # Create contingency table
    table = create_contingency_table(probeset_data)
    print("\nContingency table (batches × genotypes):")
    print(f"{'Batch':<10} {'AA':>8} {'AB':>8} {'BB':>8} {'NC':>8} {'Total':>10}")
    print("-" * 55)
    for i, batch in enumerate(probeset_data['batch_name'].unique().sort()):
        row_total = table[i].sum()
        print(f"{batch:<10} {table[i,0]:>8.0f} {table[i,1]:>8.0f} {table[i,2]:>8.0f} {table[i,3]:>8.0f} {row_total:>10.0f}")
    print("-" * 55)
    print(f"{'Total':<10} {table[:,0].sum():>8.0f} {table[:,1].sum():>8.0f} {table[:,2].sum():>8.0f} {table[:,3].sum():>8.0f} {table.sum():>10.0f}")

    # Check minimum expected frequency
    use_standard = check_minimum_expected_frequency(table, min_count=5)
    print(f"\nAll expected frequencies ≥ 5: {use_standard}")
    print(f"Recommended method: {'Standard Chi-squared' if use_standard else 'Monte Carlo'}")

    # Perform the test
    if use_standard:
        result = chi2_standard_test(table)
    else:
        result = chi2_monte_carlo_test(table, n_simulations=10000, random_seed=42)

    print(f"\nTest Results:")
    print(f"  Method: {result['method']}")
    print(f"  Chi-squared statistic: {result['statistic']:.4f}")
    print(f"  Degrees of freedom: {result['dof']}")
    print(f"  P-value: {result['pvalue']:.6f}")
    if 'n_simulations' in result:
        print(f"  Simulations: {result['n_simulations']}")

    # Interpretation
    alpha = 0.05
    print(f"\nInterpretation (α = {alpha}):")
    if result['pvalue'] < alpha:
        print(f"  ✗ Reject null hypothesis (p = {result['pvalue']:.6f} < {alpha})")
        print(f"  → Genotype frequencies differ significantly across batches")
    else:
        print(f"  ✓ Fail to reject null hypothesis (p = {result['pvalue']:.6f} ≥ {alpha})")
        print(f"  → No significant difference in genotype frequencies across batches")


def example_3_compare_methods():
    """Example 3: Compare standard and Monte Carlo methods on the same data."""
    print("\n\n" + "=" * 80)
    print("Example 3: Comparing Standard vs Monte Carlo Methods")
    print("=" * 80)

    # Load data and get a probeset with large enough samples
    df = load_data('workspace/eval/before.tsv')
    probeset_data = df.filter(pl.col('probeset_id') == 'PS004')
    table = create_contingency_table(probeset_data)

    print("\nTesting probeset PS004 with both methods:")

    # Standard method
    result_standard = chi2_standard_test(table)
    print(f"\nStandard Chi-squared Test:")
    print(f"  Chi-squared: {result_standard['statistic']:.4f}")
    print(f"  P-value: {result_standard['pvalue']:.6f}")

    # Monte Carlo method
    result_mc = chi2_monte_carlo_test(table, n_simulations=10000, random_seed=42)
    print(f"\nMonte Carlo Chi-squared Test (10,000 simulations):")
    print(f"  Chi-squared: {result_mc['statistic']:.4f}")
    print(f"  P-value: {result_mc['pvalue']:.6f}")

    print(f"\nDifference in p-values: {abs(result_standard['pvalue'] - result_mc['pvalue']):.6f}")
    print("\nNote: For large samples, both methods should give similar results.")
    print("Monte Carlo is more accurate for small samples where chi-squared")
    print("approximation may not be valid.")


def example_4_batch_processing():
    """Example 4: Batch processing and result filtering."""
    print("\n\n" + "=" * 80)
    print("Example 4: Batch Processing and Filtering")
    print("=" * 80)

    # Run the test
    results = run_two_tier_chi2_test(
        'workspace/eval/before.tsv',
        min_expected_count=5,
        n_simulations=10000,
        random_seed=42
    )

    print("\nSummary by Method:")
    print(results.group_by('method').count())

    print("\n\nTop 5 most significant results:")
    top_results = results.sort('pvalue').head(5)
    print(top_results.select(['probeset_id', 'chi2_statistic', 'pvalue', 'method']))

    print("\n\nProbesets using Monte Carlo (due to small expected frequencies):")
    mc_results = results.filter(pl.col('method') == 'monte_carlo')
    print(mc_results.select(['probeset_id', 'total_samples', 'pvalue']))

    # Multiple testing correction (Bonferroni)
    n_tests = len(results)
    alpha = 0.05
    bonferroni_threshold = alpha / n_tests
    print(f"\n\nBonferroni correction for multiple testing:")
    print(f"  Number of tests: {n_tests}")
    print(f"  Original α: {alpha}")
    print(f"  Bonferroni-corrected α: {bonferroni_threshold:.6f}")
    significant_bonferroni = results.filter(pl.col('pvalue') < bonferroni_threshold)
    print(f"  Significant after correction: {len(significant_bonferroni)}")


def example_5_custom_data():
    """Example 5: Create and test custom data."""
    print("\n\n" + "=" * 80)
    print("Example 5: Testing Custom Data")
    print("=" * 80)

    # Create custom data showing batch effect
    custom_data = pl.DataFrame({
        'probeset_id': ['CUSTOM1'] * 3,
        'batch_name': ['batchA', 'batchB', 'batchC'],
        'n_AA': [100, 80, 110],  # Similar proportions
        'n_AB': [150, 140, 145],
        'n_BB': [50, 60, 45],
        'n_NC': [10, 8, 12]
    })

    print("\nCustom data (no batch effect):")
    print(custom_data)

    # Save to temporary file
    import tempfile
    import os
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        custom_data.write_csv(f.name, separator='\t')
        temp_file = f.name

    # Test it
    results = run_two_tier_chi2_test(temp_file, random_seed=42)
    print("\nResults:")
    print(results.select(['probeset_id', 'chi2_statistic', 'pvalue', 'method']))

    # Clean up
    os.unlink(temp_file)

    # Create data with obvious batch effect
    custom_data_effect = pl.DataFrame({
        'probeset_id': ['CUSTOM2'] * 3,
        'batch_name': ['batchA', 'batchB', 'batchC'],
        'n_AA': [200, 50, 100],   # Very different proportions
        'n_AB': [50, 200, 100],
        'n_BB': [50, 50, 100],
        'n_NC': [10, 10, 10]
    })

    print("\n\nCustom data (with batch effect):")
    print(custom_data_effect)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        custom_data_effect.write_csv(f.name, separator='\t')
        temp_file = f.name

    results_effect = run_two_tier_chi2_test(temp_file, random_seed=42)
    print("\nResults:")
    print(results_effect.select(['probeset_id', 'chi2_statistic', 'pvalue', 'method']))

    os.unlink(temp_file)


if __name__ == '__main__':
    """Run all examples."""

    print("\n")
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 15 + "Two-Tier Chi-Squared Homogeneity Test Examples" + " " * 16 + "║")
    print("╚" + "=" * 78 + "╝")

    # Run all examples
    example_1_basic_usage()
    example_2_manual_analysis()
    example_3_compare_methods()
    example_4_batch_processing()
    example_5_custom_data()

    print("\n\n" + "=" * 80)
    print("All examples completed!")
    print("=" * 80)
