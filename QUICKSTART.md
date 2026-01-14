# Quick Start Guide

## Using with Your Data

Run the analysis on your data:

```bash
python chi2_homogeneity.py input_file.tsv output_file.tsv
```

For example:

```bash
python chi2_homogeneity.py workspace/eval/before.tsv results.tsv
```

## File Format

Your TSV file must have these exact column names:
- `probeset_id`
- `batch_name`
- `n_AA`
- `n_AB`
- `n_BB`
- `n_NC`

## Quick Commands

```bash
# Run the analysis with default settings
python chi2_homogeneity.py input.tsv output.tsv

# Run with custom settings
python chi2_homogeneity.py input.tsv output.tsv --n-simulations 50000 --n-workers 4

# See all available options
python chi2_homogeneity.py --help

# Run all tests to verify installation
pytest test_chi2_homogeneity.py -v

# See examples
python example_usage.py
```

## Python API

```python
from chi2_homogeneity import run_two_tier_chi2_test
import polars as pl

# Analyze your data
results = run_two_tier_chi2_test('your_data.tsv')

# Save results
results.write_csv('results.tsv', separator='\t')

# Find significant probesets (p < 0.05)
significant = results.filter(pl.col('pvalue') < 0.05)
print(significant)
```

## Output Files

Running the script generates two files:

1. **Results file** (TSV): Main results with one row per probeset
2. **Summary report** (TXT): Detailed summary with all excluded genotypes and significant results

For example, running:
```bash
python chi2_homogeneity.py input.tsv results.tsv
```

Produces:
- `results.tsv` - Main results file
- `results_summary.txt` - Detailed summary report

## Interpreting Results

The main results file includes these columns:

| Column | Description |
|--------|-------------|
| `probeset_id` | Your probeset identifier |
| `chi2_statistic` | Test statistic (higher = more difference) |
| `pvalue` | Significance level (< 0.05 = significant) |
| `dof` | Degrees of freedom |
| `method` | 'standard', 'monte_carlo', or 'skipped' |
| `n_batches` | Number of batches tested |
| `total_samples` | Total sample count |
| `n_simulations` | Number of simulations (Monte Carlo only) |
| `error` | Error message (only for skipped probesets) |
| `excluded_genotypes` | Genotypes excluded from test (e.g., 'NC' if all zeros) |

### What the p-value means:

- **p < 0.05**: Genotype frequencies differ significantly across batches
- **p ≥ 0.05**: No significant difference in genotype frequencies

### Which method is used:

- **Standard**: All expected cell frequencies are ≥ 5 (more common for large samples)
- **Monte Carlo**: Some expected frequencies are < 5 (used for small samples or sparse data)
- **Skipped**: Chi-squared test cannot be computed (e.g., entire genotype column is zero). See `error` column for details.

## Adjusting Parameters

### Command Line

```bash
# Use different threshold for method selection
python chi2_homogeneity.py data.tsv results.tsv --min-expected-count 3

# More Monte Carlo simulations for better precision
python chi2_homogeneity.py data.tsv results.tsv --n-simulations 100000

# Set random seed for reproducibility
python chi2_homogeneity.py data.tsv results.tsv --random-seed 12345

# Control number of parallel workers
python chi2_homogeneity.py data.tsv results.tsv --n-workers 4
```

### Python API

```python
# Use different threshold for method selection
results = run_two_tier_chi2_test(
    'data.tsv',
    min_expected_count=3  # Use Monte Carlo if any expected count < 3
)

# More Monte Carlo simulations for better precision
results = run_two_tier_chi2_test(
    'data.tsv',
    n_simulations=100000  # Default is 10,000
)

# Set random seed for reproducibility
results = run_two_tier_chi2_test(
    'data.tsv',
    random_seed=42
)
```

## Common Issues

### Missing columns error
**Error**: `Missing required columns: ['column_name']`
**Solution**: Check your TSV has all required columns with exact names

### All zeros in a column
**Behavior**: If an entire genotype column (e.g., all n_NC = 0) is zero across all batches, that column is automatically excluded from the test.

**Output**: The test runs on the remaining genotypes (e.g., AA, AB, BB only), and the result shows `excluded_genotypes='NC'`

**Example**: If NC is always 0, the test checks for homogeneity across AA, AB, and BB only.

**Note**: Only if ALL genotype columns are zero will the probeset be skipped (which shouldn't happen with real data).

### Large p-values everywhere
If all p-values are close to 1.0, this suggests:
- Genotype frequencies are very consistent across batches (good!)
- No batch effects detected

### Very small p-values
If many p-values are near 0, this suggests:
- Significant batch effects in genotype calling
- May need to investigate batch-specific issues

## Multiple Testing Correction

When testing many probesets, apply Bonferroni correction:

```python
import polars as pl
results = run_two_tier_chi2_test('data.tsv')

# Bonferroni correction
alpha = 0.05
n_tests = len(results)
bonferroni_threshold = alpha / n_tests

# Filter using corrected threshold
significant = results.filter(pl.col('pvalue') < bonferroni_threshold)
print(f"Significant after Bonferroni correction: {len(significant)}")
```

## Getting Help

- See `README.md` for detailed documentation
- Run `python example_usage.py` for comprehensive examples
- All code is test-driven - see `test_chi2_homogeneity.py` for usage patterns

## Citation

If using this for publications, you may want to cite:
- Standard chi-squared test: Pearson (1900)
- Monte Carlo approach: Hope (1968) or Patefield (1981) for contingency table generation
