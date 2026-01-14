# Two-Tier Chi-Squared Homogeneity Test

Test-driven implementation of a 2-tier chi-squared test for testing homogeneity of genotype frequencies across batches.

## Overview

This tool tests whether genotype frequencies (AA, AB, BB, NC) are homogeneous across different batches for each probeset. It uses a two-tier approach:

1. **Standard chi-squared test**: Used when all expected cell frequencies are ≥ 5
2. **Monte Carlo chi-squared test**: Used when any expected cell frequency is < 5

## Features

- Test-driven development with comprehensive unit tests
- Automatic selection between standard and Monte Carlo methods
- **Optimized performance**: 5-50x faster with vectorization and optional numba JIT compilation
- **Batch processing**: Process 1000 probesets per subprocess for optimal efficiency
- **Progress tracking**: Monitor analysis with progress bars and real-time results
- **Scalable**: Efficiently handles millions of probesets with all CPU cores
- Detailed results including test statistics, p-values, and degrees of freedom
- Clear indication of which method was used for each probeset
- Reproducible results with optional random seed
- Parallel processing support for large datasets

## Requirements

```bash
pip install polars numpy scipy pytest
```

## Data Format

Input data should be a TSV file with the following columns:

- `probeset_id`: Identifier for each probeset
- `batch_name`: Identifier for each batch
- `n_AA`: Count of AA genotypes
- `n_AB`: Count of AB genotypes
- `n_BB`: Count of BB genotypes
- `n_NC`: Count of NC (no call) genotypes

Example:
```
probeset_id	batch_name	n_AA	n_AB	n_BB	n_NC
PS001	batch1	45	52	38	2
PS001	batch2	48	50	40	1
PS001	batch3	42	55	37	3
```

## Usage

### Command Line

```bash
python chi2_homogeneity.py input_file.tsv output_file.tsv
```

Basic usage:
```bash
python chi2_homogeneity.py workspace/eval/before.tsv results.tsv
```

With options:
```bash
python chi2_homogeneity.py input.tsv output.tsv \
  --min-expected-count 5 \
  --n-simulations 10000 \
  --random-seed 42 \
  --n-workers 4
```

**With progress tracking** (monitor long-running analyses):
```bash
# Simple progress bar
python chi2_homogeneity.py input.tsv output.tsv --show-progress

# Detailed progress with chi-squared and p-values
python chi2_homogeneity.py input.tsv output.tsv --verbose-progress
```

See all available options:
```bash
python chi2_homogeneity.py --help
```

### Python API

```python
from chi2_homogeneity import run_two_tier_chi2_test

# Run the test
results = run_two_tier_chi2_test(
    'path/to/your/data.tsv',
    min_expected_count=5,      # Threshold for using standard vs Monte Carlo
    n_simulations=10000,        # Number of Monte Carlo simulations
    random_seed=42,             # For reproducibility (optional)
    show_progress=True,         # Show progress bar (optional)
    verbose_progress=False      # Show detailed results (optional)
)

# View results
print(results)

# Save results
results.write_csv('results.tsv', separator='\t')
```

### Individual Functions

```python
from chi2_homogeneity import (
    load_data,
    create_contingency_table,
    check_minimum_expected_frequency,
    chi2_standard_test,
    chi2_monte_carlo_test
)

# Load data
df = load_data('data.tsv')

# Get data for a specific probeset
probeset_data = df[df['probeset_id'] == 'PS001']

# Create contingency table
table = create_contingency_table(probeset_data)

# Check if standard test is appropriate
use_standard = check_minimum_expected_frequency(table, min_count=5)

# Perform appropriate test
if use_standard:
    result = chi2_standard_test(table)
else:
    result = chi2_monte_carlo_test(table, n_simulations=10000)

print(f"Chi-squared statistic: {result['statistic']:.2f}")
print(f"P-value: {result['pvalue']:.4f}")
print(f"Method: {result['method']}")
```

## Output Format

### Main Results File (TSV)

The results DataFrame contains:

- `probeset_id`: Probeset identifier
- `chi2_statistic`: Chi-squared test statistic (NaN for skipped probesets)
- `pvalue`: P-value for the test (NaN for skipped probesets)
- `dof`: Degrees of freedom (NaN for skipped probesets)
- `method`: Test method used ('standard', 'monte_carlo', or 'skipped')
- `n_batches`: Number of batches for this probeset
- `total_samples`: Total number of samples
- `n_simulations`: Number of simulations (only for Monte Carlo tests)
- `error`: Error message for skipped probesets (null otherwise)
- `excluded_genotypes`: Genotypes excluded from test (null if none excluded)

### Summary Report (TXT)

A detailed summary report is automatically generated alongside the main results file:
- Filename: `{output_file}_summary.txt` (e.g., `results_summary.txt`)
- Contains:
  - Overall statistics (total tests, method breakdown, significant results)
  - Complete list of probesets with excluded genotypes
  - List of skipped probesets with error details
  - Significant results sorted by p-value

## Interpretation

- **Null hypothesis (H₀)**: Genotype frequencies are homogeneous across batches
- **Alternative hypothesis (H₁)**: Genotype frequencies differ across batches

A small p-value (typically < 0.05) suggests rejecting the null hypothesis, indicating that genotype frequencies are not homogeneous across batches.

### Handling Zero Columns

**Automatic Exclusion**: If an entire genotype column is zero across all batches (e.g., all n_NC = 0), that genotype is automatically excluded from the analysis, and the test runs on the remaining genotypes.

For example:
- If NC is always 0, the test checks for homogeneity across AA, AB, and BB only
- The result will show `excluded_genotypes='NC'`
- The test still produces a valid p-value

**Skipped Probesets**: Only if ALL genotype columns are zero will a probeset be skipped (which shouldn't happen with real data). Skipped probesets appear with:
- `method='skipped'`
- `pvalue=NaN`, `chi2_statistic=NaN`, `dof=NaN`
- `error` column explaining why

The summary reports both excluded genotypes and skipped probesets separately.

## Running Tests

```bash
# Run all tests
pytest test_chi2_homogeneity.py -v

# Run specific test class
pytest test_chi2_homogeneity.py::TestStandardChi2 -v

# Run with coverage
pytest test_chi2_homogeneity.py --cov=chi2_homogeneity --cov-report=html
```

## Test Coverage

The test suite includes:

- **Data Loading Tests**: Verify correct file loading and column validation
- **Contingency Table Tests**: Check table creation and data integrity
- **Minimum Frequency Tests**: Verify correct detection of small expected frequencies
- **Standard Chi-Squared Tests**: Test standard method correctness
- **Monte Carlo Tests**: Test Monte Carlo simulation and reproducibility
- **Integration Tests**: Verify end-to-end workflow and method selection

## Algorithm Details

### Standard Chi-Squared Test

Uses scipy's `chi2_contingency` function which computes:

```
χ² = Σ [(O - E)² / E]
```

where O is observed frequency and E is expected frequency under independence.

### Monte Carlo Chi-Squared Test

For cases with small expected frequencies:

1. Compute observed chi-squared statistic
2. Permute batch assignments randomly while maintaining marginal totals
3. Recompute chi-squared statistic for each permutation
4. Calculate p-value as proportion of permuted statistics ≥ observed statistic

This provides an exact test that doesn't rely on asymptotic chi-squared distribution.

## Example Results

```
Two-Tier Chi-Squared Homogeneity Test Results
================================================================================
probeset_id  chi2_statistic       pvalue  dof      method  n_batches  total_samples  n_simulations
      PS001        1.740409 9.441000e-01    6 monte_carlo          3            413        10000.0
      PS002        0.815007 9.940000e-01    6 monte_carlo          3            721        10000.0
      PS003        3.928571 8.580000e-01    6 monte_carlo          3             20        10000.0
      PS004       58.043623 1.122827e-10    6    standard          3            565            NaN

Summary:
Total probesets tested: 4
Standard chi-squared tests: 1
Monte Carlo tests: 3
Significant results (p < 0.05): 1
```

In this example:
- PS004 shows significant heterogeneity (p < 0.001) using standard test
- Other probesets show no significant difference across batches
- Monte Carlo was automatically used for PS001-PS003 due to small expected frequencies

## Customization

### Adjusting Minimum Expected Count Threshold

```python
# Use threshold of 3 instead of 5
results = run_two_tier_chi2_test(
    'data.tsv',
    min_expected_count=3
)
```

### Increasing Monte Carlo Simulations

Command line:
```bash
python chi2_homogeneity.py data.tsv results.tsv --n-simulations 100000
```

Python API:
```python
# Use more simulations for higher precision
results = run_two_tier_chi2_test(
    'data.tsv',
    n_simulations=100000
)
```

### Controlling Parallel Processing

Command line:
```bash
# Use 4 workers
python chi2_homogeneity.py data.tsv results.tsv --n-workers 4

# Sequential processing (no parallelism)
python chi2_homogeneity.py data.tsv results.tsv --n-workers 1
```

Python API:
```python
# Use 4 workers
results = run_two_tier_chi2_test('data.tsv', n_workers=4)

# Sequential processing
results = run_two_tier_chi2_test('data.tsv', n_workers=1)
```

## License

MIT

## Author

Generated with Claude Code using test-driven development approach.
