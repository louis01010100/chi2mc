# Two-Tier Chi-Squared Homogeneity Test

A fast, native Rust implementation of a 2-tier chi-squared test for testing homogeneity of genotype frequencies across batches.

## Overview

This tool tests whether genotype frequencies (AA, AB, BB, NC) are homogeneous across different batches for each probeset. It uses a two-tier approach:

1. **Standard chi-squared test**: Used when all expected cell frequencies are >= 5
2. **Monte Carlo chi-squared test**: Used when any expected cell frequency is < 5

## Features

- Native Rust binary for maximum performance
- Automatic selection between standard and Monte Carlo methods
- Parallel processing with [rayon](https://docs.rs/rayon) across all CPU cores
- Progress tracking with [indicatif](https://docs.rs/indicatif) progress bars
- Detailed results including test statistics, p-values, and degrees of freedom
- Clear indication of which method was used for each probeset
- Reproducible results with configurable random seed
- Summary report generation

## Building

Requires [Rust](https://www.rust-lang.org/tools/install) (edition 2021).

```bash
cargo build --release
```

The binary will be at `target/release/chi2mc`.

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

Basic usage:
```bash
chi2mc input_file.tsv output_file.tsv
```

With options:
```bash
chi2mc input.tsv output.tsv \
  --min-expected-count 5 \
  --n-simulations 10000 \
  --random-seed 42 \
  --n-workers 4
```

With progress tracking:
```bash
# Simple progress bar
chi2mc input.tsv output.tsv --show-progress

# Detailed progress with chi-squared and p-values
chi2mc input.tsv output.tsv --verbose-progress
```

See all available options:
```bash
chi2mc --help
```

### CLI Options

| Option | Default | Description |
|---|---|---|
| `--min-expected-count` | 5 | Threshold for standard vs Monte Carlo |
| `--n-simulations` | 10000 | Number of Monte Carlo simulations |
| `--random-seed` | 42 | Random seed for reproducibility |
| `--n-workers` | all cores | Number of parallel workers |
| `--show-progress` | off | Show progress bar |
| `--verbose-progress` | off | Show per-probeset results during processing |
| `--batch-size` | 1000 | Number of probesets per worker batch |

## Output Format

### Main Results File (TSV)

The results file contains:

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

- **Null hypothesis (H0)**: Genotype frequencies are homogeneous across batches
- **Alternative hypothesis (H1)**: Genotype frequencies differ across batches

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
cargo test

# Run with output
cargo test -- --nocapture

# Run integration tests only
cargo test --test integration
```

## Project Structure

```
src/
  main.rs          # Entry point, CLI orchestration, parallel processing
  cli.rs           # Command-line argument parsing (clap)
  io.rs            # TSV input/output
  stats.rs         # Chi-squared statistic computation
  monte_carlo.rs   # Monte Carlo simulation engine
  probeset.rs      # Per-probeset processing logic
  report.rs        # Summary report generation
tests/
  integration.rs   # End-to-end integration tests
```

## Algorithm Details

### Standard Chi-Squared Test

Computes the Pearson chi-squared statistic:

```
chi2 = sum[(O - E)^2 / E]
```

where O is observed frequency and E is expected frequency under independence. The p-value is computed from the chi-squared CDF using the [statrs](https://docs.rs/statrs) crate.

### Monte Carlo Chi-Squared Test

For cases with small expected frequencies:

1. Compute observed chi-squared statistic
2. Permute batch assignments randomly while maintaining marginal totals
3. Recompute chi-squared statistic for each permutation
4. Calculate p-value as proportion of permuted statistics >= observed statistic

This provides an exact test that doesn't rely on the asymptotic chi-squared distribution. Random number generation uses [rand](https://docs.rs/rand) with the ChaCha8 PRNG for reproducibility.

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

## License

MIT

## Author

Generated with Claude Code using test-driven development approach.
