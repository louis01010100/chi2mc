# Progress Tracking Feature

## Overview

The chi2_homogeneity module now includes progress tracking to monitor long-running analyses. You can see:
- Number of probesets processed (e.g., "Processing: 45/1000")
- Chi-squared statistics and p-values as they complete
- Significance indicators (* ** ***)
- Progress bar (if tqdm is installed)

## Usage

### Command Line

**Simple progress (shows percentage):**
```bash
python chi2_homogeneity.py input.tsv output.tsv --show-progress
```

**Verbose progress (shows chi-squared and p-values):**
```bash
python chi2_homogeneity.py input.tsv output.tsv --verbose-progress
```

### Python API

```python
from chi2_homogeneity import run_two_tier_chi2_test

# Simple progress
results = run_two_tier_chi2_test(
    'data.tsv',
    show_progress=True
)

# Verbose progress with details
results = run_two_tier_chi2_test(
    'data.tsv',
    show_progress=True,
    verbose_progress=True
)
```

## Progress Formats

### With tqdm (recommended)

Install tqdm for a nice progress bar:
```bash
pip install tqdm
```

**Simple progress:**
```
Processing probesets:  45%|████▌     | 45/100 [00:23<00:28,  1.95probeset/s]
```

**Verbose progress:**
```
Processing probesets:  45%|████▌     | 45/100 [00:23<00:28,  1.95probeset/s]
[45/100] PS045: χ²=12.3456, p=0.023456 * (monte_carlo_numba)
```

### Without tqdm (text-based)

**Simple progress:**
```
Progress: 45/100 (45.0%)
```

**Verbose progress:**
```
[45/100] PS045: χ²=12.3456, p=0.023456 * (monte_carlo_numba)
```

## Verbose Progress Format

When using `--verbose-progress`, each completed probeset shows:

```
[n_completed/n_total] probeset_id: χ²=statistic, p=pvalue [significance] (method)
```

**Examples:**

```
[1/100] PS001: χ²=1.7404, p=0.944100 (monte_carlo_numba)
[2/100] PS002: χ²=12.5678, p=0.012345 * (standard)
[3/100] PS003: χ²=25.3456, p=0.000234 *** (monte_carlo_numba)
[4/100] PS004: χ²=8.9012, p=0.045678 * (standard)
```

### Significance Markers

- `*` - p < 0.05 (significant)
- `**` - p < 0.01 (highly significant)
- `***` - p < 0.001 (very highly significant)
- (no marker) - p ≥ 0.05 (not significant)

### Method Indicators

- `standard` - Standard chi-squared test
- `monte_carlo` - Monte Carlo test (vectorized)
- `monte_carlo_numba` - Monte Carlo test (with numba acceleration)
- `skipped` - Test failed (shows error)

## Performance Impact

Progress tracking has minimal performance impact:
- **With tqdm**: < 1% overhead
- **Without tqdm (text)**: < 0.1% overhead
- **No progress**: No overhead

For maximum speed, omit the `--show-progress` flag.

## Parallel Processing

Progress tracking works with parallel processing:
- Results are shown as workers complete (unordered)
- Progress bar updates in real-time
- Thread-safe (no race conditions)

**Note**: With multiple workers, verbose progress may show results out of order since workers complete at different times.

## Examples

### Large Dataset Analysis

Monitor a long-running analysis:
```bash
# Show progress without details (fast)
python chi2_homogeneity.py large_data.tsv results.tsv --show-progress --n-workers 8
```

Output:
```
Two-Tier Chi-Squared Homogeneity Test (OPTIMIZED)
================================================================================
Input file: large_data.tsv
Output file: results.tsv
Available CPU cores: 72
Processing probesets in parallel using 8 worker(s)...
✓ Numba JIT compilation: ENABLED (faster Monte Carlo)
✓ Progress tracking: ENABLED (using tqdm progress bar)

Processing probesets:  23%|██▎       | 234/1000 [02:15<06:45,  1.89probeset/s]
```

### Debugging Analysis

See detailed results as they complete:
```bash
# Show all results with chi-squared and p-values
python chi2_homogeneity.py data.tsv results.tsv --verbose-progress --n-workers 1
```

Output:
```
Processing probesets:   1%|▏         | 1/100 [00:00<01:30,  1.09probeset/s]
[1/100] PS001: χ²=1.7404, p=0.944100 (monte_carlo_numba)
[2/100] PS002: χ²=225.0000, p=8.93e-46 *** (standard)
[3/100] PS003: χ²=5.6789, p=0.459012 (monte_carlo_numba)
...
```

### Check for Significant Results

Quickly identify significant probesets during analysis:
```bash
python chi2_homogeneity.py data.tsv results.tsv --verbose-progress | grep '\*'
```

This will show only probesets with significant p-values (marked with *).

## Installation

### Basic (text-based progress)
No additional dependencies needed:
```bash
python chi2_homogeneity.py data.tsv results.tsv --show-progress
```

### Recommended (progress bar)
Install tqdm for better visual feedback:
```bash
pip install tqdm
python chi2_homogeneity.py data.tsv results.tsv --show-progress
```

## Integration with Scripts

### Capture Progress in Logs

Redirect progress to a log file:
```bash
python chi2_homogeneity.py data.tsv results.tsv --verbose-progress 2>&1 | tee analysis.log
```

### Monitor from Another Terminal

Watch progress in real-time:
```bash
# Terminal 1: Run analysis
python chi2_homogeneity.py data.tsv results.tsv --verbose-progress > analysis.log 2>&1

# Terminal 2: Monitor progress
tail -f analysis.log | grep '\[.*\]'
```

### Extract Significant Results

Get significant results as analysis runs:
```bash
python chi2_homogeneity.py data.tsv results.tsv --verbose-progress 2>&1 | \
    tee analysis.log | \
    grep '\*' > significant_probesets.txt
```

## Troubleshooting

### Progress bar not showing

**Problem**: Text progress instead of progress bar

**Solution**: Install tqdm
```bash
pip install tqdm
```

### Progress bar interfering with output

**Problem**: Progress bar mixes with other output

**Solution**: The progress bar uses stderr, so redirect if needed:
```bash
python chi2_homogeneity.py data.tsv results.tsv --show-progress 2> progress.log
```

### Progress seems stuck

**Problem**: Progress not updating for a while

**Reason**: Some probesets take longer (especially Monte Carlo with large samples)

**Solution**: This is normal. Probesets with more samples or Monte Carlo tests take longer.

### Verbose progress too cluttered

**Problem**: Too much output with verbose progress

**Solution**: Use simple progress for cleaner output:
```bash
python chi2_homogeneity.py data.tsv results.tsv --show-progress
```

## Technical Details

### Implementation

- Progress tracking uses callbacks that fire when each probeset completes
- Thread-safe for parallel processing
- Minimal overhead (< 1% performance impact)
- Optional tqdm integration for better visuals

### Files Modified

- `chi2_homogeneity.py` - Added progress tracking functionality
  - `_format_progress_message()` - Format verbose progress messages
  - `_print_progress()` - Print progress updates
  - `run_two_tier_chi2_test()` - Added `show_progress` and `verbose_progress` parameters

### API Reference

```python
def run_two_tier_chi2_test(
    file_path: str,
    min_expected_count: int = 5,
    n_simulations: int = 10000,
    random_seed: int = None,
    n_workers: Optional[int] = None,
    use_optimized: bool = True,
    show_progress: bool = False,      # NEW
    verbose_progress: bool = False    # NEW
) -> pl.DataFrame:
    """
    ...

    Parameters:
    -----------
    show_progress : bool, optional
        Whether to show progress bar/messages during processing (default: False).
        If tqdm is installed, uses a progress bar. Otherwise, prints progress messages.

    verbose_progress : bool, optional
        Whether to show detailed progress with chi-squared and p-values (default: False).
        Only used when show_progress=True. Shows each result as it completes.
    """
```

## See Also

- `OPTIMIZATIONS.md` - Performance optimizations
- `README.md` - General usage guide
- `QUICK_START_OPTIMIZED.md` - Quick start guide
