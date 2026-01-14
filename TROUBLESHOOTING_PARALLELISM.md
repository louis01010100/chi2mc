# Troubleshooting Parallelism and CPU Usage

## Verifying Full CPU Usage

The code is configured to use **all available CPU cores** by default. Here's how to verify and troubleshoot.

###  Quick Check

The code prints the number of workers it will use:

```bash
python chi2_homogeneity.py input.tsv output.tsv
```

Output shows:
```
Available CPU cores: 72
Processing probesets in parallel using 72 worker(s)...
Note: Using 72 parallel workers for processing
```

### Test Parallelism

Run the parallelism test to verify your system can use all cores:

```bash
python test_parallelism.py 72 200
```

This will:
- Spawn 72 workers
- Run 200 tasks
- Report which unique worker PIDs were used
- Show if all cores are being utilized

## Why You Might See < 100% CPU Usage

### 1. Worker Spawn Time (Normal)

With `spawn` multiprocessing context, workers take time to start:
- **First few seconds**: Only 10-20 cores active (workers starting)
- **After 5-10 seconds**: All 72 cores active
- **Near end**: Cores finish at different times

**Example throughput progression:**
```
  5s:  1.17 probeset/s  (workers starting)
 10s: 13.74 probeset/s  (more workers online)
 15s: 67.95 probeset/s  (most workers active)
 20s: 234.18 probeset/s (all workers active)
```

**This is normal behavior!**

### 2. Task Duration Variability

Some probesets take longer than others:
- **Standard chi-squared**: Very fast (< 1ms)
- **Monte Carlo (small)**: Fast (~10-50ms)
- **Monte Carlo (large)**: Slower (100-1000ms)

If you have a mix, some cores will finish early and sit idle.

### 3. Dataset Size

If you have fewer probesets than workers:
- **50 probesets, 72 workers**: Only 50 cores used
- **1000 probesets, 72 workers**: All 72 cores used

### 4. I/O Bound Operations

At the start, data loading can be I/O bound:
- Reading large TSV files
- Pickl

ing data to send to workers

Once past initial load, processing becomes CPU bound.

## How to Monitor CPU Usage

### Real-time Monitoring

**During execution:**
```bash
# Terminal 1: Run analysis
python chi2_homogeneity.py data.tsv results.tsv --show-progress

# Terminal 2: Monitor CPU
htop  # or 'top'
```

Look for:
- Multiple `python` processes (one per worker)
- CPU% increasing to ~7200% total (72 cores × 100%)

### Average CPU Usage

The total execution time tells you average utilization:

```bash
time python chi2_homogeneity.py data.tsv results.tsv
```

- **72 cores at 100%**: Task completes ~72x faster than single core
- **10 cores at 100%**: Task completes ~10x faster than single core

## Common Issues and Solutions

### Issue: Only 10 Cores Active

**Symptom**: htop shows only ~10 Python processes

**Possible Causes:**

1. **Environment variable limiting workers**
   ```bash
   # Check for limits
   echo $OMP_NUM_THREADS
   echo $MKL_NUM_THREADS

   # If set to 10, unset them
   unset OMP_NUM_THREADS
   unset MKL_NUM_THREADS
   export OMP_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   ```

2. **Explicit --n-workers flag**
   ```bash
   # Check your command
   python chi2_homogeneity.py data.tsv results.tsv --n-workers 10  # ← limiting to 10!

   # Fix: Remove --n-workers or set to 72
   python chi2_homogeneity.py data.tsv results.tsv --n-workers 72
   # Or omit it to use all cores:
   python chi2_homogeneity.py data.tsv results.tsv
   ```

3. **Small dataset**
   - If you have < 20 probesets, you won't see 72 cores active
   - Solution: This is expected, not an issue

### Issue: Slow Start, Then Fast

**Symptom**: First minute shows low CPU, then ramps up

**Cause**: Normal spawn behavior - workers take time to start

**Solution**: This is expected! Wait 30-60 seconds for all workers to spawn.

**Workaround**: Use fewer workers for faster startup:
```bash
python chi2_homogeneity.py data.tsv results.tsv --n-workers 20
```

### Issue: Inconsistent CPU Usage

**Symptom**: CPU usage fluctuates between 1000% and 7000%

**Causes:**
1. **Mixed workload**: Some probesets use Standard test (fast), others use Monte Carlo (slow)
2. **Task completion timing**: Workers finish at different times
3. **I/O pauses**: Occasional disk/memory access

**Solution**: This is normal for heterogeneous workloads. Average throughput matters more than instantaneous CPU%.

## Optimization Tips

### 1. For Maximum Throughput

Use all cores and disable progress for zero overhead:

```bash
python chi2_homogeneity.py data.tsv results.tsv
# No --show-progress flag
```

### 2. For Faster Startup

Use fewer workers to reduce spawn time:

```bash
python chi2_homogeneity.py data.tsv results.tsv --n-workers 20
```

### 3. For Consistent Workload

Reduce Monte Carlo simulations if not needed:

```bash
python chi2_homogeneity.py data.tsv results.tsv --n-simulations 2000
```

(Default is 10,000 simulations)

### 4. Check NumPy/SciPy Threading

NumPy and SciPy can use multiple threads per worker, which can over-subscribe cores:

```bash
# Force single-threaded NumPy/SciPy
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

python chi2_homogeneity.py data.tsv results.tsv
```

This prevents 72 workers × 4 threads each = 288 threads competing.

## Benchmarking

### Measure Actual Speedup

```bash
# Single core
time python chi2_homogeneity.py data.tsv results1.tsv --n-workers 1

# All cores
time python chi2_homogeneity.py data.tsv results2.tsv

# Compare times
# Speedup = single_core_time / all_cores_time
# Expected: 30-60x speedup (not 72x due to overhead)
```

### Use the Benchmark Script

```bash
python benchmark_performance.py data.tsv
```

This measures:
- Legacy vs optimized performance
- Single worker vs multi-worker
- Actual speedup achieved

## Expected Performance

### Ideal Scaling (Pure CPU Bound)

- **1 core**: 1x baseline
- **10 cores**: 9-10x faster
- **72 cores**: 40-60x faster (not 72x due to overhead)

### Real-World Scaling

With real data:
- **Initial ramp-up**: 5-10 seconds at 10-20% utilization
- **Steady state**: 70-90% average utilization
- **Wind-down**: 5-10 seconds as workers finish

**Example:**
- 1000 probesets
- 72 cores
- Total time: 60 seconds
- Average throughput: 16.7 probeset/s
- Effective parallelism: ~50x (not 72x, but excellent)

## Still Having Issues?

### Run Diagnostics

```bash
# 1. Check environment
python test_parallelism.py

# 2. Check process count during execution
python chi2_homogeneity.py data.tsv results.tsv --show-progress &
PID=$!
while kill -0 $PID 2>/dev/null; do
    echo "Active Python processes: $(pgrep -f 'python.*chi2_homogeneity' | wc -l)"
    sleep 2
done

# 3. Profile CPU usage
python -m cProfile -o profile.stats chi2_homogeneity.py data.tsv results.tsv --n-workers 1
# Analyze with: python -m pstats profile.stats
```

### Report the Issue

If you confirmed < 72 cores are used despite:
- ✓ test_parallelism.py shows 72 cores work
- ✓ Environment variables are not limiting
- ✓ --n-workers not set or set to 72
- ✓ Large dataset (> 200 probesets)

Then please report with:
```bash
python test_parallelism.py 72 200 > parallelism_test.txt 2>&1
python chi2_homogeneity.py data.tsv results.tsv --show-progress 2>&1 | head -50 > chi2_test.txt
# Share both output files
```

## Summary

**The code DOES use all 72 cores** by default. You should see:

✓ "Using 72 parallel workers" in output
✓ 72-73 Python processes in htop (1 main + 72 workers)
✓ CPU% ranging from 1000-7200% during execution
✓ Throughput increasing as workers spawn
✓ 30-60x speedup vs single core

If you're seeing different behavior, use the diagnostics above to identify the specific cause.
