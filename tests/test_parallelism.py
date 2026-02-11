"""
Test script to verify parallel worker usage.

This script helps diagnose parallelism issues by:
1. Checking environment variables that might limit threading
2. Monitoring actual CPU usage
3. Verifying the number of workers spawned
"""

import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed


def worker_task(task_id):
    """
    Simulate a CPU-intensive task.
    """
    # Get worker process ID
    pid = os.getpid()

    # Simulate work
    result = 0
    for i in range(10**7):
        result += i % 7

    return {
        'task_id': task_id,
        'pid': pid,
        'result': result
    }


def test_parallelism(n_tasks=100, n_workers=None):
    """
    Test parallel execution and report statistics.
    """
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()

    print("=" * 80)
    print("PARALLELISM TEST")
    print("=" * 80)
    print(f"System CPU count: {multiprocessing.cpu_count()}")
    print(f"Requested workers: {n_workers}")
    print(f"Total tasks: {n_tasks}")
    print()

    # Check for environment variables that might limit parallelism
    print("Environment variables that affect parallelism:")
    print("-" * 80)
    env_vars = [
        'OMP_NUM_THREADS',
        'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS',
        'NUMBA_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS',
        'MAX_CONCURRENCY'
    ]
    for var in env_vars:
        value = os.environ.get(var, 'Not set')
        print(f"  {var}: {value}")
    print()

    # Run parallel tasks
    print("Running parallel tasks...")
    print("-" * 80)

    start_time = time.time()
    unique_pids = set()
    results = []

    ctx = multiprocessing.get_context('spawn')
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as executor:
        # Submit all tasks
        futures = {
            executor.submit(worker_task, i): i
            for i in range(n_tasks)
        }

        # Collect results
        n_completed = 0
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            unique_pids.add(result['pid'])
            n_completed += 1

            # Print progress every 10 tasks
            if n_completed % 10 == 0:
                percent = 100 * n_completed / n_tasks
                print(f"Progress: {n_completed}/{n_tasks} ({percent:.1f}%) - "
                      f"Unique workers used so far: {len(unique_pids)}", end='\r')

    elapsed = time.time() - start_time

    print()  # New line after progress
    print()

    # Report statistics
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)
    print(f"Total time: {elapsed:.2f} seconds")
    print(f"Tasks per second: {n_tasks / elapsed:.2f}")
    print(f"Unique worker PIDs: {len(unique_pids)}")
    print(f"Expected workers: {n_workers}")

    if len(unique_pids) < n_workers:
        print()
        print("WARNING: Fewer workers used than expected!")
        print(f"   Expected: {n_workers} workers")
        print(f"   Actual:   {len(unique_pids)} workers")
        print()
        print("Possible causes:")
        print("1. Environment variables limiting parallelism (see above)")
        print("2. System resource constraints")
        print("3. ProcessPoolExecutor not spawning all workers")
        print()
        print("To investigate:")
        print("- Check 'top' or 'htop' during execution")
        print("- Try: export OMP_NUM_THREADS=1")
        print("- Try: export MKL_NUM_THREADS=1")
    else:
        print()
        print("All expected workers were used!")

    print("=" * 80)

    return len(unique_pids) == n_workers


def main():
    if len(sys.argv) > 1:
        n_workers = int(sys.argv[1])
    else:
        n_workers = None  # Use all CPUs

    if len(sys.argv) > 2:
        n_tasks = int(sys.argv[2])
    else:
        n_tasks = 100

    success = test_parallelism(n_tasks=n_tasks, n_workers=n_workers)

    if not success:
        print("\nParallelism test FAILED - not all workers were used")
        sys.exit(1)
    else:
        print("\nParallelism test PASSED")
        sys.exit(0)


if __name__ == '__main__':
    main()
