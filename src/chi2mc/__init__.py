"""
Chi-squared homogeneity test for genotype frequencies across batches.

Two-tier approach:
1. Standard chi-squared test for tables with all expected frequencies >= 5
2. Monte Carlo chi-squared test for tables with any expected frequency < 5
"""

from chi2mc.chi2_homogeneity import (
    NUMBA_AVAILABLE,
    TQDM_AVAILABLE,
    check_minimum_expected_frequency,
    chi2_monte_carlo_test,
    chi2_standard_test,
    create_contingency_table,
    load_data,
    run_two_tier_chi2_test,
)

__all__ = [
    "load_data",
    "create_contingency_table",
    "check_minimum_expected_frequency",
    "chi2_standard_test",
    "chi2_monte_carlo_test",
    "run_two_tier_chi2_test",
    "NUMBA_AVAILABLE",
    "TQDM_AVAILABLE",
]
