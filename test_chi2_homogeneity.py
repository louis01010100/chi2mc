"""
Unit tests for chi-squared homogeneity test module.
Test-driven development: Tests written before implementation.
"""
import pytest
import polars as pl
import numpy as np
from chi2_homogeneity import (
    load_data,
    create_contingency_table,
    check_minimum_expected_frequency,
    chi2_standard_test,
    chi2_monte_carlo_test,
    run_two_tier_chi2_test,
)


@pytest.fixture
def sample_data():
    """Sample data for testing - PS001 has small n_NC, PS002 has all cells >= 5."""
    return pl.DataFrame({
        'probeset_id': ['PS001', 'PS001', 'PS001', 'PS002', 'PS002', 'PS002'],
        'batch_name': ['batch1', 'batch2', 'batch3', 'batch1', 'batch2', 'batch3'],
        'n_AA': [45, 48, 42, 120, 115, 118],
        'n_AB': [52, 50, 55, 85, 90, 88],
        'n_BB': [38, 40, 37, 30, 32, 31],
        'n_NC': [2, 1, 3, 15, 18, 17],
    })


@pytest.fixture
def small_sample_data():
    """Small sample data that requires Monte Carlo."""
    return pl.DataFrame({
        'probeset_id': ['PS003', 'PS003', 'PS003'],
        'batch_name': ['batch1', 'batch2', 'batch3'],
        'n_AA': [2, 1, 3],
        'n_AB': [3, 4, 2],
        'n_BB': [1, 2, 1],
        'n_NC': [1, 1, 1],
    })


class TestDataLoading:
    """Test data loading functionality."""

    def test_load_data_returns_dataframe(self, tmp_path):
        """Test that load_data returns a polars DataFrame."""
        # Create temporary test file
        test_file = tmp_path / "test_data.tsv"
        test_file.write_text("probeset_id\tbatch_name\tn_AA\tn_AB\tn_BB\tn_NC\n"
                            "PS001\tbatch1\t45\t52\t38\t2\n")

        df = load_data(str(test_file))
        assert isinstance(df, pl.DataFrame)

    def test_load_data_has_correct_columns(self, tmp_path):
        """Test that loaded data has required columns."""
        test_file = tmp_path / "test_data.tsv"
        test_file.write_text("probeset_id\tbatch_name\tn_AA\tn_AB\tn_BB\tn_NC\n"
                            "PS001\tbatch1\t45\t52\t38\t2\n")

        df = load_data(str(test_file))
        required_cols = ['probeset_id', 'batch_name', 'n_AA', 'n_AB', 'n_BB', 'n_NC']
        assert all(col in df.columns for col in required_cols)


class TestContingencyTable:
    """Test contingency table creation."""

    def test_contingency_table_shape(self, sample_data):
        """Test that contingency table has correct shape."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS001')
        table = create_contingency_table(probeset_data)

        # Should be (n_batches, 4 genotypes)
        assert table.shape == (3, 4)

    def test_contingency_table_values(self, sample_data):
        """Test that contingency table has correct values."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS001')
        table = create_contingency_table(probeset_data)

        # Check first batch values (after sorting by batch_name)
        assert table[0, 0] == 45  # n_AA for batch1
        assert table[0, 1] == 52  # n_AB for batch1
        assert table[0, 2] == 38  # n_BB for batch1
        assert table[0, 3] == 2   # n_NC for batch1

    def test_contingency_table_no_negative_values(self, sample_data):
        """Test that contingency table has no negative values."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS001')
        table = create_contingency_table(probeset_data)

        assert np.all(table >= 0)


class TestMinimumExpectedFrequency:
    """Test minimum expected frequency checking."""

    def test_large_sample_passes_check(self, sample_data):
        """Test that large samples pass the minimum frequency check."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS002')
        table = create_contingency_table(probeset_data)

        passes_check = check_minimum_expected_frequency(table, min_count=5)
        assert passes_check == True

    def test_small_sample_fails_check(self, small_sample_data):
        """Test that small samples fail the minimum frequency check."""
        probeset_data = small_sample_data.filter(pl.col('probeset_id') == 'PS003')
        table = create_contingency_table(probeset_data)

        passes_check = check_minimum_expected_frequency(table, min_count=5)
        assert passes_check == False

    def test_edge_case_exactly_five(self):
        """Test edge case where all cells have exactly 5."""
        table = np.array([[5, 5, 5, 5],
                         [5, 5, 5, 5]])

        passes_check = check_minimum_expected_frequency(table, min_count=5)
        assert passes_check == True


class TestStandardChi2:
    """Test standard chi-squared test."""

    def test_standard_chi2_returns_dict(self, sample_data):
        """Test that standard chi2 returns a dictionary with required keys."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS002')
        table = create_contingency_table(probeset_data)

        result = chi2_standard_test(table)

        assert isinstance(result, dict)
        assert 'statistic' in result
        assert 'pvalue' in result
        assert 'dof' in result
        assert 'method' in result

    def test_standard_chi2_method_label(self, sample_data):
        """Test that method is labeled correctly."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS002')
        table = create_contingency_table(probeset_data)

        result = chi2_standard_test(table)
        assert result['method'] == 'standard'

    def test_standard_chi2_statistic_positive(self, sample_data):
        """Test that chi2 statistic is non-negative."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS002')
        table = create_contingency_table(probeset_data)

        result = chi2_standard_test(table)
        assert result['statistic'] >= 0

    def test_standard_chi2_pvalue_range(self, sample_data):
        """Test that p-value is between 0 and 1."""
        probeset_data = sample_data.filter(pl.col('probeset_id') == 'PS002')
        table = create_contingency_table(probeset_data)

        result = chi2_standard_test(table)
        assert 0 <= result['pvalue'] <= 1


class TestMonteCarloChi2:
    """Test Monte Carlo chi-squared test."""

    def test_monte_carlo_chi2_returns_dict(self, small_sample_data):
        """Test that Monte Carlo chi2 returns a dictionary with required keys."""
        probeset_data = small_sample_data.filter(pl.col('probeset_id') == 'PS003')
        table = create_contingency_table(probeset_data)

        result = chi2_monte_carlo_test(table, n_simulations=1000)

        assert isinstance(result, dict)
        assert 'statistic' in result
        assert 'pvalue' in result
        assert 'dof' in result
        assert 'method' in result
        assert 'n_simulations' in result

    def test_monte_carlo_method_label(self, small_sample_data):
        """Test that method is labeled correctly."""
        probeset_data = small_sample_data.filter(pl.col('probeset_id') == 'PS003')
        table = create_contingency_table(probeset_data)

        result = chi2_monte_carlo_test(table, n_simulations=1000)
        # Method can be 'monte_carlo' or 'monte_carlo_numba' depending on numba availability
        assert result['method'] in ['monte_carlo', 'monte_carlo_numba']

    def test_monte_carlo_pvalue_range(self, small_sample_data):
        """Test that p-value is between 0 and 1."""
        probeset_data = small_sample_data.filter(pl.col('probeset_id') == 'PS003')
        table = create_contingency_table(probeset_data)

        result = chi2_monte_carlo_test(table, n_simulations=1000)
        assert 0 <= result['pvalue'] <= 1

    def test_monte_carlo_reproducibility(self, small_sample_data):
        """Test that Monte Carlo results are reproducible with same seed."""
        probeset_data = small_sample_data.filter(pl.col('probeset_id') == 'PS003')
        table = create_contingency_table(probeset_data)

        result1 = chi2_monte_carlo_test(table, n_simulations=1000, random_seed=42)
        result2 = chi2_monte_carlo_test(table, n_simulations=1000, random_seed=42)

        assert result1['pvalue'] == result2['pvalue']
        assert result1['statistic'] == result2['statistic']


class TestTwoTierChi2:
    """Test the two-tier chi-squared test."""

    def test_two_tier_returns_dataframe(self, tmp_path):
        """Test that two-tier test returns a DataFrame."""
        test_file = tmp_path / "test_data.tsv"
        test_file.write_text("probeset_id\tbatch_name\tn_AA\tn_AB\tn_BB\tn_NC\n"
                            "PS001\tbatch1\t45\t52\t38\t2\n"
                            "PS001\tbatch2\t48\t50\t40\t1\n")

        results = run_two_tier_chi2_test(str(test_file), n_workers=1)
        assert isinstance(results, pl.DataFrame)

    def test_two_tier_output_columns(self, tmp_path):
        """Test that output has correct columns."""
        test_file = tmp_path / "test_data.tsv"
        test_file.write_text("probeset_id\tbatch_name\tn_AA\tn_AB\tn_BB\tn_NC\n"
                            "PS001\tbatch1\t45\t52\t38\t2\n"
                            "PS001\tbatch2\t48\t50\t40\t1\n")

        results = run_two_tier_chi2_test(str(test_file), n_workers=1)

        required_cols = ['probeset_id', 'chi2_statistic', 'pvalue', 'dof', 'method']
        assert all(col in results.columns for col in required_cols)

    def test_two_tier_uses_standard_for_large_samples(self):
        """Test that standard method is used for large samples."""
        # Create data with large samples
        data = pl.DataFrame({
            'probeset_id': ['PS001', 'PS001'],
            'batch_name': ['batch1', 'batch2'],
            'n_AA': [50, 48],
            'n_AB': [50, 52],
            'n_BB': [50, 45],
            'n_NC': [50, 55],
        })

        # Save to temp file
        import tempfile
        import os
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            data.write_csv(f.name, separator='\t')
            temp_file = f.name

        results = run_two_tier_chi2_test(temp_file, n_workers=1)

        # Clean up
        os.unlink(temp_file)

        assert results.filter(pl.col('probeset_id') == 'PS001')['method'][0] == 'standard'

    def test_two_tier_uses_monte_carlo_for_small_samples(self):
        """Test that Monte Carlo method is used for small samples."""
        # Create data with small samples
        data = pl.DataFrame({
            'probeset_id': ['PS003', 'PS003'],
            'batch_name': ['batch1', 'batch2'],
            'n_AA': [2, 1],
            'n_AB': [3, 4],
            'n_BB': [1, 2],
            'n_NC': [1, 1],
        })

        # Save to temp file
        import tempfile
        import os
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            data.write_csv(f.name, separator='\t')
            temp_file = f.name

        results = run_two_tier_chi2_test(temp_file, n_simulations=1000, n_workers=1)

        # Clean up
        os.unlink(temp_file)

        # Method can be 'monte_carlo' or 'monte_carlo_numba' depending on numba availability
        method = results.filter(pl.col('probeset_id') == 'PS003')['method'][0]
        assert method in ['monte_carlo', 'monte_carlo_numba']

    def test_two_tier_handles_multiple_probesets(self, tmp_path):
        """Test that two-tier test handles multiple probesets correctly."""
        test_file = tmp_path / "test_data.tsv"
        test_file.write_text("probeset_id\tbatch_name\tn_AA\tn_AB\tn_BB\tn_NC\n"
                            "PS001\tbatch1\t45\t52\t38\t2\n"
                            "PS001\tbatch2\t48\t50\t40\t1\n"
                            "PS002\tbatch1\t120\t85\t30\t5\n"
                            "PS002\tbatch2\t115\t90\t32\t3\n")

        results = run_two_tier_chi2_test(str(test_file), n_workers=1)

        assert len(results) == 2
        assert 'PS001' in results['probeset_id'].to_list()
        assert 'PS002' in results['probeset_id'].to_list()

    def test_parallel_processing_gives_same_results(self, tmp_path):
        """Test that parallel processing gives same results as sequential."""
        test_file = tmp_path / "test_data.tsv"
        test_file.write_text("probeset_id\tbatch_name\tn_AA\tn_AB\tn_BB\tn_NC\n"
                            "PS001\tbatch1\t45\t52\t38\t2\n"
                            "PS001\tbatch2\t48\t50\t40\t1\n"
                            "PS002\tbatch1\t120\t85\t30\t15\n"
                            "PS002\tbatch2\t115\t90\t32\t18\n"
                            "PS003\tbatch1\t2\t3\t1\t1\n"
                            "PS003\tbatch2\t1\t4\t2\t1\n")

        # Sequential processing
        results_seq = run_two_tier_chi2_test(
            str(test_file),
            n_simulations=100,  # Reduced for faster testing
            random_seed=42,
            n_workers=1
        )

        # Parallel processing with 2 workers
        results_par = run_two_tier_chi2_test(
            str(test_file),
            n_simulations=100,  # Reduced for faster testing
            random_seed=42,
            n_workers=2
        )

        # Both should have same probesets
        assert set(results_seq['probeset_id'].to_list()) == set(results_par['probeset_id'].to_list())

        # Check that results match for each probeset
        for probeset_id in results_seq['probeset_id'].to_list():
            seq_row = results_seq.filter(pl.col('probeset_id') == probeset_id)
            par_row = results_par.filter(pl.col('probeset_id') == probeset_id)

            assert seq_row['method'][0] == par_row['method'][0]
            assert seq_row['chi2_statistic'][0] == par_row['chi2_statistic'][0]
            assert seq_row['pvalue'][0] == par_row['pvalue'][0]
