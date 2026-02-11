"""Shared fixtures for chi2mc tests."""

import pytest
import polars as pl


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
