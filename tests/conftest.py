"""Pytest configuration for golden master tests."""

import pytest


def pytest_addoption(parser):
    """Add --run-slow option to pytest."""
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run slow integration tests",
    )


def pytest_configure(config):
    """Register the slow marker."""
    config.addinivalue_line(
        "markers",
        "slow: mark test as slow (run with --run-slow)",
    )


def pytest_collection_modifyitems(config, items):
    """Skip slow tests unless --run-slow is specified."""
    if config.getoption("--run-slow"):
        return

    skip_slow = pytest.mark.skip(reason="Slow test, use --run-slow to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
