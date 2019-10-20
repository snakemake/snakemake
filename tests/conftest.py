import sys
import pytest

from snakemake.utils import ON_WINDOWS

skip_on_windows = pytest.mark.skipif(ON_WINDOWS, reason="Unix stuff")


@pytest.fixture(autouse=True)
def reset_paths_between_tests():
    """Ensure that changes to sys.path are reset between tests"""
    org_path = sys.path.copy()
    yield
    sys.path = org_path
