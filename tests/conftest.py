import os
import sys
import pytest
from pathlib import Path

from snakemake.common import ON_WINDOWS
from snakemake.utils import _find_bash_on_windows

skip_on_windows = pytest.mark.skipif(ON_WINDOWS, reason="Unix stuff")


@pytest.fixture(autouse=True)
def reset_paths_between_tests():
    """Ensure that changes to sys.path are reset between tests"""
    org_path = sys.path.copy()
    yield
    sys.path = org_path


bash_cmd = _find_bash_on_windows()

if ON_WINDOWS and bash_cmd:

    @pytest.fixture(autouse=True)
    def prepend_usable_bash_to_path(monkeypatch):
        monkeypatch.setenv("PATH", os.path.dirname(bash_cmd), prepend=os.pathsep)
