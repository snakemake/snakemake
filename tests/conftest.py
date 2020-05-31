import os
import sys
import pytest
from pathlib import Path

from snakemake.common import ON_WINDOWS
from snakemake.utils import find_bash_on_windows

skip_on_windows = pytest.mark.skipif(ON_WINDOWS, reason="Unix stuff")


@pytest.fixture(autouse=True)
def reset_paths_between_tests():
    """Ensure that changes to sys.path are reset between tests"""
    org_path = sys.path.copy()
    yield
    sys.path = org_path


if ON_WINDOWS:

    @pytest.fixture(autouse=True, scope="session")
    def use_good_bash_on_windows(monkeypatch):
        bash_cmd = find_bash_on_windows()
        if bash_cmd:
            bash_dir = Path(bash_cmd).parent
            envpath = os.getenv("PATH")
            monkeypatch.setenv("PATH", bash_dir + ";" + envpath)
