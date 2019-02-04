import pytest

from snakemake.utils import ON_WINDOWS

skip_on_windows = pytest.mark.skipif(ON_WINDOWS, reason="Unix stuff")
