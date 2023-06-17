from snakemake.deployment.version_compare import compare_version_geq
import pytest


def test_versiongt():
    # Test cases are based upon releases found in apptainer and singularity
    assert compare_version_geq("1.0.0", "1.0.0") == True
    assert compare_version_geq("1.0.1", "1.0.0") == True
    assert compare_version_geq("1.1", "1.0.0") == True
    assert compare_version_geq("1.0.0", "0.9.9") == True
    assert compare_version_geq("0.9.9", "1.0") == False
    assert compare_version_geq("0.1.0", "1.0") == False
    assert compare_version_geq("1.0.0-rc.1", "1.0.0") == False
    assert compare_version_geq("1.2.0-rc.1", "1.0.0") == True
    assert compare_version_geq("2.4.1", "2.4.1") == True
    assert compare_version_geq("2.4.5", "2.4.1") == True
    assert compare_version_geq("2.3.1", "2.4.1") == False
    assert compare_version_geq("2.4.1-rc1", "2.4.1") == False
    assert compare_version_geq("3.10.0-rc.1", "3.10.0") == False
    assert compare_version_geq("2023.03-0", "4.2") == True
    assert compare_version_geq("2023.03-1", "4.2") == True
    assert compare_version_geq("4.2", "2023.03-1") == False
    assert compare_version_geq("2023.03-0", "2023.03-1") == False
