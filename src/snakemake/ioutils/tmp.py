import os
from pathlib import Path
from typing import List


def get_tmp(tmpdirs: List) -> str:
    """
    Given a list of possible temp folder paths, return the first path that
    can be created and where the current user has read/write/exe permissions.
    If none, return `system_tmpdir`.
    """

    for tmpdir in tmpdirs:
        try:
            Path(tmpdir).mkdir(parents=True, exist_ok=True)
        except PermissionError:
            continue
        else:
            if (
                os.access(tmpdir, os.R_OK)
                and os.access(tmpdir, os.W_OK)
                and os.access(tmpdir, os.X_OK)
            ):
                return tmpdir
            else:
                continue
    return "system_tmpdir"
