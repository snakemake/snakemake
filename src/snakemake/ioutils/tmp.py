from typing import List
from pathlib import Path


def get_tmp(tmpdirs: List) -> str:
    """
    Given a list of possible temp folder paths, return the first path that
    can be created and where the current user has read/write/exe permissions.
    If none, return `system_tmpdir`.
    """

    def is_dir_creatable(dir: Path) -> bool:
        import os

        if dir == dir.parent:
            return False
        elif (
            Path(dir.parent).exists()
            and Path(dir.parent).is_dir()
            and os.access(dir.parent, os.R_OK)
            and os.access(dir.parent, os.W_OK)
            and os.access(dir.parent, os.X_OK)
        ):
            return True
        else:
            is_dir_creatable(dir.parent)

    for tmpdir in tmpdirs:
        if is_dir_creatable(Path(tmpdir)):
            return tmpdir
        else:
            continue
    return "system_tmpdir"
