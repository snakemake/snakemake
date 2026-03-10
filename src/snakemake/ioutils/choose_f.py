from typing import List, Union
from pathlib import Path


def choose_file(file_list: List[Union[Path, str]]) -> Path:
    return choose_f(file_list, read=True, write=True)


def choose_folder(folder_list: List[Union[Path, str]]) -> Path:
    return choose_f(
        folder_list,
        read=True,
        write=True,
        execute=True,
    )


def choose_tmp(folder_list: List[Union[Path, str]]) -> Union[Path, str]:
    tmpdir = choose_f(
        folder_list,
        read=True,
        write=True,
        execute=True,
        creatable=True,
    )
    return "system_tmpdir" if tmpdir is None else tmpdir


def choose_f(
    list: List[Union[Path, str]],
    read: bool = None,
    write: bool = None,
    execute: bool = None,
    creatable: bool = None,
) -> Path:
    """
    Given a list of files/folder paths, return the first path that meets criteria.
    """

    def _check_f(path: Path, mode: List[bool] = [None, None, None]) -> bool:
        """
        Checks if path matches provided mode in the order: read, write, execute.
        `True` specifies that specific mode is set, while `False` specifies the oposite; if `None`, that mode is ignored.
        As an example, `mode = [True, False, None]` would match files/folders that are readable, not writable and where execution mode is irrlelevant (they can be either executable or not).
        """
        import os

        modes = zip([os.R_OK, os.W_OK, os.X_OK], mode)
        return all([os.access(path, mode) for mode, val in modes if val is not None])

    def is_creatable(path: Path) -> bool:
        """
        Checks if path is creatable.
        """
        if path == Path(path.root):
            return False
        elif (
            path.parent.exists()
            and path.parent.is_dir()
            and _check_f(path.parent, [True, True, True])
        ):
            return True
        else:
            return is_creatable(path.parent)

    for i in list:
        i = Path(i).expanduser()
        if creatable:
            if is_creatable(i):
                return i

        if _check_f(i, [read, write, execute]):
            return i

    return None
