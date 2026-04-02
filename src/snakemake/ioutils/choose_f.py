from typing import List, Union
from pathlib import Path
from snakemake.io import AnnotatedString


def choose_file(
    file_list: List[Union[Path, AnnotatedString, str]],
    read: bool = True,
    write: bool = True,
    execute: bool = None,
    creatable: bool = None,
) -> Path:
    return _choose_f(
        file_list,
        read=read,
        write=write,
        execute=execute,
        creatable=creatable,
        is_dir=False,
    )


def choose_folder(
    folder_list: List[Union[Path, AnnotatedString, str]],
    read: bool = True,
    write: bool = True,
    open: bool = True,
    creatable: bool = None,
) -> Path:
    return _choose_f(
        folder_list,
        read=read,
        write=write,
        execute=open,
        creatable=creatable,
        is_dir=True,
    )


def choose_tmp(
    folder_list: List[Union[Path, AnnotatedString, str]],
    read: bool = True,
    write: bool = True,
    open: bool = True,
    creatable: bool = True,
) -> Union[Path, str]:
    tmpdir = _choose_f(
        folder_list,
        read=read,
        write=write,
        execute=open,
        creatable=creatable,
        is_dir=True,
    )
    return "system_tmpdir" if tmpdir is None else tmpdir


def _choose_f(
    list: List[Union[Path, AnnotatedString, str]],
    read: bool = None,
    write: bool = None,
    execute: bool = None,
    creatable: bool = None,
    is_dir: bool = None,
) -> Path:
    """
    Given a list of files/folder paths, return the first path that meets criteria.
    """

    def _check_mode(path: Path, mode: List[bool] = [None, None, None]) -> bool:
        """
        Checks if path matches provided mode in the order: read, write, execute.
        `True` specifies that specific mode is set, while `False` specifies the opposite; if `None`, that mode is ignored.
        As an example, `mode = [True, False, None]` would match files/folders that are readable, not writable and where execution mode is irrlelevant (they can be either executable or not).
        """
        import os

        modes = zip([os.R_OK, os.W_OK, os.X_OK], mode)
        return all([os.access(path, mode) for mode, val in modes if val is not None])

    def is_creatable(path: Path) -> bool:
        """
        Checks if path is creatable (whether it exists or not).
        """
        if path == Path(path.root):
            return False
        elif (
            path.parent.exists()
            and path.parent.is_dir()
            and _check_mode(path.parent, [True, True, True])
        ):
            return True
        else:
            return is_creatable(path.parent)

    for i in list:
        # If storage, permissions cannot be checked
        if isinstance(i, AnnotatedString):
            if i.is_storage() and is_dir == i.is_directory():
                return True
        else:
            i = Path(i).expanduser()

            # Check permissions
            if _check_mode(i, [read, write, execute]) and is_dir == i.is_dir():
                return i

            # Checks if it is creatable
            if creatable:
                if is_creatable(i):
                    return i

    return None
