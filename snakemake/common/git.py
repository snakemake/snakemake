import os
import re

from snakemake.exceptions import WorkflowError


def split_git_path(path):
    file_sub = re.sub(r"^git\+file:/+", "/", path)
    (file_path, version) = file_sub.split("@")
    file_path = os.path.realpath(file_path)
    root_path = get_git_root(file_path)
    if file_path.startswith(root_path):
        file_path = file_path[len(root_path) :].lstrip("/")
    return (root_path, file_path, version)


def get_git_root(path):
    """
    Args:
        path: (str) Path a to a directory/file that is located inside the repo
    Returns:
        path to the root folder for git repo
    """
    import git

    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        return git_repo.git.rev_parse("--show-toplevel")
    except git.exc.NoSuchPathError:
        tail, _ = os.path.split(path)
        return get_git_root_parent_directory(tail, path)


def get_git_root_parent_directory(path, input_path):
    """
    This function will recursively go through parent directories until a git
    repository is found or until no parent directories are left, in which case
    an error will be raised. This is needed when providing a path to a
    file/folder that is located on a branch/tag not currently checked out.

    Args:
        path: (str) Path a to a directory that is located inside the repo
        input_path: (str) origin path, used when raising WorkflowError
    Returns:
        path to the root folder for git repo
    """
    import git

    try:
        git_repo = git.Repo(path, search_parent_directories=True)
        return git_repo.git.rev_parse("--show-toplevel")
    except git.exc.NoSuchPathError:
        tail, _ = os.path.split(path)
        if tail is None:
            raise WorkflowError(
                f"Neither provided git path ({input_path}) "
                + "or parent directories contain a valid git repo."
            )
        else:
            return get_git_root_parent_directory(tail, input_path)


def git_content(git_file):
    """
    This function will extract a file from a git repository, one located on
    the filesystem.
    The expected format is git+file:///path/to/your/repo/path_to_file@version

    Args:
      env_file (str): consist of path to repo, @, version, and file information
                      Ex: git+file:///home/smeds/snakemake-wrappers/bio/fastqc/wrapper.py@0.19.3
    Returns:
        file content or None if the expected format isn't meet
    """
    import git

    if git_file.startswith("git+file:"):
        (root_path, file_path, version) = split_git_path(git_file)
        return git.Repo(root_path).git.show(f"{version}:{file_path}")
    else:
        raise WorkflowError(
            "Provided git path ({}) doesn't meet the "
            "expected format:".format(git_file) + ", expected format is "
            "git+file://PATH_TO_REPO/PATH_TO_FILE_INSIDE_REPO@VERSION"
        )
