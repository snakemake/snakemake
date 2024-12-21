__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import posixpath
import re
import os
import shutil
import stat
from snakemake import utils
import tempfile
import io
from abc import ABC, abstractmethod
from urllib.parse import unquote

from snakemake_interface_executor_plugins.settings import ExecMode
from snakemake.common import (
    ON_WINDOWS,
    is_local_file,
    get_appdirs,
    parse_uri,
    smart_join,
)
from snakemake.exceptions import WorkflowError, SourceFileError
from snakemake.common.git import split_git_path
from snakemake.logging import logger


def _check_git_args(tag: str = None, branch: str = None, commit: str = None):
    n_refs = sum(1 for ref in (tag, branch, commit) if ref is not None)
    if n_refs != 1:
        raise SourceFileError(
            "exactly one of tag, branch, or commit must be specified."
        )


class SourceFile(ABC):
    @abstractmethod
    def get_path_or_uri(self): ...

    @abstractmethod
    def is_persistently_cacheable(self): ...

    def get_cache_path(self):
        uri = parse_uri(self.get_path_or_uri())
        return os.path.join(uri.scheme, unquote(uri.uri_path.lstrip("/")))

    def get_basedir(self):
        path = os.path.dirname(self.get_path_or_uri())
        return self.__class__(path)

    @abstractmethod
    def get_filename(self): ...

    def join(self, path):
        if isinstance(path, SourceFile):
            path = path.get_path_or_uri()
        return self.__class__(smart_join(self.get_path_or_uri(), path))

    def mtime(self):
        """If possible, return mtime of the file. Otherwise, return None."""
        return None

    @property
    @abstractmethod
    def is_local(self): ...

    def __hash__(self):
        return self.get_path_or_uri().__hash__()

    def __eq__(self, other):
        if isinstance(other, SourceFile):
            return self.get_path_or_uri() == other.get_path_or_uri()
        return False

    def __str__(self):
        return self.get_path_or_uri()

    def simplify_path(self):
        return self


class GenericSourceFile(SourceFile):
    def __init__(self, path_or_uri):
        self.path_or_uri = path_or_uri

    def get_path_or_uri(self):
        return self.path_or_uri

    def get_filename(self):
        return os.path.basename(self.path_or_uri)

    def is_persistently_cacheable(self):
        return False

    @property
    def is_local(self):
        return False


class LocalSourceFile(SourceFile):
    def __init__(self, path):
        self.path = path

    def get_path_or_uri(self):
        return self.path

    def is_persistently_cacheable(self):
        return False

    def get_filename(self):
        return os.path.basename(self.path)

    def abspath(self):
        return LocalSourceFile(os.path.abspath(self.path))

    def isabs(self):
        return os.path.isabs(self.path)

    def simplify_path(self):
        return utils.simplify_path(self.path)

    def mtime(self):
        return os.stat(self.path).st_mtime

    def __fspath__(self):
        return self.path

    @property
    def is_local(self):
        return True


class LocalGitFile(SourceFile):
    def __init__(
        self, repo_path, path: str, tag: str = None, ref: str = None, commit: str = None
    ):
        _check_git_args(tag, ref, commit)
        self.tag = tag
        self.commit = commit
        self._ref = ref
        self.repo_path = repo_path
        self.path = path

    def get_path_or_uri(self):
        return "git+file://{}/{}@{}".format(
            os.path.abspath(self.repo_path), self.path, self.ref
        )

    def join(self, path):
        path = os.path.normpath("/".join((self.path, path)))
        if ON_WINDOWS:
            # convert back to URL separators
            # (win specific separators are introduced by normpath above)
            path = path.replace("\\", "/")
        return LocalGitFile(
            self.repo_path, path, tag=self.tag, ref=self._ref, commit=self.commit
        )

    def get_basedir(self):
        return self.__class__(
            repo_path=self.repo_path,
            path=os.path.dirname(self.path),
            tag=self.tag,
            commit=self.commit,
            ref=self._ref,
        )

    def is_persistently_cacheable(self):
        return False

    def get_filename(self):
        return posixpath.basename(self.path)

    @property
    def ref(self):
        return self.tag or self.commit or self._ref

    @property
    def is_local(self):
        return True


class HostingProviderFile(SourceFile):
    """Marker for denoting github source files from releases."""

    valid_repo = re.compile("^.+/.+$")

    def __init__(
        self,
        repo: str = None,
        path: str = None,
        tag: str = None,
        branch: str = None,
        commit: str = None,
    ):
        if repo is None:
            raise SourceFileError("repo must be given")
        if not self.__class__.valid_repo.match(repo):
            raise SourceFileError(
                "repo {} is not a valid repo specification (must be given as owner/name)."
            )

        _check_git_args(tag, branch, commit)

        if path is None:
            raise SourceFileError("path must be given")

        if not all(
            isinstance(item, str)
            for item in (repo, path, tag, branch, commit)
            if item is not None
        ):
            raise SourceFileError("arguments must be given as str.")

        self.repo = repo
        self.tag = tag
        self.commit = commit
        self.branch = branch
        self.path = path.strip("/")
        self.token = ""

    def is_persistently_cacheable(self):
        return bool(self.tag or self.commit)

    def get_filename(self):
        return os.path.basename(self.path)

    @property
    def ref(self):
        return self.tag or self.commit or self.branch

    def get_basedir(self):
        return self.__class__(
            repo=self.repo,
            path=os.path.dirname(self.path),
            tag=self.tag,
            commit=self.commit,
            branch=self.branch,
        )

    def join(self, path):
        path = os.path.normpath(f"{self.path}/{path}")
        if ON_WINDOWS:
            # convert back to URL separators
            # (win specific separators are introduced by normpath above)
            path = path.replace("\\", "/")
        return self.__class__(
            repo=self.repo,
            path=path,
            tag=self.tag,
            commit=self.commit,
            branch=self.branch,
        )

    @property
    def is_local(self):
        return False


class GithubFile(HostingProviderFile):
    def __init__(
        self,
        repo: str,
        path: str,
        tag: str = None,
        branch: str = None,
        commit: str = None,
    ):
        super().__init__(repo, path, tag, branch, commit)
        self.token = os.environ.get("GITHUB_TOKEN", "")

    def get_path_or_uri(self):
        auth = f":{self.token}@" if self.token else ""
        return "https://{}raw.githubusercontent.com/{}/{}/{}".format(
            auth, self.repo, self.ref, self.path
        )


class GitlabFile(HostingProviderFile):
    def __init__(
        self,
        repo: str,
        path: str,
        tag: str = None,
        branch: str = None,
        commit: str = None,
        host: str = None,
    ):
        super().__init__(repo, path, tag, branch, commit)
        self.host = host
        self.token = os.environ.get("GITLAB_TOKEN", "")

    def get_path_or_uri(self):
        from urllib.parse import quote

        auth = f"&private_token={self.token}" if self.token else ""
        return "https://{}/api/v4/projects/{}/repository/files/{}/raw?ref={}{}".format(
            self.host or "gitlab.com",
            quote(self.repo, safe=""),
            quote(self.path, safe=""),
            self.ref,
            auth,
        )


def infer_source_file(path_or_uri, basedir: SourceFile = None):
    if isinstance(path_or_uri, SourceFile):
        if basedir is None or isinstance(path_or_uri, HostingProviderFile):
            return path_or_uri
        else:
            path_or_uri = path_or_uri.get_path_or_uri()
    if isinstance(path_or_uri, Path):
        path_or_uri = str(path_or_uri)
    if not isinstance(path_or_uri, str):
        raise SourceFileError(
            "must be given as Python string or one of the predefined source file marker types (see docs)"
        )
    if is_local_file(path_or_uri):
        # either local file or relative to some remote basedir
        for schema in ("file://", "file:"):
            if path_or_uri.startswith(schema):
                path_or_uri = path_or_uri[len(schema) :]
                break
        if not os.path.isabs(path_or_uri) and basedir is not None:
            return basedir.join(path_or_uri)
        return LocalSourceFile(path_or_uri)
    if path_or_uri.startswith("git+file:"):
        try:
            root_path, file_path, ref = split_git_path(path_or_uri)
        except Exception as e:
            raise WorkflowError(
                f"Failed to read source {path_or_uri} from git repo.", e
            )
        return LocalGitFile(root_path, file_path, ref=ref)
    # something else
    return GenericSourceFile(path_or_uri)


class SourceCache:
    cache_whitelist = [
        r"https://raw.githubusercontent.com/snakemake/snakemake-wrappers/\d+\.\d+.\d+"
    ]  # TODO add more prefixes for uris that are save to be cached

    def __init__(self, cache_path: Path, runtime_cache_path: Path = None):
        self.cache_path = cache_path
        os.makedirs(self.cache_path, exist_ok=True)
        if runtime_cache_path is None:
            runtime_cache_parent = self.cache_path / "runtime-cache"
            os.makedirs(runtime_cache_parent, exist_ok=True)
            self.runtime_cache = tempfile.TemporaryDirectory(
                dir=runtime_cache_parent, ignore_cleanup_errors=True
            )
            self._runtime_cache_path = None
        else:
            self._runtime_cache_path = runtime_cache_path
            self.runtime_cache = None
        self.cacheable_prefixes = re.compile("|".join(self.cache_whitelist))

    @property
    def runtime_cache_path(self):
        return self._runtime_cache_path or self.runtime_cache.name

    def open(self, source_file, mode="r"):
        cache_entry = self._cache(source_file)
        return self._open_local_or_remote(
            LocalSourceFile(cache_entry), mode, encoding="utf-8"
        )

    def exists(self, source_file):
        try:
            self._cache(source_file, retries=1)
        except Exception:
            return False
        return True

    def get_path(self, source_file):
        cache_entry = self._cache(source_file)
        return str(cache_entry)

    def _cache_entry(self, source_file: SourceFile) -> Path:
        file_cache_path = source_file.get_cache_path()
        assert file_cache_path

        # TODO add git support to smart_open!
        if source_file.is_persistently_cacheable():
            # check cache
            return self.cache_path / file_cache_path
        else:
            # check runtime cache
            return Path(self.runtime_cache_path) / file_cache_path

    def _cache(self, source_file: SourceFile, retries: int = 3):
        cache_entry = self._cache_entry(source_file)
        if not cache_entry.exists():
            self._do_cache(source_file, cache_entry, retries=retries)
        return cache_entry

    def _do_cache(self, source_file, cache_entry: Path, retries: int = 3):
        # open from origin
        with self._open_local_or_remote(source_file, "rb", retries=retries) as source:
            cache_entry.parent.mkdir(parents=True, exist_ok=True)
            tmp_source = tempfile.NamedTemporaryFile(
                prefix=str(cache_entry),
                delete=False,  # no need to delete since we move it below
            )
            tmp_source.write(source.read())
            tmp_source.close()
            # ensure read and write permissions for owner and group
            os.chmod(
                tmp_source.name,
                stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP,
            )
            # Atomic move to right name.
            # This way we avoid the need to lock.
            shutil.move(tmp_source.name, cache_entry)

        mtime = source_file.mtime()
        if mtime is not None:
            # Set to mtime of original file
            # In case we don't have that mtime, it is fine
            # to just keep the time at the time of caching
            # as mtime.
            os.utime(cache_entry, times=(mtime, mtime))

    def _open_local_or_remote(
        self, source_file: SourceFile, mode, encoding=None, retries: int = 3
    ):
        from reretry.api import retry_call

        if source_file.is_local:
            return self._open(source_file, mode, encoding=encoding)
        else:
            return retry_call(
                self._open,
                [source_file, mode, encoding],
                tries=retries,
                delay=3,
                backoff=2,
                logger=logger,
            )

    def _open(self, source_file: SourceFile, mode, encoding=None):
        from smart_open import open

        if isinstance(source_file, LocalGitFile):
            import git

            return io.BytesIO(
                git.Repo(source_file.repo_path)
                .git.show(f"{source_file.ref}:{source_file.path}")
                .encode()
            )

        path_or_uri = source_file.get_path_or_uri()

        try:
            return open(path_or_uri, mode, encoding=None if "b" in mode else encoding)
        except Exception as e:
            raise WorkflowError(f"Failed to open source file {path_or_uri}", e)
