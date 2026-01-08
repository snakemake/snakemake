__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from functools import partial
from pathlib import Path
import posixpath
import re
import os
import shutil
import stat
import threading
import typing
from typing import TYPE_CHECKING, Any, ClassVar, Dict, List, Optional

import tempfile
import io
from abc import ABC, abstractmethod
from urllib.parse import unquote

from snakemake import utils
from snakemake.utils import format
from snakemake.common import (
    ON_WINDOWS,
    LockFreeWritableFile,
    is_local_file,
    parse_uri,
    smart_join,
)
from snakemake.exceptions import WorkflowError, SourceFileError
from snakemake.common.git import split_git_path
from snakemake.logging import logger
from snakemake.io import check


if TYPE_CHECKING:
    import git


def _format_or_none(pattern: Optional[str], **kwargs: Any):
    return format(pattern, **kwargs) if pattern is not None else None


def _check_git_args(
    tag: Optional[str] = None,
    branch: Optional[str] = None,
    commit: Optional[str] = None,
):
    n_refs = sum(1 for ref in (tag, branch, commit) if ref is not None)
    if n_refs != 1:
        raise SourceFileError(
            "exactly one of tag, branch, or commit must be specified."
        )


def _replace_suffix(path: str, suffix: List[str], replacement: str) -> Optional[str]:
    for suff in suffix:
        if path.endswith(suff):
            return path[: -len(suff)] + replacement
    return None


class SourceFile(ABC):
    @abstractmethod
    def check(self) -> None: ...

    @abstractmethod
    def get_path_or_uri(self, secret_free: bool) -> str: ...

    @abstractmethod
    def is_persistently_cacheable(self) -> bool: ...

    def get_cache_path(self):
        uri = parse_uri(self.get_path_or_uri(secret_free=True))
        return os.path.join(uri.scheme, unquote(uri.uri_path.lstrip("/")))

    def get_basedir(self):
        path = os.path.dirname(self.get_path_or_uri(secret_free=False))
        return self.__class__(path)

    @abstractmethod
    def format(self, **kwargs: Any) -> "typing.Self": ...

    @abstractmethod
    def get_filename(self) -> str: ...

    @abstractmethod
    def endswith(self, suffix: str) -> bool: ...

    @abstractmethod
    def replace_suffix(
        self, suffix: List[str], replacement: str
    ) -> Optional["typing.Self"]: ...

    def join(self, path):
        if isinstance(path, SourceFile):
            path = path.get_path_or_uri(secret_free=False)
        return self.__class__(
            smart_join(base=self.get_path_or_uri(secret_free=False), path=path)
        )

    def mtime(self) -> Optional[float]:
        """If possible, return mtime of the file. Otherwise, return None."""
        return None

    @property
    @abstractmethod
    def is_local(self) -> bool: ...

    def __hash__(self):
        return self.get_path_or_uri(secret_free=True).__hash__()

    def __eq__(self, other):
        if isinstance(other, SourceFile):
            return self.get_path_or_uri(secret_free=True) == other.get_path_or_uri(
                secret_free=True
            )
        return False

    def __str__(self):
        return self.get_path_or_uri(secret_free=True)

    def simplify_path(self):
        return self.get_path_or_uri(secret_free=True)


class GenericSourceFile(SourceFile):
    def __init__(self, path_or_uri):
        self.path_or_uri = path_or_uri

    def endswith(self, suffix: str) -> bool:
        return self.path_or_uri.endswith(suffix)

    def replace_suffix(
        self, suffix: List[str], replacement: str
    ) -> Optional["typing.Self"]:
        repl = _replace_suffix(self.path_or_uri, suffix, replacement)
        if repl is None:
            return None
        else:
            return self.__class__(repl)

    def check(self) -> None:
        pass

    def format(self, **kwargs: Any) -> "typing.Self":
        return self.__class__(format(self.path_or_uri, **kwargs))

    def get_path_or_uri(self, secret_free: bool) -> str:
        return self.path_or_uri

    def get_filename(self) -> str:
        return os.path.basename(self.path_or_uri)

    def is_persistently_cacheable(self) -> bool:
        return False

    @property
    def is_local(self) -> bool:
        return False


class LocalSourceFile(SourceFile):
    def __init__(self, path):
        self.path = path

    def endswith(self, suffix: str) -> bool:
        return self.path.endswith(suffix)

    def replace_suffix(
        self, suffix: List[str], replacement: str
    ) -> Optional["typing.Self"]:
        repl = _replace_suffix(self.path, suffix, replacement)
        if repl is None:
            return None
        else:
            return self.__class__(repl)

    def check(self) -> None:
        check(self.path)

    def format(self, **kwargs: Any) -> "typing.Self":
        return self.__class__(format(self.path, **kwargs))

    def get_path_or_uri(self, secret_free: bool) -> str:
        return self.path

    def is_persistently_cacheable(self):
        return False

    def get_filename(self) -> str:
        return os.path.basename(self.path)

    def abspath(self):
        return LocalSourceFile(os.path.abspath(self.path))

    def isabs(self):
        return os.path.isabs(self.path)

    def simplify_path(self):
        return utils.simplify_path(self.path)

    def mtime(self) -> Optional[float]:
        return os.stat(self.path).st_mtime

    def __fspath__(self):
        return self.path

    @property
    def is_local(self) -> bool:
        return True


class LocalGitFile(SourceFile):
    def __init__(
        self,
        repo_path,
        path: str,
        tag: Optional[str] = None,
        ref: Optional[str] = None,
        commit: Optional[str] = None,
    ):
        _check_git_args(tag, ref, commit)
        self.tag = tag
        self.commit = commit
        self._ref = ref
        self.repo_path = repo_path
        self.path = path

    def endswith(self, suffix: str) -> bool:
        return self.path.endswith(suffix)

    def replace_suffix(
        self, suffix: List[str], replacement: str
    ) -> Optional["typing.Self"]:
        repl = _replace_suffix(self.path, suffix, replacement)
        if repl is None:
            return None
        else:
            return self.__class__(
                self.repo_path,
                path=repl,
                tag=self.tag,
                ref=self._ref,
                commit=self.commit,
            )

    def check(self) -> None:
        check(self.path)

    def format(self, **kwargs: Any) -> "typing.Self":
        return self.__class__(
            repo_path=format(self.repo_path, **kwargs),
            path=format(self.path, **kwargs),
            tag=_format_or_none(self.tag, **kwargs),
            ref=_format_or_none(self._ref, **kwargs),
            commit=_format_or_none(self.commit, **kwargs),
        )

    def get_path_or_uri(self, secret_free: bool) -> str:
        return f"git+file://{os.path.abspath(self.repo_path)}/{self.path}@{self.ref}"

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

    def is_persistently_cacheable(self) -> bool:
        return False

    def get_filename(self) -> str:
        return posixpath.basename(self.path)

    @property
    def ref(self):
        return self.tag or self.commit or self._ref

    @property
    def is_local(self):
        return True


class HostedGitRepo:
    def __init__(
        self,
        repo: str,
        cache_path: Path,
        auth: str,
        host: str,
    ):
        from git import Repo

        repo_url = f"https://{auth}{host}/{repo}"

        self.host: str = host
        self.repo_name: str = repo

        self.repo_clone = cache_path / host / repo

        self._existed_before = self.repo_clone.exists()

        if self._existed_before:
            self._repo = Repo(self.repo_clone)
        else:
            # lock-free cloning of the repository
            logger.info(f"Cloning {host}/{repo} to {self.repo_clone}")
            self.repo_clone.parent.mkdir(parents=True, exist_ok=True)
            tmpdir = tempfile.mkdtemp(prefix=f"{self.repo_clone}.")
            # the clone is not atomic, hence we do that in a temporary directory
            # We only want the database, thus we create a bare clone
            try:
                Repo.clone_from(repo_url, to_path=tmpdir, bare=True)
            except Exception:
                # clean up on any exception
                shutil.rmtree(tmpdir)
                raise
            try:
                # move is atomic if repo_clone does not exist, so we can safely move
                # the directory to the final location
                os.rename(tmpdir, self.repo_clone)
            except FileExistsError:
                # another process won the race
                shutil.rmtree(tmpdir)

            self._repo = Repo(self.repo_clone)

    def ref_exists(self, ref: str):
        import git

        try:
            self.repo.git.rev_parse("--verify", ref)
            return True
        except git.GitCommandError:
            return False

    def file_exists(self, path: str, ref: str):
        try:
            # Get the tree for the specific revision
            tree = self.repo.commit(ref).tree
            # Attempt to access the file by path
            tree[path]
            return True
        except KeyError:
            return False

    def fetch(self) -> Optional[str]:
        import git
        from reretry import retry_call

        logger.info(
            f"Fetching latest changes of {self.host}/{self.repo_name} to {self.repo_clone}"
        )
        try:
            retry_call(
                # arg is needed in order to have a refspec, associated issue with
                # workaround is here: https://github.com/gitpython-developers/GitPython/issues/296
                partial(self.repo.remotes.origin.fetch, "+refs/heads/*:refs/heads/*"),
                delay=3,
                backoff=2,
                tries=3,
            )
        except git.GitCommandError as e:
            return str(e)

    @property
    def repo(self) -> "git.Repo":
        return self._repo


class HostingProviderFile(SourceFile):
    """Marker for denoting github source files from releases."""

    valid_repo: ClassVar[re.Pattern] = re.compile("^.+/.+$")
    _hosted_repos: ClassVar[Dict[str, HostedGitRepo]] = {}
    _lock: ClassVar[threading.Lock] = threading.Lock()

    def __init__(
        self,
        repo: Optional[str] = None,
        path: Optional[str] = None,
        tag: Optional[str] = None,
        branch: Optional[str] = None,
        commit: Optional[str] = None,
        host: Optional[str] = None,
        cache_path: Optional[Path] = None,
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
        self.host = host

        self._cache_path: Optional[Path] = cache_path

        # Via __post_init__ implementing subclasses can do additional things without
        # replicating the constructor args.
        self.__post_init__()

    def __post_init__(self):
        pass

    def endswith(self, suffix: str) -> bool:
        return self.path.endswith(suffix)

    def replace_suffix(
        self, suffix: List[str], replacement: str
    ) -> Optional["typing.Self"]:
        repl = _replace_suffix(self.path, suffix, replacement)
        if repl is None:
            return None
        else:
            return self.__class__(
                repo=self.repo,
                path=repl,
                tag=self.tag,
                branch=self.branch,
                commit=self.commit,
                host=self.host,
                cache_path=self._cache_path,
            )

    def check(self) -> None:
        check(self.path)

    def format(self, **kwargs: Any) -> "typing.Self":
        return self.__class__(
            repo=_format_or_none(self.repo, **kwargs),
            path=_format_or_none(self.path, **kwargs),
            tag=_format_or_none(self.tag, **kwargs),
            branch=_format_or_none(self.branch, **kwargs),
            commit=_format_or_none(self.commit, **kwargs),
            host=_format_or_none(self.host, **kwargs),
            cache_path=self._cache_path,
        )

    @property
    def auth(self) -> str:
        return ""

    @property
    def cache_path(self) -> Path:
        assert (
            self._cache_path
        ), "bug: cache_path not set, should be done by SourceCache"

        return self._cache_path

    @cache_path.setter
    def cache_path(self, cache_path: Path):
        self._cache_path = cache_path / "snakemake-git-cache"

    @property
    def hosted_repo(self) -> HostedGitRepo:
        try:
            return self._hosted_repos[self.repo]
        except KeyError:
            assert self.host is not None
            if self.host.startswith("https://") or self.host.startswith("http://"):
                raise WorkflowError(
                    "host must be given as domain name without protocol prefix "
                    f"(e.g. github.com, but found {self.host})  for git-hosted source "
                    f"file {self.path} in repo {self.repo}."
                )

            # Ensure that multiple threads don't concurrently create the same instance
            with self._lock:
                hosted_repo = HostedGitRepo(
                    self.repo, self.cache_path, self.auth, self.host
                )
                self._hosted_repos[self.repo] = hosted_repo
            return hosted_repo

    def is_persistently_cacheable(self):
        return bool(self.tag or self.commit)

    def get_filename(self):
        return os.path.basename(self.path)

    def mtime(self) -> float:
        last_commit = next(
            self.hosted_repo.repo.iter_commits(
                rev=self.ref, paths=self.path, max_count=1
            )
        )
        return last_commit.committed_date

    def open(self) -> io.BytesIO:
        import git

        fetch_error = None

        if not self.hosted_repo.ref_exists(
            self.ref
        ) or not self.hosted_repo.file_exists(path=self.path, ref=self.ref):
            fetch_error = self.hosted_repo.fetch()

        try:
            return io.BytesIO(
                self.hosted_repo.repo.git.show(f"{self.ref}:{self.path}").encode()
            )
        except git.GitCommandError as e:
            msg = f"Failed to get cached git source file {self.repo}:{self.path}: {e}. "
            if fetch_error:
                msg += f" Unable to fetch from remote: {fetch_error}."
            raise WorkflowError(msg) from e

    @property
    def ref(self) -> str:
        ref = self.tag or self.commit or self.branch
        assert ref is not None
        return ref

    def get_basedir(self):
        return self.__class__(
            repo=self.repo,
            path=os.path.dirname(self.path),
            tag=self.tag,
            commit=self.commit,
            branch=self.branch,
            host=self.host,
            cache_path=self._cache_path,
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
            host=self.host,
            cache_path=self._cache_path,
        )

    @property
    def is_local(self):
        return False

    def __str__(self) -> str:
        return f"{self.host}/{self.repo}/{self.path}@{self.ref}"


class GithubFile(HostingProviderFile):
    def __post_init__(self):
        if self.host is None:
            self.host = "github.com"
        self.token = os.environ.get("GITHUB_TOKEN", "")

    @property
    def auth(self) -> str:
        if self.token:
            return f"{self.token}@"
        return ""

    def get_path_or_uri(self, secret_free: bool) -> str:
        auth = f":{self.token}@" if self.token and not secret_free else ""
        # TODO find out how this URL looks like with Github enterprise server and support
        # self.host being not none by removing the check in __post_init__
        return f"https://{auth}raw.githubusercontent.com/{self.repo}/{self.ref}/{self.path}"


class GitlabFile(HostingProviderFile):
    def __post_init__(self):
        if self.host is None:
            self.host = "gitlab.com"
        self.token = os.environ.get("GITLAB_TOKEN", "")

    @property
    def auth(self) -> str:
        if self.token:
            return f"{self.token}@"
        return ""

    def get_path_or_uri(self, secret_free: bool) -> str:
        from urllib.parse import quote

        auth = f"&private_token={self.token}" if self.token and not secret_free else ""
        return "https://{}/api/v4/projects/{}/repository/files/{}/raw?ref={}{}".format(
            self.host,
            quote(self.repo, safe=""),
            quote(self.path, safe=""),
            self.ref,
            auth,
        )


def infer_source_file(path_or_uri, basedir: Optional[SourceFile] = None) -> SourceFile:
    if isinstance(path_or_uri, SourceFile):
        if basedir is None or isinstance(path_or_uri, HostingProviderFile):
            return path_or_uri
        else:
            path_or_uri = path_or_uri.get_path_or_uri(secret_free=True)
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

    def __init__(self, cache_path: Path, runtime_cache_path: Optional[Path] = None):
        self.runtime_cache: Optional[tempfile.TemporaryDirectory]
        self.cache_path: Path = cache_path
        os.makedirs(self.cache_path, exist_ok=True)
        if runtime_cache_path is None:
            runtime_cache_parent = self.cache_path / "snakemake-runtime-cache"
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
        if isinstance(source_file, HostingProviderFile):
            source_file.cache_path = self.cache_path

        file_cache_path = source_file.get_cache_path()
        assert file_cache_path

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
        mtime = source_file.mtime()
        # open from origin
        with self._open_local_or_remote(source_file, "rb", retries=retries) as source:
            cache_entry.parent.mkdir(parents=True, exist_ok=True)
            with LockFreeWritableFile(cache_entry, binary=True) as entryfile:
                entryfile.chmod(
                    stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP
                )
                if mtime is not None:
                    # Set to mtime of original file
                    # In case we don't have that mtime, it is fine
                    # to just keep the time at the time of caching
                    # as mtime.
                    entryfile.utime((mtime, mtime))
                entryfile.write_from_fileobj(source)

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

        log_path = source_file.get_path_or_uri(secret_free=True)

        if isinstance(source_file, LocalGitFile):
            import git

            try:
                return io.BytesIO(
                    git.Repo(source_file.repo_path)
                    .git.show(f"{source_file.ref}:{source_file.path}")
                    .encode()
                )
            except git.GitCommandError as e:
                raise WorkflowError(
                    f"Failed to get local git source file {log_path}: {e}. "
                    "Is the local git clone up to date?"
                )
        elif isinstance(source_file, HostingProviderFile):
            return source_file.open()

        path_or_uri = source_file.get_path_or_uri(secret_free=False)

        try:
            return open(path_or_uri, mode, encoding=None if "b" in mode else encoding)
        except Exception as e:
            raise WorkflowError(
                f"Failed to open source file {log_path}",
                e,
            )
