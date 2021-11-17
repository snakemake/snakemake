__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import hashlib
from pathlib import Path
import re
import os
from snakemake import utils
import tempfile
import io
from abc import ABC, abstractmethod


from snakemake.common import is_local_file, get_appdirs, parse_uri, smart_join
from snakemake.exceptions import WorkflowError, SourceFileError
from snakemake.io import git_content, split_git_path
from snakemake.logging import logger


def _check_git_args(tag: str = None, branch: str = None, commit: str = None):
    n_refs = sum(1 for ref in (tag, branch, commit) if ref is not None)
    if n_refs != 1:
        raise SourceFileError(
            "exactly one of tag, branch, or commit must be specified."
        )


class SourceFile(ABC):
    @abstractmethod
    def get_path_or_uri(self):
        ...

    @abstractmethod
    def is_persistently_cacheable(self):
        ...

    def get_uri_hash(self):
        urihash = hashlib.sha256()
        urihash.update(self.get_path_or_uri().encode())
        return urihash.hexdigest()

    def get_basedir(self):
        path = os.path.dirname(self.get_path_or_uri())
        return self.__class__(path)

    @abstractmethod
    def get_filename(self):
        ...

    def join(self, path):
        if isinstance(path, SourceFile):
            path = path.get_path_or_uri()
        return self.__class__(smart_join(self.get_path_or_uri(), path))

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

    def __fspath__(self):
        return self.path


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
        return "git+{}/{}@{}".format(self.repo_path, self.path, self.ref)

    def join(self, path):
        return LocalGitFile(
            self.repo_path,
            "/".join((self.path, path)),
            tag=self.tag,
            ref=self.ref,
            commit=self.commit,
        )

    def is_persistently_cacheable(self):
        return False

    def get_filename(self):
        return os.path.basename(self.path)

    @property
    def ref(self):
        return self.tag or self.commit or self._ref


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

    def is_persistently_cacheable(self):
        return self.tag or self.commit

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
        path = os.path.normpath("{}/{}".format(self.path, path))
        return self.__class__(
            repo=self.repo,
            path=path,
            tag=self.tag,
            commit=self.commit,
            branch=self.branch,
        )


class GithubFile(HostingProviderFile):
    def get_path_or_uri(self):
        return "https://github.com/{}/raw/{}/{}".format(self.repo, self.ref, self.path)


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

    def get_path_or_uri(self):
        return "https://{}/{}/-/raw/{}/{}".format(
            self.host or "gitlab.com", self.repo, self.ref, self.path
        )


def infer_source_file(path_or_uri, basedir: SourceFile = None):
    if isinstance(path_or_uri, SourceFile):
        if basedir is None:
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
        root_path, file_path, ref = split_git_path(path_or_uri)
        return LocalGitFile(root_path, file_path, ref=ref)
    # something else
    return GenericSourceFile(path_or_uri)


class SourceCache:
    cache_whitelist = [
        "https://raw.githubusercontent.com/snakemake/snakemake-wrappers/\d+\.\d+.\d+"
    ]  # TODO add more prefixes for uris that are save to be cached

    def __init__(self, runtime_cache_path=None):
        self.cache = Path(
            os.path.join(get_appdirs().user_cache_dir, "snakemake/source-cache")
        )
        os.makedirs(self.cache, exist_ok=True)
        if runtime_cache_path is None:
            self.runtime_cache = tempfile.TemporaryDirectory(
                suffix="snakemake-runtime-source-cache"
            )
            self._runtime_cache_path = None
        else:
            self._runtime_cache_path = runtime_cache_path
            self.runtime_cache = None
        self.cacheable_prefixes = re.compile("|".join(self.cache_whitelist))

    @property
    def runtime_cache_path(self):
        return self._runtime_cache_path or self.runtime_cache.name

    def lock_cache(self, entry):
        from filelock import FileLock

        return FileLock(entry.with_suffix(".lock"))

    def open(self, source_file, mode="r"):
        cache_entry = self._cache(source_file)
        return self._open(cache_entry, mode)

    def exists(self, source_file):
        try:
            self._cache(source_file)
        except Exception:
            return False
        return True

    def get_path(self, source_file, mode="r"):
        cache_entry = self._cache(source_file)
        return cache_entry

    def _cache_entry(self, source_file):
        urihash = source_file.get_uri_hash()

        # TODO add git support to smart_open!
        if source_file.is_persistently_cacheable():
            # check cache
            return self.cache / urihash
        else:
            # check runtime cache
            return Path(self.runtime_cache_path) / urihash

    def _cache(self, source_file):
        cache_entry = self._cache_entry(source_file)
        with self.lock_cache(cache_entry):
            if not cache_entry.exists():
                self._do_cache(source_file, cache_entry)
        return cache_entry

    def _do_cache(self, source_file, cache_entry):
        # open from origin
        with self._open(source_file.get_path_or_uri(), "rb") as source, open(
            cache_entry, "wb"
        ) as cache_source:
            cache_source.write(source.read())

    def _open(self, path_or_uri, mode):
        from smart_open import open

        if isinstance(path_or_uri, LocalGitFile):
            import git

            return io.BytesIO(
                git.Repo(path_or_uri.repo_path)
                .git.show("{}:{}".format(path_or_uri.ref, path_or_uri.path))
                .encode()
            )

        if isinstance(path_or_uri, SourceFile):
            path_or_uri = path_or_uri.get_path_or_uri()

        try:
            return open(path_or_uri, mode)
        except Exception as e:
            raise WorkflowError("Failed to open source file {}".format(path_or_uri), e)
