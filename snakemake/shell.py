__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import _io
import sys
import os
import subprocess as sp
import inspect
import shutil
import threading

from snakemake.utils import format
from snakemake.logging import logger
from snakemake import singularity, conda
import snakemake


__author__ = "Johannes Köster"

STDOUT = sys.stdout
if not isinstance(sys.stdout, _io.TextIOWrapper):
    # workaround for nosetest since it overwrites sys.stdout
    # in a strange way that does not work with Popen
    STDOUT = None


class shell:
    _process_args = {}
    _process_prefix = ""
    _process_suffix = ""
    _lock = threading.Lock()
    _processes = {}

    @classmethod
    def get_executable(cls):
        return cls._process_args.get("executable", None)

    @classmethod
    def executable(cls, cmd):
        if os.name == "posix" and not os.path.isabs(cmd):
            # always enforce absolute path
            cmd = shutil.which(cmd)
            if not cmd:
                raise WorkflowError("Cannot set default shell {} because it "
                                    "is not available in your "
                                    "PATH.".format(cmd))
        if os.path.split(cmd)[-1] == "bash":
            cls._process_prefix = "set -euo pipefail; "
        cls._process_args["executable"] = cmd

    @classmethod
    def prefix(cls, prefix):
        cls._process_prefix = format(prefix, stepout=2)

    @classmethod
    def suffix(cls, suffix):
        cls._process_suffix = format(suffix, stepout=2)

    @classmethod
    def kill(cls, jobid):
        with cls._lock:
            if jobid in cls._processes:
                cls._processes[jobid].kill()
                del cls._processes[jobid]

    @classmethod
    def cleanup(cls):
        with cls._lock:
            cls._processes.clear()

    def __new__(cls, cmd, *args,
                iterable=False,
                read=False, bench_record=None,
                **kwargs):
        if "stepout" in kwargs:
            raise KeyError("Argument stepout is not allowed in shell command.")
        cmd = format(cmd, *args, stepout=2, **kwargs)
        context = inspect.currentframe().f_back.f_locals

        stdout = sp.PIPE if iterable or read else STDOUT

        close_fds = sys.platform != 'win32'

        jobid = context.get("jobid")
        if jobid is not None:
            if not context.get("is_shell"):
                logger.shellcmd(cmd)

        env_prefix = ""
        conda_env = context.get("conda_env", None)
        singularity_img = context.get("singularity_img", None)
        if conda_env:
            env_prefix = conda.shellcmd(conda_env)

        cmd = "{} {} {} {}".format(
                            env_prefix,
                            cls._process_prefix,
                            cmd.rstrip(),
                            cls._process_suffix)

        if singularity_img:
            args = context.get("singularity_args", "")
            cmd = singularity.shellcmd(singularity_img, cmd, args)
            logger.info("Activating singularity image {}".format(singularity_img))

        if conda_env:
            logger.info("Activating conda environment: {}".format(conda_env))

        proc = sp.Popen(cmd,
                        bufsize=-1,
                        shell=True,
                        stdout=stdout,
                        close_fds=close_fds, **cls._process_args)

        if jobid is not None:
            with cls._lock:
                cls._processes[jobid] = proc

        ret = None
        if iterable:
            return cls.iter_stdout(proc, cmd)
        if read:
            ret = proc.stdout.read()
        if bench_record is not None:
            from snakemake.benchmark import benchmarked
            with benchmarked(proc.pid, bench_record):
                retcode = proc.wait()
        else:
            retcode = proc.wait()

        if jobid is not None:
            with cls._lock:
                del cls._processes[jobid]

        if retcode:
            raise sp.CalledProcessError(retcode, cmd)
        return ret

    @staticmethod
    def iter_stdout(proc, cmd):
        for l in proc.stdout:
            yield l[:-1].decode()
        retcode = proc.wait()
        if retcode:
            raise sp.CalledProcessError(retcode, cmd)


# set bash as default shell on posix compatible OS
if os.name == "posix":
    if not shutil.which("bash"):
        logger.warning("Cannot set bash as default shell because it is not "
                       "available in your PATH. Falling back to sh.")
        if not shutil.which("sh"):
            logger.warning("Cannot fall back to sh since it seems to be not "
                           "available on this system. Using whatever is "
                           "defined as default.")
        else:
            shell.executable("sh")
    else:
        shell.executable("bash")
