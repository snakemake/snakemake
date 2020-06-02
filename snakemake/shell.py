__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import _io
import sys
import os
import subprocess as sp
import inspect
import shutil
import threading

from snakemake.utils import format, argvquote, _find_bash_on_windows
from snakemake.common import ON_WINDOWS
from snakemake.logging import logger
from snakemake.deployment import singularity
from snakemake.deployment.conda import Conda
from snakemake.exceptions import WorkflowError
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
    _win_command_prefix = ""

    @classmethod
    def get_executable(cls):
        return cls._process_args.get("executable", None)

    @classmethod
    def check_output(cls, cmd, **kwargs):
        executable = cls.get_executable()
        if ON_WINDOWS and executable:
            cmd = f'"{executable}" {cls._win_command_prefix} {argvquote(cmd)}'
            return sp.check_output(cmd, shell=False, executable=executable, **kwargs,)
        else:
            return sp.check_output(cmd, shell=True, executable=executable, **kwargs)

    @classmethod
    def executable(cls, cmd):
        if os.name in ("posix", "nt") and not os.path.isabs(cmd):
            # always enforce absolute path
            cmd = shutil.which(cmd)
            if not cmd:
                raise WorkflowError(
                    "Cannot set default shell {} because it "
                    "is not available in your "
                    "PATH.".format(cmd)
                )
        if os.path.split(cmd)[-1].lower() in ("bash", "bash.exe"):
            cls._process_prefix = "set -euo pipefail; "
            cls._win_command_prefix = "-c"
        cls._process_args["executable"] = cmd

    @classmethod
    def prefix(cls, prefix):
        cls._process_prefix = format(prefix, stepout=2)

    @classmethod
    def suffix(cls, suffix):
        cls._process_suffix = format(suffix, stepout=2)

    @classmethod
    def win_command_prefix(cls, cmd):
        """ The command prefix used on windows when specifing a explicit 
            shell executable. This would be "-c" for bash and "/C" for cmd.exe
            Note: that if no explicit executable is set commands are executed 
            with Popen(..., shell=True) which uses COMSPEC on windows where this
            is not needed.
        """
        cls._win_command_prefix = cmd

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

    def __new__(
        cls, cmd, *args, iterable=False, read=False, bench_record=None, **kwargs
    ):
        if "stepout" in kwargs:
            raise KeyError("Argument stepout is not allowed in shell command.")
        cmd = format(cmd, *args, stepout=2, **kwargs)
        context = inspect.currentframe().f_back.f_locals
        # add kwargs to context (overwriting the locals of the caller)
        context.update(kwargs)

        stdout = sp.PIPE if iterable or read else STDOUT

        close_fds = sys.platform != "win32"

        jobid = context.get("jobid")
        if not context.get("is_shell"):
            logger.shellcmd(cmd)

        env_prefix = ""
        conda_env = context.get("conda_env", None)
        container_img = context.get("container_img", None)
        env_modules = context.get("env_modules", None)
        shadow_dir = context.get("shadow_dir", None)

        cmd = "{} {} {}".format(
            cls._process_prefix, cmd.strip(), cls._process_suffix
        ).strip()

        if env_modules:
            cmd = env_modules.shellcmd(cmd)
            logger.info("Activating environment modules: {}".format(env_modules))

        if conda_env:
            cmd = Conda(container_img).shellcmd(conda_env, cmd)

        if container_img:
            args = context.get("singularity_args", "")
            cmd = singularity.shellcmd(
                container_img,
                cmd,
                args,
                shell_executable=cls._process_args["executable"],
                container_workdir=shadow_dir,
            )
            logger.info("Activating singularity image {}".format(container_img))
        if conda_env:
            logger.info("Activating conda environment: {}".format(conda_env))

        threads = str(context.get("threads", 1))
        # environment variable lists for linear algebra libraries taken from:
        # https://stackoverflow.com/a/53224849/2352071
        # https://github.com/xianyi/OpenBLAS/tree/59243d49ab8e958bb3872f16a7c0ef8c04067c0a#setting-the-number-of-threads-using-environment-variables
        envvars = dict(os.environ)
        envvars["OMP_NUM_THREADS"] = threads
        envvars["GOTO_NUM_THREADS"] = threads
        envvars["OPENBLAS_NUM_THREADS"] = threads
        envvars["MKL_NUM_THREADS"] = threads
        envvars["VECLIB_MAXIMUM_THREADS"] = threads
        envvars["NUMEXPR_NUM_THREADS"] = threads

        use_shell = True

        if ON_WINDOWS and cls.get_executable():
            # If executable is set on Windows shell mode can not be used
            # and the executable should be prepended the command together
            # with a command prefix (e.g. -c for bash).
            use_shell = False
            cmd = f'"{cls.get_executable()}" {cls._win_command_prefix} {argvquote(cmd)}'

        proc = sp.Popen(
            cmd,
            bufsize=-1,
            shell=use_shell,
            stdout=stdout,
            universal_newlines=iterable or read or None,
            close_fds=close_fds,
            **cls._process_args,
            env=envvars,
        )

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
            yield l[:-1]
        retcode = proc.wait()
        if retcode:
            raise sp.CalledProcessError(retcode, cmd)


# set bash as default shell on posix compatible OS
if os.name == "posix":
    if not shutil.which("bash"):
        logger.warning(
            "Cannot set bash as default shell because it is not "
            "available in your PATH. Falling back to sh."
        )
        if not shutil.which("sh"):
            logger.warning(
                "Cannot fall back to sh since it seems to be not "
                "available on this system. Using whatever is "
                "defined as default."
            )
        else:
            shell.executable("sh")
    else:
        shell.executable("bash")
