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

from snakemake.utils import format, ON_WINDOWS, argvquote, find_bash_on_windows
from snakemake.logging import logger
from snakemake import singularity
from snakemake.conda import Conda
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
    _quote_win_cmd = False

    @classmethod
    def get_executable(cls):
        return cls._process_args.get("executable", None)

    @classmethod
    def check_output(cls, cmd, **kwargs):
        if ON_WINDOWS and cls.get_executable():
            return sp.check_output(
                cls._process_prefix + " "+ argvquote(cmd),
                shell=False,
                executable=cls.get_executable(),
                **kwargs
            )
        else:
            return sp.check_output(
                cmd, shell=True, executable=cls.get_executable(), **kwargs)

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
        if ON_WINDOWS and cmd.lower().endswith('bash.exe'):
            cls._process_prefix = '"{}" -c'.format(cmd)
            cls.win_quote_cmd(True)
        cls._process_args["executable"] = cmd

    @classmethod
    def prefix(cls, prefix):
        cls._process_prefix = format(prefix, stepout=2)

    @classmethod
    def suffix(cls, suffix):
        cls._process_suffix = format(suffix, stepout=2)

    @classmethod
    def win_quote_cmd(cls, opt):
        """ Quote the entire shell command as a single argument on Windows.
            This can be usefull if the executable is set (e.g. to "bash.exe"), 
            and the prefix is "-c",  because the command then has to be 
            interpreted as a single argument by CreateProcess on Windows.
        """
        cls._quote_win_cmd = bool(opt)

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

    @classmethod
    def use_bash_on_win(cls, bashcmd=None):
        """ Configures the shell to use bash on Windows
            when running shell commands. If no path to 
            bash is given it will attempt to find it.
        """
        if ON_WINDOWS:
            if bashcmd is None:
                bashcmd = find_bash_on_windows()
            if bashcmd and os.path.exists(bashcmd):
                cls.executable(bashcmd)
            else:
                raise ("Could not locate bash:" + str(bashcmd))


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
        if not context.get("is_shell"):
            logger.shellcmd(cmd)

        env_prefix = ""
        conda_env = context.get("conda_env", None)
        singularity_img = context.get("singularity_img", None)
        shadow_dir = context.get("shadow_dir", None)

        if ON_WINDOWS and cls._quote_win_cmd:
            cmd = argvquote(cmd.strip())

        cmd = "{} {} {}".format(
                            cls._process_prefix,
                            cmd.strip(),
                            cls._process_suffix).strip()

        conda = None
        if conda_env:
            cmd = Conda(singularity_img).shellcmd(conda_env, cmd)

        if singularity_img:
            args = context.get("singularity_args", "")
            cmd = singularity.shellcmd(
                singularity_img, cmd, args,
                shell_executable=cls._process_args["executable"],
                container_workdir=shadow_dir)
            logger.info(
                "Activating singularity image {}".format(singularity_img))
        if conda_env:
            logger.info("Activating conda environment: {}".format(conda_env))

        # shell can't be True on Windows with an explicit executable
        use_shell = not(ON_WINDOWS and cls.get_executable())

        proc = sp.Popen(cmd,
                        bufsize=-1,
                        shell=use_shell,
                        stdout=stdout,
                        universal_newlines=iterable,
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
            yield l[:-1]
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
