__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from pathlib import Path
import _io
import sys
import os
import subprocess as sp
import inspect
import shutil
import stat
import tempfile
import threading

from snakemake.utils import format, argvquote, cmd_exe_quote
from snakemake.common import ON_WINDOWS, RULEFUNC_CONTEXT_MARKER
from snakemake.logging import logger
from snakemake.deployment import singularity
from snakemake.deployment.conda import Conda
from snakemake.exceptions import WorkflowError


__author__ = "Johannes Köster"

STDOUT = sys.stdout
if not isinstance(sys.stdout, _io.TextIOWrapper):
    # workaround for nosetest since it overwrites sys.stdout
    # in a strange way that does not work with Popen
    STDOUT = None


# There is a max length for a command executed as well as a maximum
# length for each argument passed to a command. The latter impacts us
# especially when doing `sh -c 'long script from user'`. On Linux, it's
# hardcoded in the kernel as 32 pages, or 128kB. On OSX it appears to be
# close to `getconf ARG_MAX`, about 253kb.
MAX_ARG_LEN = 16 * 4096 - 1


class shell:
    _process_args = {}
    _process_prefix = None
    _process_suffix = ""
    _win_command_prefix = ""
    _lock = threading.Lock()
    _processes = {}
    conda_block_conflicting_envvars = True

    @classmethod
    def get_executable(cls):
        return cls._process_args.get("executable", None)

    @classmethod
    def check_output(cls, cmd, **kwargs):
        executable = cls.get_executable()
        if ON_WINDOWS and executable:
            win_prefix = cls._get_win_command_prefix()
            cmd = f'"{executable}" {win_prefix} {argvquote(cmd)}'
            logger.debug(f"Executing: {cmd}")
            return sp.check_output(cmd, shell=False, executable=executable, **kwargs)
        else:
            return sp.check_output(cmd, shell=True, executable=executable, **kwargs)

    @classmethod
    def executable(cls, cmd):
        if isinstance(cmd, Path):
            cmd = str(cmd)
        if cmd and not os.path.isabs(cmd):
            # always enforce absolute path
            cmd = shutil.which(cmd)
            if not cmd:
                raise WorkflowError(
                    f"Cannot set default shell {cmd} because it is not available in your PATH."
                )
        cls._process_args["executable"] = cmd
        logger.debug(f"Setting shell executable to {cmd}.")

    @classmethod
    def _get_process_prefix(cls, shell_exec=None):
        shell_exec = cls._get_executable_name(shell_exec)
        if (
            shell_exec == "bash" or (ON_WINDOWS and shell_exec == "bash.exe")
        ) and cls._process_prefix is None:
            return "set -euo pipefail; "
        else:
            return cls._process_prefix or ""

    @classmethod
    def _get_win_command_prefix(cls, use_default=False, shell_exec=None):
        assert ON_WINDOWS
        if use_default or (cls._win_command_prefix and shell_exec is None):
            # use whatever is the default
            return cls._win_command_prefix
        shell_exec = cls._get_executable_name(shell_exec)
        if shell_exec == "bash" or shell_exec == "bash.exe":
            return "-c"
        else:
            return ""

    @classmethod
    def _check_executable(cls, shell_exec=None):
        shell_exec = shell_exec or cls.get_executable()
        if shell_exec is not None:
            if ON_WINDOWS and shell_exec == r"C:\Windows\System32\bash.exe":
                raise WorkflowError(
                    "Cannot use WSL bash.exe on Windows. Ensure that you have "
                    "a usable bash.exe available on your path."
                )
            if not os.path.isabs(shell_exec):
                path = shutil.which(shell_exec)
                if not path:
                    raise WorkflowError(
                        f"Cannot set shell to {shell_exec} because it is not "
                        "available in your PATH."
                    )
            elif not os.path.exists(shell_exec):
                raise WorkflowError(
                    f"Cannot set shell to {shell_exec} because it does not exist."
                )

    @classmethod
    def _get_executable_name(cls, shell_exec=None):
        shell_exec = shell_exec or cls.get_executable()
        if shell_exec:
            return os.path.split(shell_exec)[-1].lower()
        else:
            return None

    @classmethod
    def prefix(cls, prefix):
        cls._process_prefix = format(prefix, stepout=2)

    @classmethod
    def suffix(cls, suffix):
        cls._process_suffix = format(suffix, stepout=2)

    @classmethod
    def win_command_prefix(cls, cmd):
        """The command prefix used on windows when specifying a explicit
        shell executable. This would be "-c" for bash.
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
    def terminate(cls, jobid):
        with cls._lock:
            if jobid in cls._processes:
                cls._processes[jobid].terminate()
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

        if ON_WINDOWS and not cls.get_executable():
            # If bash is not used on Windows quoting must be handled in a special way
            kwargs["quote_func"] = cmd_exe_quote

        cmd = format(cmd, *args, stepout=2, **kwargs)

        stdout = sp.PIPE if iterable or read else STDOUT

        close_fds = sys.platform != "win32"

        func_context = inspect.currentframe().f_back.f_locals

        if func_context.get(RULEFUNC_CONTEXT_MARKER):
            # If this comes from a rule, we expect certain information to be passed
            # implicitly via the rule func context, which is added here.
            context = func_context
        else:
            # Otherwise, context is just filled via kwargs.
            context = dict()
        # add kwargs to context (overwriting the locals of the caller)
        context.update(kwargs)

        jobid = context.get("jobid")
        if not context.get("is_shell") and jobid is not None:
            logger.shellcmd(cmd)

        conda_env = context.get("conda_env", None)
        conda_base_path = context.get("conda_base_path", None)
        container_img = context.get("container_img", None)
        env_modules = context.get("env_modules", None)
        shadow_dir = context.get("shadow_dir", None)
        resources = context.get("resources", {})
        singularity_args = context.get("singularity_args", "")
        threads = context.get("threads", 1)

        shell_executable = resources.get("shell_exec")
        if shell_executable is not None:
            process_args = dict(cls._process_args)
            process_args["executable"] = shell_executable
        else:
            shell_executable = cls.get_executable()
            process_args = cls._process_args

        cls._check_executable(shell_executable)

        cmd = " ".join(
            (cls._get_process_prefix(shell_executable), cmd, cls._process_suffix)
        ).strip()

        # If the executor is the submit executor or the jobstep executor for the SLURM
        # backend, we do not want the environment modules to be activated:
        # if the rule requires a Python module, snakemake's environment might be
        # incompatible with the module's environment.
        if env_modules and "slurm" not in (item.filename for item in inspect.stack()):
            cmd = env_modules.shellcmd(cmd)
            logger.info(f"Activating environment modules: {env_modules}")

        if conda_env:
            if ON_WINDOWS and not cls.get_executable():
                # If we use cmd.exe directly on windows we need to prepend batch activation script.
                cmd = Conda(
                    container_img=container_img, prefix_path=conda_base_path
                ).shellcmd_win(conda_env, cmd)
            else:
                cmd = Conda(
                    container_img=container_img, prefix_path=conda_base_path
                ).shellcmd(conda_env, cmd)

        tmpdir = None
        if len(cmd.replace("'", r"'\''")) + 2 > MAX_ARG_LEN:
            tmpdir = tempfile.mkdtemp(dir=".snakemake", prefix="shell_tmp.")
            script = os.path.join(os.path.abspath(tmpdir), "script.sh")
            with open(script, "w") as script_fd:
                print(cmd, file=script_fd)
            os.chmod(script, os.stat(script).st_mode | stat.S_IXUSR | stat.S_IRUSR)
            cmd = '"{}" "{}"'.format(cls.get_executable() or "/bin/sh", script)

        if container_img:
            cmd = singularity.shellcmd(
                container_img,
                cmd,
                singularity_args,
                envvars=None,
                shell_executable=shell_executable,
                container_workdir=shadow_dir,
                is_python_script=context.get("is_python_script", False),
            )
            logger.info(f"Activating singularity image {container_img}")
        if conda_env:
            logger.info(f"Activating conda environment: {os.path.relpath(conda_env)}")

        tmpdir_resource = resources.get("tmpdir", None)
        # environment variable lists for linear algebra libraries taken from:
        # https://stackoverflow.com/a/53224849/2352071
        # https://github.com/xianyi/OpenBLAS/tree/59243d49ab8e958bb3872f16a7c0ef8c04067c0a#setting-the-number-of-threads-using-environment-variables
        envvars = dict(os.environ)
        threads = str(threads)
        envvars["OMP_NUM_THREADS"] = threads
        envvars["GOTO_NUM_THREADS"] = threads
        envvars["OPENBLAS_NUM_THREADS"] = threads
        envvars["MKL_NUM_THREADS"] = threads
        envvars["VECLIB_MAXIMUM_THREADS"] = threads
        envvars["NUMEXPR_NUM_THREADS"] = threads

        if tmpdir_resource:
            envvars["TMPDIR"] = tmpdir_resource
            envvars["TMP"] = tmpdir_resource
            envvars["TEMPDIR"] = tmpdir_resource
            envvars["TEMP"] = tmpdir_resource

        if "additional_envvars" in kwargs:
            env = kwargs["additional_envvars"]
            if not isinstance(env, dict) or not all(
                isinstance(v, str) for v in env.values()
            ):
                raise WorkflowError(
                    "Given environment variables for shell command have to be a dict of strings, "
                    "but the following was provided instead:\n{}".format(env)
                )
            envvars.update(env)

        if conda_env and cls.conda_block_conflicting_envvars:
            # remove envvars that conflict with conda
            for var in ["R_LIBS", "PYTHONPATH", "PERLLIB", "PERL5LIB"]:
                try:
                    del envvars[var]
                except KeyError:
                    pass

        use_shell = True
        if ON_WINDOWS and shell_executable:
            # If executable is set on Windows shell mode can not be used
            # and the executable should be prepended the command together
            # with a command prefix (e.g. -c for bash).
            use_shell = False
            win_prefix = cls._get_win_command_prefix(
                use_default=False, shell_exec=shell_executable
            )
            cmd = '"{}" {} {}'.format(
                shell_executable,
                win_prefix,
                argvquote(cmd),
            )

        proc = sp.Popen(
            cmd,
            bufsize=-1,
            shell=use_shell,
            stdout=stdout,
            universal_newlines=iterable or read or None,
            close_fds=close_fds,
            **process_args,
            env=envvars,
        )

        if jobid is not None:
            with cls._lock:
                cls._processes[jobid] = proc

        ret = None
        if iterable:
            return cls.iter_stdout(proc, cmd, tmpdir)
        if read:
            ret = proc.stdout.read()
        if bench_record is not None:
            from snakemake.benchmark import benchmarked

            with benchmarked(proc.pid, bench_record):
                retcode = proc.wait()
        else:
            retcode = proc.wait()

        if tmpdir:
            shutil.rmtree(tmpdir)

        if jobid is not None:
            with cls._lock:
                try:
                    del cls._processes[jobid]
                except KeyError:
                    pass

        if retcode:
            raise sp.CalledProcessError(retcode, cmd)
        return ret

    @staticmethod
    def iter_stdout(proc, cmd, tmpdir):
        for l in proc.stdout:
            yield l[:-1]
        retcode = proc.wait()
        if tmpdir:
            shutil.rmtree(tmpdir)
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
elif ON_WINDOWS:
    shell.executable(None)
