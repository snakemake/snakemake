# -*- coding: utf-8 -*-

import _io
import sys
import os
import subprocess as sp

from snakemake.utils import format
from snakemake.logging import logger


__author__ = "Johannes Köster"


STDOUT = sys.stdout
if not isinstance(sys.stdout, _io.TextIOWrapper):
    # workaround for nosetest since it overwrites sys.stdout
    # in a strange way that does not work with Popen
    STDOUT = None


class shell:
    _process_args = {}
    _process_prefix = ""

    @classmethod
    def executable(cls, cmd):
        cls._process_args["executable"] = cmd

    @classmethod
    def prefix(cls, prefix):
        cls._process_prefix = format(prefix, stepout=2)

    def __new__(
        cls, cmd, *args, async=False, iterable=False, read=False, **kwargs):
        if "stepout" in kwargs:
            raise KeyError("Argument stepout is not allowed in shell command.")
        cmd = format(cmd, *args, stepout=2, **kwargs)

        logger.shellcmd(cmd)

        stdout = sp.PIPE if iterable or async or read else STDOUT

        close_fds = sys.platform != 'win32'
        proc = sp.Popen(cls._process_prefix + cmd, bufsize=-1, shell=True, stdout=stdout,
            close_fds=close_fds, **cls._process_args)

        ret = None
        if iterable:
            return cls.iter_stdout(proc, cmd)
        if read:
            ret = proc.stdout.read()
        elif async:
            return proc
        retcode = proc.wait()
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


if "SHELL" in os.environ:
    shell.executable(os.environ["SHELL"])
