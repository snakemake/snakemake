__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import logging as _logging
import platform
import time
import sys
import os
import json
from multiprocessing import Lock
import tempfile
from functools import partial

from snakemake.common import DYNAMIC_FILL


class ColorizingStreamHandler(_logging.StreamHandler):
    _output_lock = Lock()

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[%dm"
    BOLD_SEQ = "\033[1m"

    colors = {
        'WARNING': YELLOW,
        'INFO': GREEN,
        'DEBUG': BLUE,
        'CRITICAL': RED,
        'ERROR': RED
    }

    def __init__(self, nocolor=False, stream=sys.stderr, timestamp=False):
        super().__init__(stream=stream)
        self.nocolor = nocolor or not self.can_color_tty
        self.timestamp = timestamp

    @property
    def can_color_tty(self):
        if 'TERM' in os.environ and os.environ['TERM'] == 'dumb':
            return False
        return self.is_tty and not platform.system() == 'Windows'

    @property
    def is_tty(self):
        isatty = getattr(self.stream, 'isatty', None)
        return isatty and isatty()

    def emit(self, record):
        with self._output_lock:
            try:
                self.format(record)  # add the message to the record
                self.stream.write(self.decorate(record))
                self.stream.write(getattr(self, 'terminator', '\n'))
                self.flush()
            except BrokenPipeError as e:
                raise e
            except (KeyboardInterrupt, SystemExit):
                # ignore any exceptions in these cases as any relevant messages have been printed before
                pass
            except Exception as e:
                self.handleError(record)

    def decorate(self, record):
        message = [record.message]
        if self.timestamp:
            message.insert(0, "[{}] ".format(time.asctime()))
        if not self.nocolor and record.levelname in self.colors:
            message.insert(0, self.COLOR_SEQ %
                           (30 + self.colors[record.levelname]))
            message.append(self.RESET_SEQ)
        return "".join(message)


class Logger:
    def __init__(self):
        self.logger = _logging.getLogger(__name__)
        self.log_handler = [self.text_handler]
        self.stream_handler = None
        self.printshellcmds = False
        self.printreason = False
        self.quiet = False
        self.logfile = None

    def setup(self):
        # logfile output is done always
        self.logfile_fd, self.logfile = tempfile.mkstemp(
            prefix="",
            suffix=".snakemake.log")
        self.logfile_handler = _logging.FileHandler(self.logfile)
        self.logger.addHandler(self.logfile_handler)

    def cleanup(self):
        self.logger.removeHandler(self.logfile_handler)
        self.logfile_handler.close()
        os.close(self.logfile_fd)
        os.remove(self.logfile)

    def get_logfile(self):
        if self.logfile is not None:
            self.logfile_handler.flush()
        return self.logfile

    def handler(self, msg):
        for handler in self.log_handler:
            handler(msg)

    def set_stream_handler(self, stream_handler):
        if self.stream_handler is not None:
            self.logger.removeHandler(self.stream_handler)
        self.stream_handler = stream_handler
        self.logger.addHandler(stream_handler)

    def set_level(self, level):
        self.logger.setLevel(level)

    def info(self, msg):
        self.handler(dict(level="info", msg=msg))

    def warning(self, msg):
        self.handler(dict(level="warning", msg=msg))

    def debug(self, msg):
        self.handler(dict(level="debug", msg=msg))

    def error(self, msg):
        self.handler(dict(level="error", msg=msg))

    def progress(self, done=None, total=None):
        self.handler(dict(level="progress", done=done, total=total))

    def resources_info(self, msg):
        self.handler(dict(level="resources_info", msg=msg))

    def run_info(self, msg):
        self.handler(dict(level="run_info", msg=msg))

    def job_info(self, **msg):
        msg["level"] = "job_info"
        self.handler(msg)

    def shellcmd(self, msg):
        if msg is not None:
            self.handler(dict(level="shellcmd", msg=msg))

    def job_finished(self, **msg):
        msg["level"] = "job_finished"
        self.handler(msg)

    def rule_info(self, **msg):
        msg["level"] = "rule_info"
        self.handler(msg)

    def d3dag(self, **msg):
        msg["level"] = "d3dag"
        self.handler(msg)

    def text_handler(self, msg):
        """The default snakemake log handler.

        Prints the output to the console.

        Args:
            msg (dict):     the log message dictionary
        """

        def job_info(msg):
            def format_item(item, omit=None, valueformat=str):
                value = msg[item]
                if value != omit:
                    return "\t{}: {}".format(item, valueformat(value))

            yield "{}rule {}:".format("local" if msg["local"] else "",
                                      msg["name"])
            for item in ["input", "output", "log"]:
                fmt = format_item(item, omit=[], valueformat=", ".join)
                if fmt != None:
                    yield fmt

            singleitems = ["benchmark"]
            for item in singleitems:
                fmt = format_item(item, omit=None)
                if fmt != None:
                    yield fmt

            wildcards = format_wildcards(msg["wildcards"])
            if wildcards:
                yield "\twildcards: " + wildcards

            for item, omit in zip("priority threads".split(), [0, 1]):
                fmt = format_item(item, omit=omit)
                if fmt != None:
                    yield fmt

            resources = format_resources(msg["resources"])
            if resources:
                yield "\tresources: " + resources

            if self.printreason:
                singleitems.append("reason")

        level = msg["level"]
        if level == "info":
            self.logger.warning(msg["msg"])
        if level == "warning":
            self.logger.warning(msg["msg"])
        elif level == "error":
            self.logger.error(msg["msg"])
        elif level == "debug":
            self.logger.debug(msg["msg"])
        elif level == "resources_info":
            self.logger.warning(msg["msg"])
        elif level == "run_info":
            self.logger.warning(msg["msg"])
        elif level == "progress" and not self.quiet:
            done = msg["done"]
            total = msg["total"]
            self.logger.info("{} of {} steps ({:.0%}) done".format(
                done, total, done / total))
        elif level == "job_info":
            if not self.quiet:
                if msg["msg"] is not None:
                    self.logger.info(msg["msg"])
                else:
                    self.logger.info("\n".join(job_info(msg)))
        elif level == "shellcmd":
            if self.printshellcmds:
                self.logger.warning(msg["msg"])
        elif level == "job_finished":
            # do not display this on the console for now
            pass
        elif level == "rule_info":
            self.logger.info(msg["name"])
            if msg["docstring"]:
                self.logger.info("\t" + msg["docstring"])
        elif level == "d3dag":
            print(json.dumps({"nodes": msg["nodes"], "links": msg["edges"]}))


def format_dict(dict, omit_keys=[], omit_values=[]):
    return ", ".join("{}={}".format(name, value)
                     for name, value in dict.items()
                     if name not in omit_keys and value not in omit_values)


format_resources = partial(format_dict, omit_keys={"_cores", "_nodes"})
format_wildcards = partial(format_dict, omit_values={DYNAMIC_FILL})


def format_resource_names(resources, omit_resources="_cores _nodes".split()):
    return ", ".join(name for name in resources if name not in omit_resources)


logger = Logger()


def setup_logger(handler=None,
                 quiet=False,
                 printshellcmds=False,
                 printreason=False,
                 nocolor=False,
                 stdout=False,
                 debug=False,
                 timestamp=False):
    logger.setup()
    if handler is not None:
        # custom log handler
        logger.log_handler.append(handler)
    else:
        # console output only if no custom logger was specified
        stream_handler = ColorizingStreamHandler(
            nocolor=nocolor,
            stream=sys.stdout if stdout else sys.stderr,
            timestamp=timestamp)
        logger.set_stream_handler(stream_handler)

    logger.set_level(_logging.DEBUG if debug else _logging.INFO)
    logger.quiet = quiet
    logger.printshellcmds = printshellcmds
    logger.printreason = printreason
