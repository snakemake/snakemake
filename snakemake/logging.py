__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


import logging
import logging.handlers
import platform
import time
import datetime
import sys
import os
import json
import threading
from queue import Queue
from functools import partial
from typing import TYPE_CHECKING
import textwrap
from typing import List, Optional
from snakemake_interface_logger_plugins.base import LogHandlerBase
from snakemake_interface_logger_plugins.settings import OutputSettingsLoggerInterface


if TYPE_CHECKING:
    from snakemake_interface_executor_plugins.settings import ExecMode
    from snakemake.settings.enums import Quietness


def timestamp():
    """Helper method to format the timestamp."""
    return f"[{time.asctime()}]"


def show_logs(logs):
    """Helper method to show logs."""
    for f in logs:
        try:
            with open(f, "r") as log_file:
                content = log_file.read()
        except FileNotFoundError:
            yield f"Logfile {f} not found."
            return
        except UnicodeDecodeError:
            yield f"Logfile {f} is not a text file."
            return
        lines = content.splitlines()
        logfile_header = f"Logfile {f}:"
        if not lines:
            logfile_header += " empty file"
            yield logfile_header
            return
        yield logfile_header
        max_len = min(max(max(len(line) for line in lines), len(logfile_header)), 80)
        yield "=" * max_len
        yield from lines
        yield "=" * max_len


def format_dict(dict_like, omit_keys=None, omit_values=None):
    from snakemake.io import Namedlist

    omit_keys = omit_keys or []
    omit_values = omit_values or []

    if isinstance(dict_like, (Namedlist, dict)):
        items = dict_like.items()

    else:
        raise ValueError(
            "bug: format_dict applied to something neither a dict nor a Namedlist"
        )
    return ", ".join(
        f"{name}={value}"
        for name, value in items
        if name not in omit_keys and value not in omit_values
    )


format_resources = partial(format_dict, omit_keys={"_cores", "_nodes"})
format_wildcards = format_dict


def format_resource_names(resources, omit_resources="_cores _nodes".split()):
    return ", ".join(name for name in resources if name not in omit_resources)


def format_percentage(done, total):
    """Format percentage from given fraction while avoiding superfluous precision."""
    if done == total:
        return "100%"
    if done == 0:
        return "0%"
    precision = 0
    fraction = done / total
    fmt_precision = "{{:.{}%}}".format

    def fmt(fraction):
        return fmt_precision(precision).format(fraction)

    while fmt(fraction) == "100%" or fmt(fraction) == "0%":
        precision += 1
    return fmt(fraction)


def get_level(record: logging.LogRecord) -> str:
    """
    Gets snakemake log level from a log record. If there is no snakemake log level,
    returns the log record's level name.

    Args:
        record (logging.LogRecord)
    Returns:
        str: The log level

    """
    level = record.__dict__.get("level", None)

    if level is None:
        level = record.levelname

    return level.lower()


def is_quiet_about(quiet: "Quietness", msg_type: str):
    from snakemake.settings.enums import Quietness

    parsed = Quietness.parse_choice(msg_type)

    return Quietness.ALL in quiet or parsed in quiet


class DefaultFormatter(logging.Formatter):
    def __init__(
        self,
        quiet: "Quietness",
        show_failed_logs: bool = False,
        printshellcmds: bool = False,
    ):
        self.quiet = set() if quiet is None else quiet
        self.show_failed_logs = show_failed_logs
        self.printshellcmds = printshellcmds
        self.last_msg_was_job_info = False

    def format(self, record):
        """
        Override format method to format Snakemake-specific log messages.
        """

        level = get_level(record)
        record_dict = record.__dict__.copy()

        if level == "info":
            return self.format_info(record_dict)
        elif level == "host":
            return self.format_host(record_dict)
        elif level == "job_info":
            return self.format_job_info(record_dict)
        elif level == "group_info":
            return self.format_group_info(record_dict)
        elif level == "job_error":
            return self.format_job_error(record_dict)
        elif level == "group_error":
            return self.format_group_error(record_dict)
        elif level == "progress":
            return self.format_progress(record_dict)
        elif level == "job_finished":
            return self.format_job_finished(record_dict)
        elif level == "shellcmd":
            return self.format_shellcmd(record_dict)
        elif level == "rule_info":
            return self.format_rule_info(record_dict)
        elif level == "d3dag":
            return self.format_d3dag(record_dict)
        elif level == "dag_debug":
            return self.format_dag_debug(record_dict)
        elif level == "run_info":
            return self.format_run_info(record_dict)
        else:
            return record_dict["msg"]

    def format_info(self, msg):
        """
        Format 'info' level messages.
        """
        output = []

        # Check if 'indent' is specified
        indent = "    " if msg.get("indent", False) else ""

        # Split the message by lines in case it's multiline
        lines = msg["msg"].split("\n")

        # Apply indentation to each line
        for line in lines:
            output.append(f"{indent}{line}")

        # Return the formatted message as a single string with newlines
        return "\n".join(output)

    def format_run_info(self, msg):
        """Format the run_info log messages."""
        return msg["msg"]  # Log the message directly

    def format_host(self, msg):
        """Format for host log."""
        return f"host: {platform.node()}"

    def format_job_info(self, msg):
        """Format for job_info log."""
        output = []

        output.append(timestamp())
        if msg["rule_msg"]:
            output.append(f"Job {msg['jobid']}: {msg['rule_msg']}")
            if not is_quiet_about(self.quiet, "reason"):
                output.append(f"Reason: {msg['reason']}")
        else:
            output.append("\n".join(self._format_job_info(msg)))

        if msg["is_checkpoint"]:
            output.append("DAG of jobs will be updated after completion.")
        if msg["is_handover"]:
            output.append("Handing over execution to foreign system...")

        return "\n".join(output)

    def format_group_info(self, msg):
        """Format for group_info log."""
        msg = f"{timestamp()} group job {msg['groupid']} (jobs in lexicogr. order):"

        return msg

    def format_job_error(self, msg):
        """Format for job_error log."""
        output = []
        output.append(timestamp())
        output.append("\n".join(self._format_job_error(msg)))
        return "\n".join(output)

    def format_group_error(self, msg):
        """Format for group_error log."""
        output = []
        output.append(timestamp())
        output.append("\n".join(self._format_group_error(msg)))
        return "\n".join(output)

    def format_progress(self, msg):
        """Format for progress log."""
        done = msg["done"]
        total = msg["total"]
        return f"{done} of {total} steps ({format_percentage(done, total)}) done"

    def format_job_finished(self, msg):
        """Format for job_finished log."""
        return f"{timestamp()} {msg['msg']}"

    def format_shellcmd(self, msg):
        """Format for shellcmd log."""
        if self.printshellcmds:
            return msg["msg"]
        return ""

    def format_rule_info(self, msg):
        """Format for rule_info log."""
        output = [msg["rule_name"]]
        if msg["docstring"]:
            output.append(f"    {msg['docstring']}")
        return "\n".join(output)

    def format_d3dag(self, msg):
        """Format for d3dag log."""

        return json.dumps({"nodes": msg["nodes"], "links": msg["edges"]})

    def format_dag_debug(self, msg):
        """Format for dag_debug log."""
        output = []

        if "file" in msg:
            output.append(
                f"file {msg['file']}:\n    {msg['msg']}\n{textwrap.indent(str(msg['exception']), '    ')}"
            )
        else:
            job = msg["job"]
            output.append(
                f"{msg['status']} job {job.rule.name}\n    wildcards: {format_wildcards(job.wildcards)}"
            )
        return "\n".join(output)

    def _format_job_info(self, msg):
        """Helper method to format job info details."""

        def format_item(item, omit=None, valueformat=str):
            value = msg[item]
            if value != omit:
                return f"    {item}: {valueformat(value)}"

        output = [
            f"{'local' if msg['local'] else ''}{'checkpoint' if msg['is_checkpoint'] else 'rule'} {msg['rule_name']}:"
        ]
        for item in ["input", "output", "log"]:
            fmt = format_item(item, omit=[], valueformat=", ".join)
            if fmt:
                output.append(fmt)

        singleitems = ["jobid", "benchmark"]
        if not is_quiet_about(self.quiet, "reason"):
            singleitems.append("reason")
        for item in singleitems:
            fmt = format_item(item, omit=None)
            if fmt:
                output.append(fmt)

        wildcards = format_wildcards(msg["wildcards"])
        if wildcards:
            output.append(f"    wildcards: {wildcards}")

        for item, omit in zip("priority threads".split(), [0, 1]):
            fmt = format_item(item, omit=omit)
            if fmt:
                output.append(fmt)

        resources = format_resources(msg["resources"])
        if resources:
            output.append(f"    resources: {resources}")

        return output

    def _format_job_error(self, msg):
        """Helper method to format job error details."""
        output = [f"Error in rule {msg['rule_name']}:"]

        if msg["msg"]:
            output.append(f"    message: {msg['rule_msg']}")
        output.append(f"    jobid: {msg['jobid']}")
        if msg["input"]:
            output.append(f"    input: {', '.join(msg['input'])}")
        if msg["output"]:
            output.append(f"    output: {', '.join(msg['output'])}")
        if msg["log"]:
            output.append(
                f"    log: {', '.join(msg['log'])} (check log file(s) for error details)"
            )
        if msg["conda_env"]:
            output.append(f"    conda-env: {msg['conda_env']}")
        if msg["shellcmd"]:
            output.append(
                f"    shell:\n        {msg['shellcmd']}\n        (command exited with non-zero exit code)"
            )

        for item in msg["aux"].items():
            output.append(f"    {item[0]}: {item[1]}")

        if self.show_failed_logs and msg["log"]:
            output.extend(show_logs(msg["log"]))

        return output

    def _format_group_error(self, msg):
        """Helper method to format group error details."""
        output = [f"Error in group {msg['groupid']}:"]

        if msg["msg"]:
            output.append(f"    message: {msg['msg']}")
        if msg["aux_logs"]:
            output.append(
                f"    log: {', '.join(msg['aux_logs'])} (check log file(s) for error details)"
            )

        output.append("    jobs:")
        for info in msg["job_error_info"]:
            output.append(f"        rule {info['name']}:")
            output.append(f"            jobid: {info['jobid']}")
            if info["output"]:
                output.append(f"            output: {', '.join(info['output'])}")
            if info["log"]:
                output.append(
                    f"            log: {', '.join(info['log'])} (check log file(s) for error details)"
                )

        logs = msg["aux_logs"] + [
            f for info in msg["job_error_info"] for f in info["log"]
        ]
        if self.show_failed_logs and logs:
            output.extend(self.show_logs(logs))

        return output


class DefaultFilter:
    def __init__(self, quiet, debug_dag, dryrun) -> None:
        if quiet is None:
            quiet = set()
        self.quiet = quiet
        self.debug_dag = debug_dag
        self.dryrun = dryrun

    def filter(self, record):
        from snakemake.settings.enums import Quietness

        level = get_level(record)
        if self.dryrun and level == "run_info":
            return True

        if Quietness.ALL in self.quiet and not self.dryrun:
            return False

        quietness_map = {
            "job_info": Quietness.RULES,
            "group_info": Quietness.RULES,
            "job_error": Quietness.RULES,
            "group_error": Quietness.RULES,
            "progress": Quietness.PROGRESS,
            "shellcmd": Quietness.PROGRESS,
            "job_finished": Quietness.PROGRESS,
            "resources_info": Quietness.PROGRESS,
            "run_info": Quietness.PROGRESS,
            "host": Quietness.HOST,
            "info": Quietness.PROGRESS,
        }

        # Check quietness for specific levels
        if level in quietness_map:
            if quietness_map[level] in self.quiet:
                return False

        # Handle dag_debug specifically
        if level == "dag_debug" and not self.debug_dag:
            return False

        return True


class ColorizingTextHandler(logging.StreamHandler):
    """
    Custom handler that combines colorization and Snakemake-specific formatting.
    """

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[%dm"
    BOLD_SEQ = "\033[1m"

    colors = {
        "WARNING": YELLOW,
        "INFO": GREEN,
        "DEBUG": BLUE,
        "CRITICAL": MAGENTA,
        "ERROR": RED,
    }

    def __init__(
        self,
        nocolor=False,
        stream=sys.stderr,
        mode=None,
        formatter: Optional[logging.Formatter] = None,
        filter: Optional[logging.Filter] = None,
    ):
        super().__init__(stream=stream)
        self.last_msg_was_job_info = False
        self._output_lock = threading.Lock()
        self.nocolor = nocolor or not self.can_color_tty(mode)
        self.mode = mode

        if formatter:
            self.setFormatter(formatter)
        if filter:
            self.addFilter(filter)

    def can_color_tty(self, mode):
        """
        Colors are supported when:
        1. Terminal is not "dumb"
        2. Running in subprocess mode
        3. Using a TTY on non-Windows systems
        """
        from snakemake_interface_executor_plugins.settings import ExecMode

        # Case 1: Check if terminal is "dumb"
        if "TERM" in os.environ:
            if os.environ["TERM"] == "dumb":
                return False

        # Case 2: Always support colors in subprocess mode
        if mode == ExecMode.SUBPROCESS:
            return True

        # Case 3: Support colors on TTY except for Windows
        is_windows = platform.system() == "Windows"
        has_tty = self.is_tty

        if has_tty and not is_windows:
            return True

        return False

    @property
    def is_tty(self):
        isatty = getattr(self.stream, "isatty", None)
        return isatty and isatty()

    def emit(self, record):
        """
        Emit a log message with custom formatting and color.
        """

        with self._output_lock:
            try:
                level = get_level(record)

                if level == "job_info":
                    if not self.last_msg_was_job_info:
                        self.stream.write(
                            "\n"
                        )  # Add a blank line before a new job_info message
                    self.last_msg_was_job_info = True
                else:
                    # Reset flag if the message is not a 'job_info'
                    self.last_msg_was_job_info = False
                formatted_message = self.format(record)

                # Apply color to the formatted message
                self.stream.write(self.decorate(record, formatted_message))
                self.stream.write(getattr(self, "terminator", "\n"))
                self.flush()
            except BrokenPipeError:
                raise
            except (KeyboardInterrupt, SystemExit):
                pass  # Ignore exceptions for these cases, all errors have been handled before.
            except Exception:
                self.handleError(record)

    def decorate(self, record, message):
        """
        Add color to the log message based on its level.
        """
        message = [message]

        if record.levelname == "INFO" and not hasattr(record, "level"):
            record.levelname = "WARNING"  # mimic old snakemake logging color

        if not self.nocolor and record.levelname in self.colors:
            message.insert(0, self.COLOR_SEQ % (30 + self.colors[record.levelname]))
            message.append(self.RESET_SEQ)

        return "".join(message)


class LoggerManager:
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.initialized = False
        self.queue_listener = None
        self.mode = None

        self.logfile_handlers = {}
        self.settings: OutputSettingsLoggerInterface = None

    def setup(
        self,
        mode: "ExecMode",
        handlers: List[LogHandlerBase],
        settings: OutputSettingsLoggerInterface,
    ):
        from snakemake_interface_executor_plugins.settings import ExecMode

        self.mode = mode
        self.settings = settings
        self.initialized = True

        stream_handlers = []
        other_handlers = []

        if self.mode == ExecMode.SUBPROCESS:
            handler = self._default_streamhandler()
            handler.setLevel(logging.ERROR)
            stream_handlers.append(handler)
        elif self.mode == ExecMode.REMOTE:
            stream_handlers.append(self._default_streamhandler())
        elif handlers:
            for handler in handlers:
                configured_handler = self._configure_plugin_handler(handler)
                if configured_handler.writes_to_file:
                    self.logfile_handlers[configured_handler] = (
                        configured_handler.baseFilename
                    )
                elif configured_handler.writes_to_stream:
                    stream_handlers.append(configured_handler)
                else:
                    other_handlers.append(configured_handler)
        else:
            stream_handlers.append(self._default_streamhandler())
        if len(stream_handlers) > 1:
            raise ValueError("More than 1 stream log handler specified!")
        all_handlers = (
            stream_handlers + other_handlers + list(self.logfile_handlers.keys())
        )

        q = Queue(-1)
        self.queue_listener = logging.handlers.QueueListener(
            q,
            *all_handlers,
            respect_handler_level=True,
        )
        self.queue_listener.start()
        self.logger.setLevel(logging.DEBUG if settings.verbose else logging.INFO)
        self.logger.addHandler(logging.handlers.QueueHandler(q))

    def _configure_plugin_handler(self, plugin):
        if not plugin.has_filter:
            plugin.addFilter(self._default_filter())
        if not plugin.has_formatter:
            plugin.setFormatter(self._default_formatter())
        return plugin

    def _default_filter(self):
        return DefaultFilter(
            self.settings.quiet, self.settings.debug_dag, self.settings.dryrun
        )

    def _default_formatter(self):
        return DefaultFormatter(
            self.settings.quiet,
            self.settings.show_failed_logs,
            self.settings.printshellcmds,
        )

    def _default_filehandler(self, logfile):
        logfile_handler = logging.FileHandler(logfile)
        logfile_handler.setFormatter(self._default_formatter())
        logfile_handler.addFilter(self._default_filter())
        logfile_handler.setLevel(
            logging.DEBUG if self.settings.verbose else logging.INFO
        )
        logfile_handler.name = "DefaultLogFileHandler"
        return logfile_handler

    def _default_streamhandler(self):
        stream_handler = ColorizingTextHandler(
            nocolor=self.settings.nocolor,
            stream=sys.stdout if self.settings.stdout else sys.stderr,
            mode=self.mode,
        )
        stream_handler.addFilter(self._default_filter())
        stream_handler.setFormatter(self._default_formatter())
        stream_handler.name = "DefaultStreamHandler"
        return stream_handler

    def get_logfile(self):
        return self.logfile_handlers.values()

    def logfile_hint(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        """Log the logfile location if applicable."""
        logfiles = self.logfile_handlers.values()
        if self.mode == ExecMode.DEFAULT and not self.settings.dryrun and logfiles:
            log_paths = ", ".join([os.path.abspath(p) for p in logfiles])
            self.logger.info(f"Complete log(s): {log_paths}")
        return logfiles

    def cleanup_logfile(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT:
            for handler in self.logfile_handlers.keys():
                handler.close()

    def setup_logfile(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT and not self.settings.dryrun:
            os.makedirs(os.path.join(".snakemake", "log"), exist_ok=True)
            logfile = os.path.abspath(
                os.path.join(
                    ".snakemake",
                    "log",
                    datetime.datetime.now().isoformat().replace(":", "")
                    + ".snakemake.log",
                )
            )
            handler = self._default_filehandler(logfile)
            self.logfile_handlers[handler] = logfile

    def stop(self):
        if self.queue_listener is not None and self.queue_listener._thread is not None:
            self.queue_listener.stop()


# Global logger instance
logger = logging.getLogger(__name__)
logger_manager = LoggerManager(logger)
