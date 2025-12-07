from __future__ import annotations

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
import threading
from queue import Queue
from functools import partial
from typing import TYPE_CHECKING, Any
import textwrap
from typing import List, Optional, Collection, TextIO
from snakemake_interface_logger_plugins.base import LogHandlerBase
from snakemake_interface_logger_plugins.common import LogEvent

if TYPE_CHECKING:
    from snakemake.settings.enums import Quietness
    from snakemake.settings.types import OutputSettings


def timestamp() -> str:
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


def format_dict(dict_like, omit_keys=None, omit_values=None) -> str:
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


def format_percentage(done: int, total: int) -> str:
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


def get_event(record: logging.LogRecord) -> Optional[LogEvent]:
    """Get snakemake log event from a log record."""
    return getattr(record, "event", None)


def is_quiet_about(quiet: Collection["Quietness"], msg_type: str) -> bool:
    from snakemake.settings.enums import Quietness

    parsed = Quietness.parse_choice(msg_type)

    return Quietness.ALL in quiet or parsed in quiet


class DefaultFormatter(logging.Formatter):
    quiet: Collection["Quietness"]
    show_failed_logs: bool

    def __init__(
        self,
        quiet: Optional[Collection["Quietness"]],
        show_failed_logs: bool = False,
    ):
        self.quiet = set() if quiet is None else quiet
        self.show_failed_logs = show_failed_logs
        self.last_msg_was_job_info = False

    def format(self, record: logging.LogRecord) -> str:
        """
        Override format method to format Snakemake-specific log messages.
        """
        event = get_event(record)
        record_dict = record.__dict__.copy()

        def default_formatter(rd):
            return rd["msg"] or ""  # handle None

        formatters = {
            None: default_formatter,
            LogEvent.JOB_INFO: self.format_job_info,
            LogEvent.JOB_ERROR: self.format_job_error,
            LogEvent.JOB_FINISHED: self.format_job_finished,
            LogEvent.GROUP_INFO: self.format_group_info,
            LogEvent.GROUP_ERROR: self.format_group_error,
            LogEvent.SHELLCMD: self.format_shellcmd,
            LogEvent.RUN_INFO: self.format_run_info,
            LogEvent.DEBUG_DAG: self.format_dag_debug,
            LogEvent.PROGRESS: self.format_progress,
        }

        formatter = formatters.get(event, default_formatter)
        return formatter(record_dict)

    def format_run_info(self, msg: dict[str, Any]):
        """Format the run_info log messages."""
        return msg["msg"]  # Log the message directly

    def format_host(self, msg: dict[str, Any]):
        """Format for host log."""
        return f"host: {platform.node()}"

    def format_job_info(self, msg: dict[str, Any]):
        """Format for job_info log."""
        output = []

        output.append(timestamp())
        if msg["rule_msg"]:
            output.append(f"Job {msg['jobid']}: {msg['rule_msg']}")
            if not is_quiet_about(self.quiet, "reason"):
                output.append(f"Reason: {msg['reason']}")
        else:
            output.append("\n".join(self._format_job_info(msg)))

        if msg.get("indent", False):
            return textwrap.indent("\n".join(output), "    ")
        return "\n".join(output)

    def format_group_info(self, msg: dict[str, Any]):
        """Format for group_info log."""
        msg = f"{timestamp()} {msg['msg']}"

        return msg

    def format_job_error(self, msg: dict[str, Any]):
        """Format for job_error log."""
        output = []
        output.append(timestamp())
        output.append("\n".join(self._format_job_error(msg)))

        if msg.get("indent", False):
            return textwrap.indent("\n".join(output), "    ")
        return "\n".join(output)

    def format_group_error(self, msg: dict[str, Any]):
        """Format for group_error log."""
        output = []
        output.append(timestamp())
        output.append("\n".join(self._format_group_error(msg)))
        return "\n".join(output)

    def format_progress(self, msg: dict[str, Any]):
        """Format for progress log."""
        done = msg["done"]
        total = msg["total"]
        return f"{done} of {total} steps ({format_percentage(done, total)}) done"

    def format_job_finished(self, msg: dict[str, Any]):
        """Format for job_finished log."""
        return f"{timestamp()}\n{msg['msg']}"

    def format_shellcmd(self, msg: dict[str, Any]):
        """Format for shellcmd log."""
        return msg["msg"]

    def format_dag_debug(self, msg: dict[str, Any]):
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

    def _format_job_info(self, msg: dict[str, Any]):
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

    def _format_job_error(self, msg: dict[str, Any]):
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

    def _format_group_error(self, msg: dict[str, Any]):
        """Helper method to format group error details."""
        output = []

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
            output.extend(show_logs(logs))

        return output


class DefaultFilter:
    """Default log filter.

    Attributes
    ----------
    quiet
        Quietness values to filter out.
    debug_dag
        Whether to allow DEBUG_DAG events.
    dryrun
    printshellcmds
        Whether to allow SHELLCMD events.
    """

    quiet: Collection["Quietness"]
    debug_dag: bool
    dryrun: bool
    printshellcmds: bool

    def __init__(
        self,
        quiet: Optional[Collection["Quietness"]],
        debug_dag: bool,
        dryrun: bool,
        printshellcmds: bool,
    ) -> None:
        if quiet is None:
            quiet = set()
        self.quiet = quiet
        self.debug_dag = debug_dag
        self.dryrun = dryrun
        self.printshellcmds = printshellcmds

    def filter(self, record: logging.LogRecord) -> bool:
        from snakemake.settings.enums import Quietness

        event = get_event(record)

        if Quietness.ALL in self.quiet and not self.dryrun:
            return False

        if hasattr(record, "quietness"):
            if record.quietness in self.quiet:
                return False

        quietness_map = {
            LogEvent.JOB_INFO: Quietness.RULES,
            LogEvent.GROUP_INFO: Quietness.RULES,
            LogEvent.JOB_ERROR: Quietness.RULES,
            LogEvent.GROUP_ERROR: Quietness.RULES,
            LogEvent.PROGRESS: Quietness.PROGRESS,
            LogEvent.SHELLCMD: Quietness.RULES,
            LogEvent.JOB_FINISHED: Quietness.PROGRESS,
            LogEvent.RESOURCES_INFO: Quietness.PROGRESS,
            LogEvent.RUN_INFO: Quietness.PROGRESS,
        }

        # Handle shell commands
        if event == LogEvent.SHELLCMD and not self.printshellcmds:
            return False

        # Check quietness for specific levels
        if event in quietness_map:
            if quietness_map[event] in self.quiet:
                return False

        # Handle dag_debug specifically
        if event == LogEvent.DEBUG_DAG and not self.debug_dag:
            return False
        if event == LogEvent.WORKFLOW_STARTED:
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

    yellow_info_events = [
        LogEvent.RUN_INFO,
        LogEvent.SHELLCMD,
        LogEvent.JOB_STARTED,
        None,  # To mimic old coloring where log.info was mapped to log.warn
    ]

    nocolor: bool
    stream: TextIO

    def __init__(
        self,
        nocolor: bool = False,
        stream: TextIO = sys.stderr,
        formatter: Optional[logging.Formatter] = None,
        filter: Optional[logging.Filter] = None,
    ):
        super().__init__(stream=stream)
        self.last_msg_was_job_info = False
        self._output_lock = threading.Lock()
        self.nocolor = nocolor or not self.can_color_tty()

        if formatter:
            self.setFormatter(formatter)
        if filter:
            self.addFilter(filter)

    def can_color_tty(self) -> bool:
        """
        Colors are supported when:
        1. Terminal is not "dumb"
        2. Using a TTY on non-Windows systems
        """
        # Case 1: Check if terminal is "dumb"
        if os.environ.get("TERM") == "dumb":
            return False

        # Case 2: Support colors on TTY except for Windows
        is_windows = platform.system() == "Windows"
        has_tty = self.is_tty

        if has_tty and not is_windows:
            return True

        return False

    @property
    def is_tty(self) -> bool:
        isatty = getattr(self.stream, "isatty", None)
        return bool(isatty and isatty())

    def emit(self, record: logging.LogRecord) -> None:
        """
        Emit a log message with custom formatting and color.
        """

        with self._output_lock:
            try:
                event = get_event(record)

                if event == LogEvent.JOB_INFO:
                    if not self.last_msg_was_job_info:
                        self.stream.write(
                            "\n"
                        )  # Add a blank line before a new job_info message
                    self.last_msg_was_job_info = True
                else:
                    # Reset flag if the message is not a 'job_info'
                    self.last_msg_was_job_info = False
                formatted_message = self.format(record)
                if formatted_message == "None" or formatted_message == "":
                    return
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

    def decorate(self, record: logging.LogRecord, message: str) -> str:
        """
        Add color to the log message based on its level.
        """
        message = [message]

        event = get_event(record)

        if not self.nocolor and record.levelname in self.colors:
            if record.levelno == logging.INFO and event in self.yellow_info_events:
                color = self.colors["WARNING"]
            else:
                color = self.colors[record.levelname]

            message.insert(0, self.COLOR_SEQ % (30 + color))
            message.append(self.RESET_SEQ)

        return "".join(message)


class LoggerManager:
    """Sets up and manages workflow logging system.

    Attributes
    ----------
    logger
        Logger object all handlers are attached to.
    initialized
        Whether :meth:`setup` has been called.
    queue_listener
        Queue listener used to process all plugin handlers in the main thread. An associated
        :class:`logging.handlers.QueueHandler` is attached to :attr:`logger`.
    needs_rulegraph
        Whether any plugin requested a RULEGRAPH event be logged.
    logfile_handlers
        Mapping from :class:`logging.Handler` instances to their associated output files. Used to
        report all log files at the end of the run.
    settings
        Global logging settings. This is used to configure the default stream/file handlers and is
        also passed to all plugins.
    """

    logger: logging.Logger
    initialized: bool
    queue_listener: Optional[logging.handlers.QueueListener]
    needs_rulegraph: bool
    logfile_handlers: dict[logging.Handler, str]
    settings: Optional["OutputSettings"]

    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.queue_listener = None

        self.needs_rulegraph = False
        self.logfile_handlers = {}
        self.settings = None

    def setup(
        self,
        handlers: List[LogHandlerBase],
        settings: "OutputSettings",
    ) -> None:
        """Set up the logging system based on settings and handlers."""
        # Clear any existing handlers to prevent duplicates
        self.logger.handlers.clear()
        self.logfile_handlers.clear()
        self.settings = settings

        # Skip plugin handlers if requested and return early
        if settings.skip_plugin_handlers or not handlers:
            handler = self._default_streamhandler()
            if settings.log_level_override is not None:
                handler.setLevel(settings.log_level_override)
            self.logger.addHandler(handler)

            # Set logger level
            if settings.log_level_override is not None:
                self.logger.setLevel(settings.log_level_override)
            else:
                self.logger.setLevel(
                    logging.DEBUG if settings.verbose else logging.INFO
                )
            return

        # Configure plugin handlers
        plugin_handlers = []
        has_stream_handler = False

        for handler in handlers:
            if handler.needs_rulegraph:
                self.needs_rulegraph = True

            configured_handler = self._configure_plugin_handler(handler)

            # Check for stream handlers (only one allowed)
            if configured_handler.writes_to_stream:
                if has_stream_handler:
                    raise ValueError("More than 1 stream log handler specified!")
                has_stream_handler = True

            # Track file handlers for later reference
            if configured_handler.writes_to_file:
                self.logfile_handlers[configured_handler] = (
                    configured_handler.baseFilename
                )

            plugin_handlers.append(configured_handler)

        # Set up queue for thread-safe plugin handlers
        if plugin_handlers:
            self._queue = Queue(-1)
            self.queue_listener = logging.handlers.QueueListener(
                self._queue,
                *plugin_handlers,
                respect_handler_level=True,
            )
            self.queue_listener.start()
            self.logger.addHandler(logging.handlers.QueueHandler(self._queue))

        # Ensure we have console output
        if not has_stream_handler:
            self.logger.addHandler(self._default_streamhandler())

        # Set final logger level
        if settings.log_level_override is not None:
            self.logger.setLevel(settings.log_level_override)
        else:
            self.logger.setLevel(logging.DEBUG if settings.verbose else logging.INFO)

        self.initialized = True

    def _configure_plugin_handler(self, plugin: LogHandlerBase) -> LogHandlerBase:
        if not plugin.has_filter:
            plugin.addFilter(self._default_filter())
        if not plugin.has_formatter:
            plugin.setFormatter(self._default_formatter())
        return plugin

    def _default_filter(self) -> DefaultFilter:
        return DefaultFilter(
            self.settings.quiet,
            self.settings.debug_dag,
            self.settings.dryrun,
            self.settings.printshellcmds,
        )

    def _default_formatter(self) -> DefaultFormatter:
        return DefaultFormatter(
            self.settings.quiet,
            self.settings.show_failed_logs,
        )

    def _default_filehandler(self, logfile) -> logging.Handler:
        logfile_handler = logging.FileHandler(logfile)
        logfile_handler.setFormatter(self._default_formatter())
        logfile_handler.addFilter(self._default_filter())
        logfile_handler.setLevel(
            logging.DEBUG if self.settings.verbose else logging.INFO
        )
        logfile_handler.name = "DefaultLogFileHandler"
        return logfile_handler

    def _default_streamhandler(self) -> logging.Handler:
        stream_handler = ColorizingTextHandler(
            nocolor=self.settings.nocolor,
            stream=sys.stdout if self.settings.stdout else sys.stderr,
        )
        stream_handler.addFilter(self._default_filter())
        stream_handler.setFormatter(self._default_formatter())
        stream_handler.name = "DefaultStreamHandler"
        return stream_handler

    def get_logfile(self) -> List[str]:
        return list(self.logfile_handlers.values())

    def logfile_hint(self):
        """Log the logfile location if applicable."""
        logfiles = self.logfile_handlers.values()
        if self.settings.enable_file_logging and not self.settings.dryrun and logfiles:
            log_paths = ", ".join([os.path.abspath(p) for p in logfiles])
            self.logger.info(f"Complete log(s): {log_paths}")
        return logfiles

    def cleanup_logfile(self) -> None:
        if self.settings.enable_file_logging:
            for handler in self.logfile_handlers.keys():
                self.logger.removeHandler(handler)
                handler.close()
            self.logfile_handlers.clear()

    def setup_logfile(self, workdir: Optional[os.PathLike] = None) -> None:
        if self.settings.enable_file_logging and not self.settings.dryrun:
            if workdir:
                logdir = os.path.join(workdir, ".snakemake", "log")
            else:
                logdir = os.path.join(".snakemake", "log")
            try:
                os.makedirs(logdir, exist_ok=True)
                logfile = os.path.abspath(
                    os.path.join(
                        logdir,
                        datetime.datetime.now().isoformat().replace(":", "")
                        + ".snakemake.log",
                    )
                )
                handler = self._default_filehandler(logfile)
                self.logger.addHandler(handler)
                self.logfile_handlers[handler] = logfile

            except OSError as e:
                self.logger.error(f"Failed to setup log file: {e}")

    def stop(self) -> None:
        if self.queue_listener is not None and self.queue_listener._thread is not None:
            self.queue_listener.stop()


# Global logger instance
logger = logging.getLogger(__name__)
logger_manager = LoggerManager(logger)
