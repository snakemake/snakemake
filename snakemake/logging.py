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

import textwrap
from typing import Optional
from snakemake_interface_logger_plugins.base import LoggerPluginBase


def get_default_exec_mode():
    from snakemake_interface_executor_plugins.settings import ExecMode

    return ExecMode.DEFAULT


# Helper functions for logging


def timestamp():
    """Helper method to format the timestamp."""
    return f"[{time.asctime()}]"


def show_logs(logs):
    """Helper method to show logs."""
    for f in logs:
        try:
            content = open(f, "r").read()
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
        max_len = min(max(max(len(l) for l in lines), len(logfile_header)), 80)
        yield "=" * max_len
        yield from lines
        yield "=" * max_len


def format_dict(dict_like, omit_keys=None, omit_values=None):
    from snakemake.io import Namedlist

    omit_keys = omit_keys or []
    omit_values = omit_values or []

    if isinstance(dict_like, Namedlist):
        items = dict_like.items()
    elif isinstance(dict_like, dict):
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
    fmt = lambda fraction: fmt_precision(precision).format(fraction)
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


class DefaultFormatter(logging.Formatter):
    def __init__(self, printreason=False, show_failed_logs=False, printshellcmds=False):
        self.printreason = printreason
        self.show_failed_logs = show_failed_logs
        self.printshellcmds = printshellcmds
        self.last_msg_was_job_info = False

    def format(self, record):
        """
        Override format method to format Snakemake-specific log messages.
        """

        level = get_level(record)
        record_dict = record.__dict__.copy()

        # Call specific formatters based on the log level
        match level:
            case "info":
                return self.format_info(record_dict)
            case "host":
                return self.format_host(record_dict)
            case "job_info":
                return self.format_job_info(record_dict)
            case "group_info":
                return self.format_group_info(record_dict)
            case "job_error":
                return self.format_job_error(record_dict)
            case "group_error":
                return self.format_group_error(record_dict)
            case "progress":
                return self.format_progress(record_dict)
            case "job_finished":
                return self.format_job_finished(record_dict)
            case "shellcmd":
                return self.format_shellcmd(record_dict)
            case "rule_info":
                return self.format_rule_info(record_dict)
            case "d3dag":
                return self.format_d3dag(record_dict)
            case "dag_debug":
                return self.format_dag_debug(record_dict)
            case "run_info":
                return self.format_run_info(record_dict)
            case _:
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
            if self.printreason:
                output.append(f"Reason: {msg['reason']}")
        else:
            output.append("\n".join(self._format_job_info(msg)))

        if msg["is_checkpoint"]:
            output.append("DAG of jobs will be updated after completion.")
        if msg["is_handover"]:
            output.append("Handing over execution to foreign system...")
        output.append("")
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
        return f"{timestamp()} Finished job {msg['jobid']}."

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
        if self.printreason:
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
    def __init__(
        self,
        quiet,
        debug_dag,
    ) -> None:
        if quiet is None:
            quiet = set()
        self.quiet = quiet
        self.debug_dag = debug_dag

    def is_quiet_about(self, msg_type: str):
        from snakemake.settings.enums import Quietness

        return (
            Quietness.ALL in self.quiet
            or Quietness.parse_choice(msg_type) in self.quiet
        )

    def filter(self, record):
        if self.is_quiet_about("all"):
            return False

        level = get_level(record)

        # Respect quiet mode filtering
        quietness_map = {
            "job_info": "rules",
            "group_info": "rules",
            "job_error": "rules",
            "group_error": "rules",
            "progress": "progress",
            "shellcmd": "progress",  # or other quietness types
            "job_finished": "progress",
            "resources_info": "progress",
            "run_info": "progress",
            "host": "host",
            "info": "progress",
        }

        # Check quietness for specific levels and skip accordingly
        if level in quietness_map:
            if self.is_quiet_about(quietness_map[level]):
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
        Check if the terminal supports colors.
        """
        from snakemake_interface_executor_plugins.settings import ExecMode

        if "TERM" in os.environ and os.environ["TERM"] == "dumb":
            return False
        if mode == ExecMode.SUBPROCESS:
            return True
        return self.is_tty and not platform.system() == "Windows"

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
            except BrokenPipeError as e:
                raise e
            except (KeyboardInterrupt, SystemExit):
                pass  # Ignore exceptions for these cases
            except Exception as e:
                self.handleError(record)

    def is_quiet_about(self, msg_type: str):
        from snakemake.settings.enums import Quietness

        return (
            Quietness.ALL in self.quiet
            or Quietness.parse_choice(msg_type) in self.quiet
        )

    def decorate(self, record, message):
        """
        Add color to the log message based on its level.
        """
        message = [message]  # Treat the formatted message as a list

        # Use record.levelname to apply color if applicable
        if not self.nocolor and record.levelname in self.colors:
            message.insert(0, self.COLOR_SEQ % (30 + self.colors[record.levelname]))
            message.append(self.RESET_SEQ)

        return "".join(message)


class LoggerManager:
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.quiet = None
        self.printshellcmds = False
        self.printreason = True
        self.debug_dag = False
        self.nocolor = False
        self.stdout = False
        self.debug = False
        self.mode = None
        self.show_failed_logs = False
        self.dryrun = False
        self.handlers = []
        self.initialized = False
        self.queue_listener = None
        self.formatter = DefaultFormatter(
            printreason=self.printreason,
            show_failed_logs=self.show_failed_logs,
            printshellcmds=self.printshellcmds,
        )
        self.filter = DefaultFilter(quiet=self.quiet, debug_dag=self.debug_dag)

    def _get_handlers_of_type(self, handler_type: type):
        """Helper function to get all handlers of a specified type."""
        return [
            h
            for h in self.handlers
            if isinstance(h, handler_type) or issubclass(type(h), handler_type)
        ]

    def get_logfile(self) -> Optional[list[str]]:
        """Retrieve the log file paths from the file handlers in the logger."""
        logfiles = [
            h.baseFilename for h in self._get_handlers_of_type(logging.FileHandler)
        ]
        return logfiles if logfiles else None

    def logfile_hint(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        """Log the logfile location if applicable."""
        logfiles = self.get_logfile()
        if self.mode == ExecMode.DEFAULT and not self.dryrun and logfiles:
            log_paths = ", ".join([os.path.relpath(p) for p in logfiles])
            self.logger.info(f"Complete log(s): {log_paths}")
        return logfiles

    def cleanup_logfile(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT:
            for handler in self._get_handlers_of_type(logging.FileHandler):
                handler.close()

    def setup_logfile(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT and not self.dryrun:
            os.makedirs(os.path.join(".snakemake", "log"), exist_ok=True)
            logfile = os.path.abspath(
                os.path.join(
                    ".snakemake",
                    "log",
                    datetime.datetime.now().isoformat().replace(":", "")
                    + ".snakemake.log",
                )
            )

            logfile_handler = logging.FileHandler(logfile)
            logfile_handler.setFormatter(self.formatter)
            logfile_handler.addFilter(self.filter)
            logfile_handler.setLevel(logging.DEBUG if self.debug else logging.INFO)
            logfile_handler.name = "DefaultLogFileHandler"
            self.handlers.append(logfile_handler)

    def add_stream_handler(
        self,
        stream_handler: Optional[logging.StreamHandler] = None,
        use_default_filter: bool = True,
        use_default_formatter: bool = True,
    ):
        """
        Adds a stream handler. Checks if one already exists by the same name before adding.
        """
        if stream_handler is None:
            stream_handler = ColorizingTextHandler(
                nocolor=self.nocolor,
                stream=sys.stdout if self.stdout else sys.stderr,
                mode=self.mode,
            )
        if not stream_handler.filters and use_default_filter:
            stream_handler.addFilter(self.filter)
        if not stream_handler.formatter and use_default_formatter:
            stream_handler.setFormatter(self.formatter)

        if not stream_handler.name:
            stream_handler.name = "DefaultStreamHandler"

        # try to prevent multiple stream handlers. may need to revisit this.
        if not any([getattr(h, "stream", False) for h in self.handlers]):
            self.handlers.append(stream_handler)

    def stop(self):
        if self.queue_listener is not None:
            self.queue_listener.stop()

    def configure_logger(
        self,
        quiet=None,
        printshellcmds: Optional[bool] = None,
        printreason: Optional[bool] = None,
        debug_dag: Optional[bool] = None,
        nocolor: Optional[bool] = None,
        stdout: Optional[bool] = None,
        debug: Optional[bool] = None,
        mode=None,
        show_failed_logs: Optional[bool] = None,
        dryrun: Optional[bool] = None,
        plugins: list[LoggerPluginBase] = None,
    ):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.initialized:
            return

        for key, value in locals().items():
            if key != "self" and value is not None:
                setattr(self, key, value)

        # Update the logger settings based on the current mode
        # Subproecces execMode gets nullhandler so nothing gets logged
        # TODO: Think about if remote should be null handler too. Don't want remote logging to HTTP logger for example.
        if self.mode == ExecMode.SUBPROCESS:
            self.add_stream_handler(logging.NullHandler())
        elif self.mode == ExecMode.REMOTE:
            self.add_stream_handler()
        elif plugins:
            for plugin in plugins:
                handler = plugin.create_handler(
                    quiet=quiet,
                    printreason=printreason,
                    printshellcmds=printshellcmds,
                    debug=debug,
                    debug_dag=debug_dag,
                    mode=mode,
                    show_failed_logs=show_failed_logs,
                    dryrun=dryrun,
                    nocolor=nocolor,
                    stdout=stdout,
                )
                # logging.Streamhandlers (and subclasses of it) should have attr called stream
                if getattr(handler, "stream", False):
                    self.add_stream_handler(
                        stream_handler=handler,
                        use_default_filter=False,
                        use_default_formatter=False,
                    )
                else:
                    self.handlers.append(handler)
        else:
            self.add_stream_handler()

        q = Queue(-1)
        self.queue_listener = logging.handlers.QueueListener(
            q,
            *self.handlers,
            respect_handler_level=True,
        )
        self.queue_listener.start()
        self.logger.addHandler(logging.handlers.QueueHandler(q))
        self.initialized = True

        self.logger.setLevel(logging.DEBUG if self.debug else logging.INFO)


# Global logger instance
logger = logging.getLogger(__name__)
logger_manager = LoggerManager(logger)
