__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import logging as _logging
import platform
import time
import datetime
import sys
import os
import json
import threading
from functools import partial
import inspect
import textwrap


def get_default_exec_mode():
    from snakemake_interface_executor_plugins.settings import ExecMode

    return ExecMode.DEFAULT


class SlackLogger:
    def __init__(self):
        from slack_sdk import WebClient

        self.token = os.getenv("SLACK_TOKEN")
        if not self.token:
            print(
                "The use of slack logging requires the user to set a user specific slack User OAuth token to the SLACK_TOKEN environment variable. Set this variable by 'export SLACK_TOKEN=your_token'. To generate your token please visit https://api.slack.com/authentication/token-types#user."
            )
            exit(-1)
        self.slack = WebClient(self.token)
        # Check for success
        try:
            auth = self.slack.auth_test().data
        except Exception:
            print(
                "Slack connection failed. Please compare your provided slack token exported in the SLACK_TOKEN environment variable with your online token (app). This token can be tested at https://api.slack.com/methods/auth.test/test. A different token can be set up by 'export SLACK_TOKEN=your_token'."
            )
            exit(-1)
        self.own_id = auth["user_id"]
        self.error_occured = False
        self.slack.chat_postMessage(
            channel=self.own_id, text="Snakemake has connected."
        )

    def log_handler(self, msg):
        if "error" in msg["level"] and not self.error_occured:
            self.slack.chat_postMessage(
                channel=self.own_id, text="At least one error occurred."
            )
            self.error_occured = True

        if msg["level"] == "progress" and msg["done"] == msg["total"]:
            # workflow finished
            self.slack.chat_postMessage(channel=self.own_id, text="Workflow complete.")


class WMSLogger:
    def __init__(self, address=None, args=None, metadata=None):
        """A WMS monitor is a workflow management system logger to enable
        monitoring with something like Panoptes. The address corresponds to
        the --wms-monitor argument, and args should be a list of key/value
        pairs with extra arguments to send to identify the workflow. We require
        the logging server to exist and receive creating a workflow to start
        the run, but we don't exit with error if any updates fail, as the
        workflow will already be running and it would not be worth stopping it.
        """

        from snakemake.resources import DefaultResources

        self.address = address or "http:127.0.0.1:5000"
        self.args = list(map(DefaultResources.decode_arg, args)) if args else []
        self.args = {item[0]: item[1] for item in list(self.args)}

        self.metadata = metadata or {}

        # A token is suggested but not required, depends on server
        self.token = os.getenv("WMS_MONITOR_TOKEN")
        self.service_info()

        # Create or retrieve the existing workflow
        self.create_workflow()

    def service_info(self):
        """Service Info ensures that the server is running. We exit on error
        if this isn't the case, so the function can be called in init.
        """
        import requests

        # We first ensure that the server is running, period
        response = requests.get(
            f"{self.address}/api/service-info", headers=self._headers
        )
        if response.status_code != 200:
            sys.stderr.write(f"Problem with server: {self.address} {os.linesep}")
            sys.exit(-1)

        # And then that it's ready to be interacted with
        if response.json().get("status") != "running":
            sys.stderr.write(
                f"The status of the server {self.address} is not in 'running' mode {os.linesep}"
            )
            sys.exit(-1)

    def create_workflow(self):
        """Creating a workflow means pinging the wms server for a new id, or
        if providing an argument for an existing workflow, ensuring that
        it exists and receiving back the same identifier.
        """
        import requests

        # Send the working directory to the server
        workdir = (
            os.getcwd()
            if not self.metadata.get("directory")
            else os.path.abspath(self.metadata["directory"])
        )

        # Prepare a request that has metadata about the job
        metadata = {
            "command": self.metadata.get("command"),
            "workdir": workdir,
        }

        response = requests.get(
            f"{self.address}/create_workflow",
            headers=self._headers,
            params=self.args,
            data=metadata,
        )

        # Extract the id from the response
        id = response.json()["id"]

        # Check the response, will exit on any error
        self.check_response(response, "/create_workflow")

        # Provide server parameters to the logger
        headers = (
            {"Content-Type": "application/json"}
            if self._headers is None
            else {**self._headers, **{"Content-Type": "application/json"}}
        )

        # Send the workflow name to the server
        response_change_workflow_name = requests.put(
            f"{self.address }/api/workflow/{id}",
            headers=headers,
            data=json.dumps(self.args),
        )
        # Check the response, will exit on any error
        self.check_response(response_change_workflow_name, f"/api/workflow/{id}")

        # Provide server parameters to the logger
        self.server = {"url": self.address, "id": id}

    def check_response(self, response, endpoint="wms monitor request"):
        """A helper function to take a response and check for an expected set of
        error codes, 404 (not found), 401 (requires authorization), 403 (permission
        denied), 500 (server error) and 200 (success).
        """
        status_code = response.status_code
        # Cut out early on success
        if status_code == 200:
            return

        if status_code == 404:
            sys.stderr.write(f"The wms {endpoint} endpoint was not found")
            sys.exit(-1)
        elif status_code == 401:
            sys.stderr.write(
                "Authorization is required with a WMS_MONITOR_TOKEN in the environment"
            )
            sys.exit(-1)
        elif status_code == 500:
            sys.stderr.write(
                f"There was a server error when trying to access {endpoint}"
            )
            sys.exit(-1)
        elif status_code == 403:
            sys.stderr.write("Permission is denied to %s." % endpoint)
            sys.exit(-1)

        # Any other response code is not acceptable
        sys.stderr.write(
            f"The {endpoint} response code {response.status_code} is not recognized."
        )

    @property
    def _headers(self):
        """return authenticated headers if the user has provided a token"""
        headers = None
        if self.token:
            headers = {"Authorization": "Bearer %s" % self.token}
        return headers

    def _parse_message(self, msg):
        """Given a message dictionary, we want to loop through the key, value
        pairs and convert some attributes to strings (e.g., jobs are fine to be
        represented as names) and return a dictionary.
        """
        result = {}
        for key, value in msg.items():
            # For a job, the name is sufficient
            if key == "job":
                result[key] = str(value)

            # For an exception, return the name and a message
            elif key == "exception":
                result[key] = "{}: {}".format(
                    msg["exception"].__class__.__name__,
                    msg["exception"] or "Exception",
                )

            # All other fields are json serializable
            else:
                result[key] = value

        # Return a json dumped string
        return json.dumps(result)

    def log_handler(self, msg):
        """Custom wms server log handler.

        Sends the log to the server.

        Args:
            msg (dict):    the log message dictionary
        """
        import requests

        url = self.server["url"] + "/update_workflow_status"
        server_info = {
            "msg": self._parse_message(msg),
            "timestamp": time.asctime(),
            "id": self.server["id"],
        }
        response = requests.post(url, data=server_info, headers=self._headers)
        self.check_response(response, "/update_workflow_status")


class DefaultFormatter(_logging.Formatter):

    def __init__(self, printreason=False, show_failed_logs=False, printshellcmds=False):
        self.printreason = printreason
        self.show_failed_logs = show_failed_logs
        self.printshellcmds = printshellcmds
        self.last_msg_was_job_info = False

    def format(self, record):
        """
        Override format method to handle Snakemake-specific log messages.
        """
        msg = record.msg
        level = msg.get("level", "INFO")  # Default to "INFO" if level not in message

        # Call specific handlers based on the log level
        if level == "info":
            return self.handle_info(msg)
        elif level == "host":
            return self.handle_host(msg)
        elif level == "job_info":
            return self.handle_job_info(msg)
        elif level == "group_info":
            return self.handle_group_info(msg)
        elif level == "job_error":
            return self.handle_job_error(msg)
        elif level == "group_error":
            return self.handle_group_error(msg)
        elif level == "progress":
            return self.handle_progress(msg)
        elif level == "job_finished":
            return self.handle_job_finished(msg)
        elif level == "shellcmd":
            return self.handle_shellcmd(msg)
        elif level == "rule_info":
            return self.handle_rule_info(msg)
        elif level == "d3dag":
            return self.handle_d3dag(msg)
        elif level == "dag_debug":
            return self.handle_dag_debug(msg)
        elif level == "run_info":
            return self.handle_run_info(msg)
        else:
            return msg["msg"]

    def handle_info(self, msg):
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

    def handle_run_info(self, msg):
        """Format the run_info log messages."""
        return msg["msg"]  # Log the message directly

    def handle_host(self, msg):
        """Format for host log."""
        return f"host: {platform.node()}"

    def handle_job_info(self, msg):
        """Format for job_info log."""
        output = []

        output.append(self.timestamp())
        if msg["msg"]:
            output.append(f"Job {msg['jobid']}: {msg['msg']}")
            if self.printreason:
                output.append(f"Reason: {msg['reason']}")
        else:
            output.append("\n".join(self.format_job_info(msg)))

        if msg["is_checkpoint"]:
            output.append("DAG of jobs will be updated after completion.")
        if msg["is_handover"]:
            output.append("Handing over execution to foreign system...")
        output.append("")
        return "\n".join(output)

    def handle_group_info(self, msg):
        """Format for group_info log."""
        msg = (
            f"{self.timestamp()} group job {msg['groupid']} (jobs in lexicogr. order):"
        )

        return msg

    def handle_job_error(self, msg):
        """Format for job_error log."""
        output = []
        output.append(self.timestamp())
        output.append("\n".join(self.format_job_error(msg)))
        return "\n".join(output)

    def handle_group_error(self, msg):
        """Format for group_error log."""
        output = []
        output.append(self.timestamp())
        output.append("\n".join(self.format_group_error(msg)))
        return "\n".join(output)

    def handle_progress(self, msg):
        """Format for progress log."""
        done = msg["done"]
        total = msg["total"]
        return f"{done} of {total} steps ({self.format_percentage(done, total)}) done"

    def handle_job_finished(self, msg):
        """Format for job_finished log."""
        return f"{self.timestamp()} Finished job {msg['jobid']}."

    def handle_shellcmd(self, msg):
        """Format for shellcmd log."""
        if self.printshellcmds:
            return msg["msg"]
        return ""

    def handle_rule_info(self, msg):
        """Format for rule_info log."""
        output = [msg["name"]]
        if msg["docstring"]:
            output.append(f"    {msg['docstring']}")
        return "\n".join(output)

    def handle_d3dag(self, msg):
        """Format for d3dag log."""
        return json.dumps({"nodes": msg["nodes"], "links": msg["edges"]})

    def handle_dag_debug(self, msg):
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

    def format_job_info(self, msg):
        """Helper method to format job info details."""

        def format_item(item, omit=None, valueformat=str):
            value = msg[item]
            if value != omit:
                return f"    {item}: {valueformat(value)}"

        output = [
            f"{'local' if msg['local'] else ''}{'checkpoint' if msg['is_checkpoint'] else 'rule'} {msg['name']}:"
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

    def format_job_error(self, msg):
        """Helper method to format job error details."""
        output = [f"Error in rule {msg['name']}:"]

        if msg["msg"]:
            output.append(f"    message: {msg['msg']}")
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

    def format_group_error(self, msg):
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

    def timestamp(self):
        """Helper method to format the timestamp."""
        return f"[{time.asctime()}]"

    def format_percentage(self, done, total):
        """Helper method to format percentage."""
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


class DefaultFilter:
    def __init__(
        self,
        quiet,
        debug_dag,
    ) -> None:
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

        msg = record.msg
        level = msg.get("level", "INFO").lower()

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


class ColorizingTextHandler(_logging.StreamHandler):
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
    ):
        super().__init__(stream=stream)
        self.last_msg_was_job_info = False
        self._output_lock = threading.Lock()
        self.nocolor = nocolor or not self.can_color_tty(mode)
        self.mode = mode

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

                level = record.msg.get("level", "INFO").lower()

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


class Logger:

    def __init__(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        self.logger = _logging.getLogger(__name__)
        self.stream_handler = None
        self.printshellcmds = False
        self.printreason = False
        self.debug_dag = False
        self.quiet = set()
        self.logfile = None
        self.last_msg_was_job_info = False
        self.mode = ExecMode.DEFAULT
        self.show_failed_logs = False
        self.logfile_handler = None
        self.dryrun = False
        self.level = 0
        self.default_formatter = None
        self.default_filter = None

    def setup_logfile(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT and not self.dryrun:
            os.makedirs(os.path.join(".snakemake", "log"), exist_ok=True)
            self.logfile = os.path.abspath(
                os.path.join(
                    ".snakemake",
                    "log",
                    datetime.datetime.now().isoformat().replace(":", "")
                    + ".snakemake.log",
                )
            )

            self.logfile_handler = _logging.FileHandler(self.logfile)
            self.logfile_handler.setFormatter(self.default_formatter)
            self.logfile_handler.addFilter(self.default_filter)
            self.logfile_handler.setLevel(self.level)
            self.logger.addHandler(self.logfile_handler)

    def cleanup(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT and self.logfile_handler is not None:
            self.logger.removeHandler(self.logfile_handler)
            self.logfile_handler.close()

    def set_stream_handler(self, stream_handler: _logging.Handler):
        """Set the stream handler, replacing any existing one."""
        if self.stream_handler is not None:
            self.logger.removeHandler(self.stream_handler)
        self.stream_handler = stream_handler
        stream_handler.setFormatter(self.default_formatter)
        stream_handler.addFilter(self.default_filter)
        stream_handler.setLevel(self.level)
        self.logger.addHandler(stream_handler)

    def set_level(self, level):
        """Set the logging level."""
        self.level = level
        self.logger.setLevel(level)

    def handler(self, msg):
        """Send the message to the logger."""
        # Get the custom log level from the message
        custom_level = msg.get("level", "info").lower()

        # Map custom levels to standard Python logging levels
        log_level = {
            "debug": _logging.DEBUG,
            "info": _logging.WARNING,
            "warning": _logging.WARNING,
            "error": _logging.ERROR,
            "critical": _logging.CRITICAL,
            "run_info": _logging.WARNING,
            "job_info": _logging.INFO,
            "group_info": _logging.INFO,
            "job_error": _logging.ERROR,
            "group_error": _logging.ERROR,
            "progress": _logging.INFO,
            "job_finished": _logging.INFO,
            "shellcmd": _logging.INFO,
            "rule_info": _logging.INFO,
            "dag_debug": (
                _logging.WARNING if self.debug_dag else _logging.DEBUG
            ),  # bump debug_dag to warning if we want it
            "resources_info": _logging.WARNING,
            "host": _logging.INFO,
            "job_stats": _logging.WARNING,
        }.get(
            custom_level, _logging.INFO
        )  # Default to INFO if not recognized

        record = self.logger.makeRecord(
            name=self.logger.name,
            level=log_level,
            fn="",
            lno=0,
            msg=msg,
            args=(),
            exc_info=None,
        )

        self.logger.handle(record)

    # Logging methods for various log levels
    def info(self, msg, indent=False):
        self.handler(dict(level="info", msg=msg, indent=indent))

    def warning(self, msg, *fmt_items):
        if fmt_items:
            msg = msg % fmt_items
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

    def group_info(self, **msg):
        msg["level"] = "group_info"
        self.handler(msg)

    def job_info(self, **msg):
        msg["level"] = "job_info"
        self.handler(msg)

    def job_error(self, **msg):
        msg["level"] = "job_error"
        self.handler(msg)

    def group_error(self, **msg):
        msg["level"] = "group_error"
        self.handler(msg)

    def shellcmd(self, msg, indent=False):
        if msg is not None:
            msg = dict(level="shellcmd", msg=msg)
            msg["indent"] = indent
            self.handler(msg)

    def job_finished(self, **msg):
        msg["level"] = "job_finished"
        self.handler(msg)

    def rule_info(self, **msg):
        msg["level"] = "rule_info"
        self.handler(msg)

    def d3dag(self, **msg):
        msg["level"] = "d3dag"
        self.handler(msg)

    def dag_debug(self, msg):
        self.handler(dict(level="dag_debug", **msg))

    def host_info(self):
        self.handler(dict(level="host"))

    def logfile_hint(self):
        """Log the logfile location if applicable."""
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT and not self.dryrun:
            logfile = self.get_logfile()
            self.info(f"Complete log: {os.path.relpath(logfile)}")

    def get_logfile(self):
        """Get the path to the logfile."""
        if self.logfile is not None:
            self.logfile_handler.flush()
        return self.logfile


logger = Logger()

def setup_logger(
    handler=[],
    quiet=False,
    printshellcmds=False,
    printreason=True,
    debug_dag=False,
    nocolor=False,
    stdout=False,
    debug=False,
    mode=None,
    show_failed_logs=False,
    dryrun=False,
):
    from snakemake.settings.types import Quietness

    if mode is None:
        mode = get_default_exec_mode()

    if quiet is None:
        # not quiet at all
        quiet = set()
    elif isinstance(quiet, bool):
        if quiet:
            quiet = {Quietness.PROGRESS, Quietness.RULES}
        else:
            quiet = set()
    elif not isinstance(quiet, set):
        raise ValueError(
            "Unsupported value provided for quiet mode (either bool, None or set allowed)."
        )

    stream_handler = ColorizingTextHandler(
        nocolor=nocolor,
        stream=sys.stdout if stdout else sys.stderr,
        mode=mode,
    )
    formatter = DefaultFormatter(
        printreason=printreason,
        show_failed_logs=show_failed_logs,
        printshellmds=printshellcmds,
    )
    filter = DefaultFilter(quiet=quiet, debug_dag=debug_dag)
    logger.default_formatter = formatter
    logger.default_filter = filter
    logger.set_stream_handler(stream_handler)
    logger.set_level(_logging.DEBUG if debug else _logging.INFO)
    logger.quiet = quiet
    logger.printshellcmds = printshellcmds
    logger.printreason = printreason
    logger.debug_dag = debug_dag
    logger.mode = mode
    logger.dryrun = dryrun
    logger.show_failed_logs = show_failed_logs
