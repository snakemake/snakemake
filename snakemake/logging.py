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


class ColorizingStreamHandler(_logging.StreamHandler):
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

    def __init__(self, nocolor=False, stream=sys.stderr, mode=None):
        super().__init__(stream=stream)

        if mode is None:
            mode = get_default_exec_mode()

        self._output_lock = threading.Lock()

        self.nocolor = nocolor or not self.can_color_tty(mode)

    def can_color_tty(self, mode):
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
        with self._output_lock:
            try:
                self.format(record)  # add the message to the record
                self.stream.write(self.decorate(record))
                self.stream.write(getattr(self, "terminator", "\n"))
                self.flush()
            except BrokenPipeError as e:
                raise e
            except (KeyboardInterrupt, SystemExit):
                # ignore any exceptions in these cases as any relevant messages have been printed before
                pass
            except Exception as e:
                self.handleError(record)

    def decorate(self, record):
        message = record.message
        message = [message]
        if not self.nocolor and record.levelname in self.colors:
            message.insert(0, self.COLOR_SEQ % (30 + self.colors[record.levelname]))
            message.append(self.RESET_SEQ)
        return "".join(message)


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


class Logger:
    def __init__(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        self.logger = _logging.getLogger(__name__)
        self.log_handler = [self.text_handler]
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
            self.logger.addHandler(self.logfile_handler)

    def cleanup(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT and self.logfile_handler is not None:
            self.logger.removeHandler(self.logfile_handler)
            self.logfile_handler.close()
        self.log_handler = [self.text_handler]

    def get_logfile(self):
        if self.logfile is not None:
            self.logfile_handler.flush()
        return self.logfile

    def handler(self, msg):
        msg["timestamp"] = time.time()
        for handler in self.log_handler:
            handler(msg)

    def set_stream_handler(self, stream_handler):
        if self.stream_handler is not None:
            self.logger.removeHandler(self.stream_handler)
        self.stream_handler = stream_handler
        self.logger.addHandler(stream_handler)

    def set_level(self, level):
        self.logger.setLevel(level)

    def logfile_hint(self):
        from snakemake_interface_executor_plugins.settings import ExecMode

        if self.mode == ExecMode.DEFAULT and not self.dryrun:
            logfile = self.get_logfile()
            self.info(f"Complete log: {os.path.relpath(logfile)}")

    def location(self, msg):
        callerframerecord = inspect.stack()[1]
        frame = callerframerecord[0]
        info = inspect.getframeinfo(frame)
        self.debug(
            "{}: {info.filename}, {info.function}, {info.lineno}".format(msg, info=info)
        )

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

    def dag_debug(self, msg):
        self.handler(dict(level="dag_debug", **msg))

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

    def host_info(self):
        self.handler(dict(level="host"))

    def is_quiet_about(self, msg_type: str):
        from snakemake.settings.enums import Quietness

        return (
            Quietness.ALL in self.quiet
            or Quietness.parse_choice(msg_type) in self.quiet
        )

    def text_handler(self, msg):
        """The default snakemake log handler.

        Prints the output to the console.

        Args:
            msg (dict):     the log message dictionary
        """
        if self.is_quiet_about("all"):
            # do not log anything
            return

        def job_info(msg):
            def format_item(item, omit=None, valueformat=str):
                value = msg[item]
                if value != omit:
                    return f"    {item}: {valueformat(value)}"

            yield "{}{} {}:".format(
                "local" if msg["local"] else "",
                "checkpoint" if msg["is_checkpoint"] else "rule",
                msg["name"],
            )
            for item in ["input", "output", "log"]:
                fmt = format_item(item, omit=[], valueformat=", ".join)
                if fmt != None:
                    yield fmt

            singleitems = ["jobid", "benchmark"]
            if self.printreason:
                singleitems.append("reason")
            for item in singleitems:
                fmt = format_item(item, omit=None)
                if fmt != None:
                    yield fmt

            wildcards = format_wildcards(msg["wildcards"])
            if wildcards:
                yield f"    wildcards: {wildcards}"

            for item, omit in zip("priority threads".split(), [0, 1]):
                fmt = format_item(item, omit=omit)
                if fmt != None:
                    yield fmt

            resources = format_resources(msg["resources"])
            if resources:
                yield f"    resources: {resources}"

        def show_logs(logs):
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
                # take the length of the longest line, but limit to max 80
                max_len = min(max(max(len(l) for l in lines), len(logfile_header)), 80)
                yield "=" * max_len
                yield from lines
                yield "=" * max_len

        def indent(item):
            if msg.get("indent"):
                return f"    {item}"
            else:
                return item

        def timestamp():
            self.logger.info(indent(f"[{time.asctime()}]"))

        level = msg["level"]

        if level == "host" and not self.is_quiet_about("host"):
            self.logger.info(f"host: {platform.node()}")
        elif level == "job_info" and not self.is_quiet_about("rules"):
            if not self.last_msg_was_job_info:
                self.logger.info("")
            timestamp()
            if msg["msg"] is not None:
                self.logger.info(indent("Job {}: {}".format(msg["jobid"], msg["msg"])))
                if self.printreason:
                    self.logger.info(indent("Reason: {}".format(msg["reason"])))
            else:
                self.logger.info("\n".join(map(indent, job_info(msg))))
            if msg["is_checkpoint"]:
                self.logger.warning(
                    indent("DAG of jobs will be updated after completion.")
                )
            if msg["is_handover"]:
                self.logger.warning("Handing over execution to foreign system...")
            self.logger.info("")

            self.last_msg_was_job_info = True
        elif level == "group_info" and not self.is_quiet_about("rules"):
            timestamp()
            msg = "group job {} (jobs in lexicogr. order):".format(msg["groupid"])
            if not self.last_msg_was_job_info:
                msg = "\n" + msg
            self.logger.info(msg)
        elif level == "job_error":

            def job_error():
                yield "Error in rule {}:".format(msg["name"])
                if msg["msg"]:
                    yield "    message: {}".format(msg["msg"])
                yield "    jobid: {}".format(msg["jobid"])
                if msg["input"]:
                    yield "    input: {}".format(", ".join(msg["input"]))
                if msg["output"]:
                    yield "    output: {}".format(", ".join(msg["output"]))
                if msg["log"]:
                    yield "    log: {} (check log file(s) for error details)".format(
                        ", ".join(msg["log"])
                    )
                if msg["conda_env"]:
                    yield "    conda-env: {}".format(msg["conda_env"])
                if msg["shellcmd"]:
                    yield "    shell:\n        {}\n        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)".format(
                        msg["shellcmd"]
                    )

                for item in msg["aux"].items():
                    yield "    {}: {}".format(*item)

                if self.show_failed_logs and msg["log"]:
                    yield from show_logs(msg["log"])

                yield ""

            timestamp()
            self.logger.error("\n".join(map(indent, job_error())))
        elif level == "group_error":

            def group_error():
                yield f"Error in group {msg['groupid']}:"
                if msg["msg"]:
                    yield f"    message: {msg['msg']}"
                if msg["aux_logs"]:
                    yield f"    log: {', '.join(msg['aux_logs'])} (check log file(s) for error details)"
                yield "    jobs:"
                for info in msg["job_error_info"]:
                    yield f"        rule {info['name']}:"
                    yield f"            jobid: {info['jobid']}"
                    if info["output"]:
                        yield f"            output: {', '.join(info['output'])}"
                    if info["log"]:
                        yield f"            log: {', '.join(info['log'])} (check log file(s) for error details)"
                logs = msg["aux_logs"] + [
                    f for info in msg["job_error_info"] for f in info["log"]
                ]
                if self.show_failed_logs and logs:
                    yield from show_logs(logs)
                yield ""

            timestamp()
            self.logger.error("\n".join(group_error()))
        else:
            if level == "info" and not self.is_quiet_about("progress"):
                self.logger.warning(msg["msg"])
            if level == "warning":
                self.logger.critical(msg["msg"])
            elif level == "error":
                self.logger.error(msg["msg"])
            elif level == "debug":
                self.logger.debug(msg["msg"])
            elif level == "resources_info" and not self.is_quiet_about("progress"):
                self.logger.warning(msg["msg"])
            elif level == "run_info" and not self.is_quiet_about("progress"):
                self.logger.warning(msg["msg"])
            elif level == "progress" and not self.is_quiet_about("progress"):
                done = msg["done"]
                total = msg["total"]
                self.logger.info(
                    "{} of {} steps ({}) done".format(
                        done, total, format_percentage(done, total)
                    )
                )
            elif level == "shellcmd":
                if self.printshellcmds:
                    self.logger.warning(indent(msg["msg"]))
            elif level == "job_finished" and not self.is_quiet_about("progress"):
                timestamp()
                self.logger.info("Finished job {}.".format(msg["jobid"]))
                pass
            elif level == "rule_info":
                self.logger.info(msg["name"])
                if msg["docstring"]:
                    self.logger.info("    " + msg["docstring"])
            elif level == "d3dag":
                print(json.dumps({"nodes": msg["nodes"], "links": msg["edges"]}))
            elif level == "dag_debug":
                if self.debug_dag:
                    if "file" in msg:
                        self.logger.warning(
                            "file {file}:\n    {msg}\n{exception}".format(
                                file=msg["file"],
                                msg=msg["msg"],
                                exception=textwrap.indent(
                                    str(msg["exception"]), "    "
                                ),
                            )
                        )
                    else:
                        job = msg["job"]
                        self.logger.warning(
                            "{status} job {name}\n    wildcards: {wc}".format(
                                status=msg["status"],
                                name=job.rule.name,
                                wc=format_wildcards(job.wildcards),
                            )
                        )

            self.last_msg_was_job_info = False


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

    logger.log_handler.extend(handler)

    # console output only if no custom logger was specified
    stream_handler = ColorizingStreamHandler(
        nocolor=nocolor,
        stream=sys.stdout if stdout else sys.stderr,
        mode=mode,
    )
    logger.set_stream_handler(stream_handler)
    logger.set_level(_logging.DEBUG if debug else _logging.INFO)
    logger.quiet = quiet
    logger.printshellcmds = printshellcmds
    logger.printreason = printreason
    logger.debug_dag = debug_dag
    logger.mode = mode
    logger.dryrun = dryrun
    logger.show_failed_logs = show_failed_logs
