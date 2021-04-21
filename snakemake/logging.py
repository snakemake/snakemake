__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
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
import tempfile
from functools import partial
import inspect
import traceback
import textwrap

from snakemake.common import DYNAMIC_FILL
from snakemake.common import Mode


class ColorizingStreamHandler(_logging.StreamHandler):

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[%dm"
    BOLD_SEQ = "\033[1m"

    colors = {
        "WARNING": YELLOW,
        "INFO": GREEN,
        "DEBUG": BLUE,
        "CRITICAL": RED,
        "ERROR": RED,
    }

    def __init__(
        self, nocolor=False, stream=sys.stderr, use_threads=False, mode=Mode.default
    ):
        super().__init__(stream=stream)

        self._output_lock = threading.Lock()

        self.nocolor = nocolor or not self.can_color_tty(mode)

    def can_color_tty(self, mode):
        if "TERM" in os.environ and os.environ["TERM"] == "dumb":
            return False
        if mode == Mode.subprocess:
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
        from slacker import Slacker

        self.token = os.getenv("SLACK_TOKEN")
        if not self.token:
            print(
                "The use of slack logging requires the user to set a user specific slack legacy token to the SLACK_TOKEN environment variable. Set this variable by 'export SLACK_TOKEN=your_token'. To generate your token please visit https://api.slack.com/custom-integrations/legacy-tokens."
            )
            exit(-1)
        self.slack = Slacker(self.token)
        # Check for success
        try:
            auth = self.slack.auth.test().body
        except Exception:
            print(
                "Slack connection failed. Please compare your provided slack token exported in the SLACK_TOKEN environment variable with your online token at https://api.slack.com/custom-integrations/legacy-tokens. A different token can be set up by 'export SLACK_TOKEN=your_token'."
            )
            exit(-1)
        self.own_id = auth["user_id"]
        self.error_occured = False

    def log_handler(self, msg):
        if msg["level"] == "error" and not self.error_occured:
            self.slack.chat.post_message(
                self.own_id, text="At least one error occured.", username="snakemake"
            )
            self.error_occured = True

        if msg["level"] == "progress" and msg["done"] == msg["total"]:
            # workflow finished
            self.slack.chat.post_message(
                self.own_id, text="Workflow complete.", username="snakemake"
            )


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

        from snakemake.resources import parse_resources

        self.address = address or "http:127.0.0.1:5000"
        self.args = parse_resources(args) or []
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
            self.address + "/api/service-info", headers=self._headers
        )
        if response.status_code != 200:
            sys.stderr.write(
                "Problem with server: {} {}".format(self.address, os.linesep)
            )
            sys.exit(-1)

        # And then that it's ready to be interacted with
        if response.json().get("status") != "running":
            sys.stderr.write(
                "The status of the server {} is not in 'running' mode {}".format(
                    self.address, os.linesep
                )
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
            "snakefile": os.path.join(workdir, self.metadata.get("snakefile")),
            "command": self.metadata.get("command"),
            "workdir": workdir,
        }

        response = requests.get(
            self.address + "/create_workflow",
            headers=self._headers,
            params=self.args,
            data=metadata,
        )

        # Check the response, will exit on any error
        self.check_response(response, "/create_workflow")

        # Provide server parameters to the logger
        self.server = {"url": self.address, "id": response.json()["id"]}

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
            sys.stderr.write("The wms %s endpoint was not found" % endpoint)
            sys.exit(-1)
        elif status_code == 401:
            sys.stderr.write(
                "Authorization is required with a WMS_MONITOR_TOKEN in the environment"
            )
            sys.exit(-1)
        elif status_code == 500:
            sys.stderr.write(
                "There was a server error when trying to access %s" % endpoint
            )
            sys.exit(-1)
        elif status_code == 403:
            sys.stderr.write("Permission is denied to %s." % endpoint)
            sys.exit(-1)

        # Any other response code is not acceptable
        sys.stderr.write(
            "The %s response code %s is not recognized."
            % (endpoint, response.status_code)
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
        represnted as names) and return a dictionary.
        """
        result = {}
        for key, value in msg.items():

            # For a job, the name is sufficient
            if key == "job":
                result[key] = str(value)

            # For an exception, return the name and a message
            elif key == "exception":
                result[key] = "%s: %s" % (
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
            msg (dict):     the log message dictionary
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
        self.logger = _logging.getLogger(__name__)
        self.log_handler = [self.text_handler]
        self.stream_handler = None
        self.printshellcmds = False
        self.printreason = False
        self.debug_dag = False
        self.quiet = False
        self.logfile = None
        self.last_msg_was_job_info = False
        self.mode = Mode.default
        self.show_failed_logs = False
        self.logfile_handler = None

    def setup_logfile(self):
        if self.mode == Mode.default:
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
        if self.mode == Mode.default and self.logfile_handler is not None:
            self.logger.removeHandler(self.logfile_handler)
            self.logfile_handler.close()
        self.log_handler = [self.text_handler]

    def get_logfile(self):
        if self.logfile is not None:
            self.logfile_handler.flush()
        return self.logfile

    def remove_logfile(self):
        if self.mode == Mode.default:
            self.logfile_handler.close()
            os.remove(self.logfile)

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

    def logfile_hint(self):
        if self.mode == Mode.default:
            logfile = self.get_logfile()
            self.info("Complete log: {}".format(logfile))

    def location(self, msg):
        callerframerecord = inspect.stack()[1]
        frame = callerframerecord[0]
        info = inspect.getframeinfo(frame)
        self.debug(
            "{}: {info.filename}, {info.function}, {info.lineno}".format(msg, info=info)
        )

    def info(self, msg, indent=False):
        self.handler(dict(level="info", msg=msg, indent=indent))

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
                    return "    {}: {}".format(item, valueformat(value))

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
                yield "    wildcards: " + wildcards

            for item, omit in zip("priority threads".split(), [0, 1]):
                fmt = format_item(item, omit=omit)
                if fmt != None:
                    yield fmt

            resources = format_resources(msg["resources"])
            if resources:
                yield "    resources: " + resources

        def indent(item):
            if msg.get("indent"):
                return "    " + item
            else:
                return item

        def timestamp():
            self.logger.info(indent("[{}]".format(time.asctime())))

        level = msg["level"]

        if level == "job_info" and not self.quiet:
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
                    indent("Downstream jobs will be updated " "after completion.")
                )
            if msg["is_handover"]:
                self.logger.warning("Handing over execution to foreign system...")
            self.logger.info("")

            self.last_msg_was_job_info = True
        elif level == "group_info" and not self.quiet:
            timestamp()
            if not self.last_msg_was_job_info:
                self.logger.info("")
            self.logger.info(
                "group job {} (jobs in lexicogr. order):".format(msg["groupid"])
            )
        elif level == "job_error":
            timestamp()
            self.logger.error(indent("Error in rule {}:".format(msg["name"])))
            self.logger.error(indent("    jobid: {}".format(msg["jobid"])))
            if msg["output"]:
                self.logger.error(
                    indent("    output: {}".format(", ".join(msg["output"])))
                )
            if msg["log"]:
                self.logger.error(
                    indent(
                        "    log: {} (check log file(s) for error message)".format(
                            ", ".join(msg["log"])
                        )
                    )
                )
            if msg["conda_env"]:
                self.logger.error(indent("    conda-env: {}".format(msg["conda_env"])))
            if msg["shellcmd"]:
                self.logger.error(
                    indent(
                        "    shell:\n        {}\n        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)".format(
                            msg["shellcmd"]
                        )
                    )
                )

            for item in msg["aux"].items():
                self.logger.error(indent("    {}: {}".format(*item)))

            if self.show_failed_logs and msg["log"]:
                for f in msg["log"]:
                    try:
                        self.logger.error("Logfile {}:\n{}".format(f, open(f).read()))
                    except FileNotFoundError:
                        self.logger.error("Logfile {} not found.".format(f))

            self.logger.error("")
        elif level == "group_error":
            timestamp()
            self.logger.error("Error in group job {}:".format(msg["groupid"]))
        else:
            if level == "info" and not self.quiet:
                self.logger.warning(msg["msg"])
            if level == "warning":
                self.logger.warning(msg["msg"])
            elif level == "error":
                self.logger.error(msg["msg"])
            elif level == "debug":
                self.logger.debug(msg["msg"])
            elif level == "resources_info" and not self.quiet:
                self.logger.warning(msg["msg"])
            elif level == "run_info":
                self.logger.warning(msg["msg"])
            elif level == "progress" and not self.quiet:
                done = msg["done"]
                total = msg["total"]
                p = done / total
                percent_fmt = ("{:.2%}" if p < 0.01 else "{:.0%}").format(p)
                self.logger.info(
                    "{} of {} steps ({}) done".format(done, total, percent_fmt)
                )
            elif level == "shellcmd":
                if self.printshellcmds:
                    self.logger.warning(indent(msg["msg"]))
            elif level == "job_finished" and not self.quiet:
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


def format_dict(dict_like, omit_keys=[], omit_values=[]):
    from snakemake.io import Namedlist

    if isinstance(dict_like, Namedlist):
        items = dict_like.items()
    elif isinstance(dict_like, dict):
        items = dict_like.items()
    else:
        raise ValueError(
            "bug: format_dict applied to something neither a dict nor a Namedlist"
        )
    return ", ".join(
        "{}={}".format(name, str(value))
        for name, value in items
        if name not in omit_keys and value not in omit_values
    )


format_resources = partial(format_dict, omit_keys={"_cores", "_nodes"})
format_wildcards = partial(format_dict, omit_values={DYNAMIC_FILL})


def format_resource_names(resources, omit_resources="_cores _nodes".split()):
    return ", ".join(name for name in resources if name not in omit_resources)


logger = Logger()


def setup_logger(
    handler=[],
    quiet=False,
    printshellcmds=False,
    printreason=False,
    debug_dag=False,
    nocolor=False,
    stdout=False,
    debug=False,
    use_threads=False,
    mode=Mode.default,
    show_failed_logs=False,
):
    logger.log_handler.extend(handler)

    # console output only if no custom logger was specified
    stream_handler = ColorizingStreamHandler(
        nocolor=nocolor,
        stream=sys.stdout if stdout else sys.stderr,
        use_threads=use_threads,
        mode=mode,
    )
    logger.set_stream_handler(stream_handler)
    logger.set_level(_logging.DEBUG if debug else _logging.INFO)
    logger.quiet = quiet
    logger.printshellcmds = printshellcmds
    logger.printreason = printreason
    logger.debug_dag = debug_dag
    logger.mode = mode
    logger.show_failed_logs = show_failed_logs
