# -*- coding: utf-8 -*-

import logging
import platform
import time
import sys
from multiprocessing import Lock

__author__ = "Johannes Köster"


class ColorizingStreamHandler(logging.StreamHandler):
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
        self.nocolor = nocolor or not self.is_tty or platform.system() == 'Windows'
        self.timestamp = timestamp

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
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                self.handleError(record)

    def decorate(self, record):
        message = [record.message]
        if self.timestamp:
            message.insert(0, "[{}] ".format(time.asctime()))
        if not self.nocolor and record.levelname in self.colors:
            message.insert(0, self.COLOR_SEQ % (30 + self.colors[record.levelname]))
            message.append(self.RESET_SEQ)
        return "".join(message)



logger = Logger()
stream_handler = None
handler = None


class Logger:
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def info(self, msg):
        handler(dict(level="info", msg=msg))

    def debug(self, msg):
        handler(dict(level="debug", msg=msg))

    def critical(self, msg):
        handler(dict(level="critical", msg=msg))

    def progress(self, count, total):
        handler(dict(level="progress", count=count, total=total))

    def job_info(self, **msg):
        msg["level"] = "job_info"
        handler(msg)
    

def quiet_console_handler(msg):
    level = msg["level"]
    if level == "info":
        logger.warning(msg["msg"])
    elif level == "critical":
        logger.critical(msg["msg"])
    elif level == "debug":
        logger.debug(msg["msg"])


def console_handler(msg):
    def job_info(msg):
        def format_item(item, omit=None, valueformat=str):
            value = msg[item]
            if value != omit:
                return "\t{}: {}".format(item, valueformat(value))
                
        yield "{}rule {}:".format("local" if msg["local"] else "", msg["name"])
        for item in "input output".split():
            yield format_item(item, omit=[], ", ".join)
        for item in "log reason".split():
            yield format_item(item, omit="")
        for item in "priority threads".split():
            yield format_item(item, omit=1)
        
    quiet_console_handler(msg)
    if level == "progress":
        done = msg["done"]
        total = msg["total"]
        logger.info("{} of {} steps ({:.0%}) done".format(done, total, done / total))
    elif level == "job_info":
        if "msg" in msg:
            logger.info(msg["msg"])
        else:
            logger.info("\n".join(job_info(msg)))


def init_logger(log_handler=console_handler, nocolor=False, stdout=False, debug=False, timestamp=False):
    global logger
    global stream_handler
    global handler
    handler = log_handler
    if stream_handler:
        logger.removeHandler(stream_handler)
    stream_handler = ColorizingStreamHandler(
        nocolor=nocolor, stream=sys.stdout if stdout else sys.stderr,
        timestamp=timestamp
    )
    logger.logger.addHandler(stream_handler)
    logger.logger.setLevel(logging.DEBUG if debug else logging.INFO)


init_logger()
