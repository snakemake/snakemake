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

logger = logging.getLogger(__name__)
handler = None


def init_logger(nocolor=False, stdout=False, debug=False, timestamp=False):
    global logger
    global handler
    if handler:
        logger.removeHandler(handler)
    handler = ColorizingStreamHandler(
        nocolor=nocolor, stream=sys.stdout if stdout else sys.stderr,
        timestamp=timestamp
    )
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

init_logger()
