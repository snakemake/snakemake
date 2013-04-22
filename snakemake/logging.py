# -*- coding: utf-8 -*-

import logging
import platform
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

    def __init__(self, nocolor=False, stream=sys.stderr):
        super().__init__(stream=stream)
        self.nocolor = nocolor

    @property
    def is_tty(self):
        isatty = getattr(self.stream, 'isatty', None)
        return isatty and isatty()

    def emit(self, record):
        try:
            self.format(record)  # add the message to the record
            self._output_lock.acquire()
            if self.is_tty:
                self.stream.write(self.colorize(record))
            else:
                self.stream.write(record.message)
            self.stream.write(getattr(self, 'terminator', '\n'))
            self.flush()
            self._output_lock.release()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)
            self._output_lock.release()

    def colorize(self, record):
        if (not self.nocolor
            and record.levelname in self.colors
            and platform.system() != 'Windows'):
            return "{color}{message}{reset}".format(
                color=self.COLOR_SEQ % (30 + self.colors[record.levelname]),
                message=record.message,
                reset=self.RESET_SEQ
            )
        return record.message

logger = logging.getLogger(__name__)
handler = None


def init_logger(nocolor=False, stdout=False, debug=False):
    global logger
    global handler
    if handler:
        logger.removeHandler(handler)
    handler = ColorizingStreamHandler(
        nocolor=nocolor, stream=sys.stdout if stdout else sys.stderr)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
