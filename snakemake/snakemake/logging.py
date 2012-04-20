# -*- coding: utf-8 -*-

import logging, string, platform
from multiprocessing import Lock

__author__ = "Johannes Köster"

class ColorizingStreamHandler(logging.StreamHandler):
	nocolor = False
	quiet = False
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

	@property
	def is_tty(self):
		isatty = getattr(self.stream, 'isatty', None)
		return isatty and isatty()

	def emit(self, record):
		try:
			if self.quiet and record.levelname in ('INFO', 'DEBUG'): return
			message = self.format(record)
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
	
	@classmethod
	def colorize(cls, record):
		if not cls.nocolor and record.levelname in cls.colors and platform.system() != 'Windows':
			return "{color}{message}{reset}".format(
				color = cls.COLOR_SEQ % (30 + cls.colors[record.levelname]),
				message = record.message,
				reset = cls.RESET_SEQ
			)
		return record.message

logger = logging.getLogger(__name__)
logger.addHandler(ColorizingStreamHandler())
logger.setLevel(logging.INFO)
	
