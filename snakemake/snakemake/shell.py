# -*- coding: utf-8 -*-

import _io
import signal
import sys, os, inspect, time, shlex
import subprocess as sp
from threading import Thread
from snakemake.exceptions import TerminatedException
import atexit

__author__ = "Johannes Köster"

class SequenceFormatter:
	def __init__(self, sequence):
		self._sequence = sequence

	def __getitem__(self, i):
		return self._sequence[i]

	def __str__(self):
		return " ".join(self._sequence)



def format(string, *args, stepout = 1, **kwargs):
	frame = inspect.currentframe().f_back
	while stepout > 1:
		if not frame.f_back:
			break
		frame = frame.f_back
		stepout -= 1
	
	variables = dict(frame.f_globals)
	# add local variables from calling rule/function
	variables.update(frame.f_locals)
	variables.update(kwargs)
	strmethods = list()
	for key, value in list(variables.items()):
		if type(value) in (list, tuple, set, frozenset):
			variables[key] = SequenceFormatter(value)
	try:
		return string.format(*args, **variables)
	except KeyError as ex:
		raise NameError("The name {} is unknown in this context.".format(str(ex)))

class PipeWriter:
	def __init__(self, towrite):
		self._towrite = towrite
		
	def write(self, line):
		self._towrite.write(line)

class ListWriter(PipeWriter):
	def write(self, line):
		self._towrite.append(line[:-1].decode("utf-8"))

class TextIOWriter(PipeWriter):
	def write(self, line):
		self._towrite.write(line.decode("utf-8"))
		
class shell(sp.Popen):
	_process_args = {}
	_process_prefix = ""
	_processes = []

	@classmethod
	def executable(cls, cmd):
		cls._process_args["executable"] = cmd

	@classmethod
	def prefix(cls, prefix):
		cls._process_prefix = format(prefix, stepout = 2)
	
	def __init__(self, cmd, *args, async = False, iterable = False, **kwargs):
		if async or iterable:
			stdout, stdin = (sp.PIPE, sp.PIPE)
		else:
			stdout, stdin = (sys.stdout, None)
			if not isinstance(sys.stdout, _io.TextIOWrapper):
				# workaround for nosetest since it overwrites sys.stdout in a strange way that does not work with Popen
				stdout = None

		self.cmd = format(cmd, *args, stepout = 2, **kwargs)
		super().__init__(self._process_prefix + self.cmd, shell=True, stdin = stdin, stdout=stdout, close_fds=True, **shell._process_args)
		
		self.async = async or iterable
		self._stdout_free = True
		self._pipethread = None
		self._stdin = self.stdin

		self._prog_term = False
		
		shell._processes.append(self)
		
		try:
			signal.signal(signal.SIGTERM, self.terminate_all)
		except ValueError:
			# if signal handling cannot be set, ignore it. Snakemake will terminate for processes that are not ill-behaving anyway.
			pass
		if not self.async:
			self.wait()

	def wait(self):
		#ret = super().wait()
		ret = self.poll()
		while ret == None:
			if self._prog_term:
				try:
					self.kill()
				except OSError:
					pass
				raise sp.CalledProcessError(-2, self.cmd)
			time.sleep(0.1)
			ret = self.poll()
		if self._pipethread:
			self._pipethread.join()
		if self._prog_term:
			raise TerminatedException()
		if ret != 0:
			raise sp.CalledProcessError(ret, self.cmd)
		
	@staticmethod
	def join_all():
		for p in shell._processes:
			if p.async:
				if p._stdout_free:
					for l in p._stdoutlines(): pass
				p.wait()
		shell._processes = []

	@staticmethod
	def terminate_all(*args):
		for p in shell._processes:
			p._prog_term = True
		
	def _stdoutlines(self):
		while True:
			o = self.stdout.readline()
			if not o:
				return
			yield o

	def _write_pipes(self, pipes, toclose):
		for l in self._stdoutlines():
			for pipe in pipes:
				pipe.write(l)
		for pipe in toclose:
			pipe.close()
				
	def __or__(self, other):
		if not self.async:
			raise SyntaxError("Pipe operator \"|\" not allowed for synchronous shell() calls (i.e. without async=True).")
		if not isinstance(other, tuple):
			other = tuple([other])

		pipes = []
		toclose = []
		for o in other:
			writer = None
			if isinstance(o, shell):
				writer = PipeWriter(o._stdin)
				toclose.append(o._stdin)
				# move all stdin to the beginning of the pipe
				o._stdin = self._stdin
			elif isinstance(o, _io._TextIOBase):
				writer = TextIOWriter(o)
			elif isinstance(o, list):
				writer = ListWriter(o)
			else:
				raise ValueError("Only shell, files, stdout or lists allowed right to a shell pipe.")
			pipes.append(writer)
		
		self._stdout_free = False
		self._pipethread = Thread(target = self._write_pipes, args = (pipes, toclose))
		self._pipethread.daemon = True
		self._pipethread.start()
			
		if len(other) == 1:
			return other[0]

	def __iter__(self):
		buf = list()
		self | buf
		shell.join_all()
		return buf.__iter__()

atexit.register(shell.terminate_all)

if "SHELL" in os.environ:
	#shell._process_args["executable"] = os.environ["SHELL"]
	shell.executable(os.environ["SHELL"])

if __name__ == "__main__":
	shell("echo b; echo a; echo c > foo", async=True)
	shell("echo b; echo a; echo c", async=True) | shell("sort", async=True) | sys.stdout
	
	shell("echo a; echo b", async=True) | (shell("sort", async=True), shell("cut -f1", async=True))
	
	x = shell("echo 2; echo 1", async=True)
	x | shell("sort", async=True) | sys.stdout
	shell.join_all() # ensure that all shells are finished before next line
	
	#import time; time.sleep(1)
	print("test")
	
	y = []
	shell("echo foo; for i in {{1..50000}}; do echo test; done; echo ''; echo bar", async=True) | y
	shell.join_all() # ensure that all shells are finished before next line
	print(y)
	
	shell("echo foo; echo bar")
	
	shell.join_all()
	
