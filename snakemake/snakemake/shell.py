import _io
import sys, os, inspect
import subprocess as sp
from threading import Thread

def format(string, *args, stepout = 1, **kwargs):
	frame = inspect.currentframe().f_back
	while stepout > 1:
		frame = frame.f_back
		stepout -= 1
	
	variables = dict(frame.f_globals)
	# add local variables from calling rule/function
	variables.update(frame.f_locals)
	variables.update(kwargs)
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
		
class Shell(sp.Popen):
	_process_args = {}
	_processes = []
	
	def __init__(self, cmd, *args, **kwargs):
		super(Shell, self).__init__(format(cmd, *args, stepout = 3, **kwargs), shell=True, stdin = sp.PIPE, stdout=sp.PIPE, close_fds=True, **Shell._process_args)
		self._stdout_free = True
		self._pipethread = None
		self._stdin = self.stdin
		Shell._processes.append(self)
	
	def wait(self):
		super(Shell, self).wait()
		if self._pipethread:
			self._pipethread.join()
			
		
	@staticmethod
	def join():
		for p in Shell._processes:
			if p._stdout_free:
				for l in p._stdoutlines(): pass
			p.wait()
		Shell._processes = []
		
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
		if not isinstance(other, tuple):
			other = tuple([other])

		pipes = []
		toclose = []
		for o in other:
			writer = None
			if isinstance(o, Shell):
				writer = PipeWriter(o._stdin)
				toclose.append(o._stdin)
				# move all stdin to the beginning of the pipe
				o._stdin = self._stdin
			elif isinstance(o, _io._TextIOBase):
				writer = TextIOWriter(o)
			elif isinstance(o, list):
				writer = ListWriter(o)
			else:
				print(o, file=sys.stderr)
				raise ValueError("Only shell, files, stdout or lists allowed right to a shell pipe.")
			pipes.append(writer)
		
		self._stdout_free = False
		self._pipethread = Thread(target = self._write_pipes, args = (pipes, toclose))
		self._pipethread.start()
			
		if len(other) == 1:
			return other[0]

if "SHELL" in os.environ:
	Shell._process_args["executable"] = os.environ["SHELL"]

def shell(cmd, *args, **kwargs):
	sp.check_call(format(cmd, *args, stepout=2, **kwargs), shell=True, **Shell._process_args)
	#p = Shell(cmd, *args, **kwargs)
	#p | sys.stdout
	#p.wait()

if __name__ == "__main__":
	Shell("echo b; echo a; echo c > foo")
	Shell("echo b; echo a; echo c") | Shell("sort") | sys.stdout
	
	Shell("echo a; echo b") | (Shell("sort"), Shell("cut -f1"))
	
	x = Shell("echo 2; echo 1")
	x | Shell("sort") | sys.stdout
	Shell.join() # ensure that all shells are finished before next line
	
	#import time; time.sleep(1)
	print("test")
	
	y = []
	Shell("echo foo; echo ''; echo bar") | y
	Shell.join() # ensure that all shells are finished before next line
	#import time; time.sleep(5)
	print(y)
	
	Shell.join()
	
