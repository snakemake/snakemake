import inspect
import sys
import subprocess as sp

class shell:
	def __init__(self, cmd):
		self.cmd = cmd
		self._process = None
	
	def _run(self, stdin = None):
		if not self._process:
			self._process = sp.Popen(self.cmd, shell=True, stdin = stdin, stdout=sp.PIPE, close_fds=True)

	def wait(self):
		self._run()
		self._process.wait()

	def __or__(self, other):
		self._run()
		if isinstance(other, shell):
			other._run(stdin = self._process.stdout)
		else:
			stdout = self._process.stdout.read().decode("utf-8")
			if isinstance(other, str):
				other += stdout
			else:
				other.write(stdout)
		return other


#shell("echo b; echo a; echo c > blub")
shell("echo b; echo a; echo c") | shell("sort") | sys.stdout
