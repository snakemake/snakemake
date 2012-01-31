import inspect
import sys
import subprocess as sp

class shell(sp.Popen):
	processes = []
	def __init__(self, cmd):
		super(shell, self).__init__(cmd, shell=True, stdin = sp.PIPE, stdout=sp.PIPE, close_fds=True)
		shell.processes.append(self)
		
	def stdoutlines(self):
		while True:
			o = self.stdout.readline()
			if o == b'' and self.poll() != None:
				return
			yield o
				
	def __or__(self, other):
		if isinstance(other, shell):
			for o in self.stdoutlines():
				other.stdin.write(o)
			other.stdin.close()
		else:
			for o in self.stdoutlines():
				other.write(o.decode("utf-8"))
		return other


shell("echo b; echo a; echo c > blub")
shell("echo b; echo a; echo c") | shell("sort") | sys.stdout

x = shell("echo 2; echo 1")
x | shell("sort") | sys.stdout

for p in shell.processes:
	p.wait()
