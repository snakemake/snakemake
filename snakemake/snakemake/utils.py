import os, subprocess, inspect
if "SHELL" in os.environ:
	def _shell(cmd):
		subprocess.check_call(cmd, shell=True, executable = os.environ["SHELL"])
else:
	def _shell(cmd):
		subprocess.check_call(cmd, shell=True)

def shell(cmd, *args, **kwargs):
	variables = dict(inspect.currentframe().f_back.f_globals)
	# add local variables from calling rule/function
	variables.update(inspect.currentframe().f_back.f_locals)
	variables.update(kwargs)
	try:
		_shell(cmd.format(*args, **variables))
	except KeyError as ex:
		raise NameError("The name {} is unknown in this context.".format(str(ex)))
