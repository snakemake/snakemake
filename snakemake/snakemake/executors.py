
import os, time, textwrap, stat, shutil, random, string
from functools import partial

from snakemake.shell import shell
from snakemake.logging import logger
from snakemake.utils import format
from snakemake.exceptions import print_exception, get_exception_origin, format_error, RuleException

if os.name == "posix":
	from multiprocessing import Event
	from concurrent.futures import ProcessPoolExecutor
	PoolExecutor = ProcessPoolExecutor
else:
	from threading import Event
	from concurrent.futures import ThreadPoolExecutor
	PoolExecutor = ThreadPoolExecutor
	
class DryrunExecutor:

	def __init__(self, workflow, printreason=False, quiet=False, printshellcmds=False):
		self.workflow = workflow
		self.quiet = quiet
		self.printreason = printreason
		self.printshellcmds = printshellcmds

	def run(self, job, callback = None, error_callback = None):
		self._run(job)
		callback(job)

	def shutdown(self):
		pass
	
	def _run(self, job):
		if job.message:
			logger.info(job.message)
		else:
			# TODO show file types
			desc = ["rule {}:".format(job.rule.name)]
			for name, value in (("input", job.input), 
			                    ("output", job.output), 
			                    ("reason", job.reason if self.printreason else None)):
				if value:
					desc.append(self.format_ruleitem(name, value))
			if self.printshellcmds and job.shellcmd:
				desc.append(job.shellcmd)
			logger.info("\n".join(desc))
		
	@staticmethod
	def format_ruleitem(name, value):
		return "" if not value else "\t{}: {}".format(name, value)

class TouchExecutor(DryrunExecutor):

	def run(self, job, callback = None, error_callback = None):
		super()._run(job)
		try:
			for f in job.expanded_output:
				f.touch()
			time.sleep(0.1)
			callback(job)
		except OSError as ex:
			print_exception(ex, self.workflow.linemaps)
			error_callback()

class CPUExecutor(DryrunExecutor):

	def __init__(self, workflow, cores, printreason=False, quiet=False, printshellcmds=False):
		super().__init__(workflow, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
		self._pool = PoolExecutor(max_workers = cores)

	def run(self, job, callback = None, error_callback = None):
		super()._run(job)
		try:
			for f in job.expanded_output:
				f.remove()
		except OSError as ex:
			print_exception("Could not remove output file {} of dynamic rule {}".format(f, job.rule), self.workflow.linemaps)
		future = self._pool.submit(run_wrapper, job.rule.run_func, job.input, job.output, job.wildcards, job.threads, job.log, self.workflow.linemaps)
		future.add_done_callback(partial(self._callback, job, callback, error_callback))
	
	def shutdown(self):
		self._pool.shutdown()
	
	def _callback(self, job, callback, error_callback, future):
		try:
			ex = future.exception()
			if ex:
				raise ex
			job.check_output()
			job.protect_output()
			# TODO handle temp and protected files
			callback(job)
		except (Exception, BaseException) as ex:
			print_exception(ex, self.workflow.linemaps)
			job.cleanup()
			error_callback()
			
class ClusterExecutor(DryrunExecutor):

	def __init__(self, workflow, cores, submitcmd="qsub", printreason=False, quiet=False, printshellcmds=False):
		if workflow.snakemakepath is None:
			raise ValueError("Cluster executor needs to know the path to the snakemake binary.")
		super().__init__(workflow, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
		self.cores = cores
		self.submitcmd = submitcmd
		self.startedjobs = 0
		self._tmpdir = None
	
	def __del__(self):
		shutil.rmtree(self.tmpdir)
	
	def run(self, job, callback = None, error_callback = None):
		super()._run(job)
		workdir = os.getcwd()
		jobid = self.startedjobs
		
		jobscript = os.path.join(self.tmpdir, "{}.sh".format(jobid))
		jobfinished = os.path.join(self.tmpdir, "{}.jobfinished".format(jobid))
		jobfailed = os.path.join(self.tmpdir, "{}.jobfailed".format(jobid))
		cores = self.cores if self.cores else ""
		with open(jobscript, "w") as f:
			print(format(textwrap.dedent("""
			                                       #!/bin/sh
			                                       #rule: {job}
			                                       #input: {job.input}
			                                       #output: {job.output}
			                                       {self.workflow.snakemakepath} --force -j{cores} --directory {workdir} --nocolor --quiet {job.output} && touch "{jobfinished}" || touch "{jobfailed}"
			                                       exit 0
			                                       """)), file=f)
		print(os.getcwd())
		os.chmod(jobscript, stat.S_IEXEC)
		shell('ls -l {jobscript}; {self.submitcmd} "{jobscript}"')
		self.startedjobs += 1
		threading.Thread(target=self._wait_for_job, args=(job, jobscript, jobfinished, jobfailed)).start()
	
	def _wait_for_job(self, job, jobfinished, jobfailed, jobscript):
		while True:
			if os.path.exists(jobfinished):
				os.remove(jobfinished)
				os.remove(jobscript)
				job.finished()
				# TODO handle temp and protected files
				return
			if os.path.exists(jobfailed):
				os.remove(jobfailed)
				os.remove(jobscript)
				print_exception(ClusterJobException(job), self.workflow.linemaps)
				self._error = True
				self._jobs = set()
				self._open_jobs.set()
				return
			time.sleep(1)
			
	@property
	def tmpdir(self):
		if self._tmpdir is None:		
			while True:
				self._tmpdir = ".snakemake.tmp." + "".join(random.sample(string.ascii_uppercase + string.digits,6))
				if not os.path.exists(self._tmpdir):
					os.mkdir(self._tmpdir)
		return self._tmpdir


def run_wrapper(run, input, output, wildcards, threads, log, linemaps):
	"""
	Wrapper around the run method that handles directory creation and
	output file deletion on error.
	
	Arguments
	run       -- the run method
	input     -- list of input files
	output    -- list of output files
	wildcards -- so far processed wildcards
	threads   -- usable threads
	log       -- path to log file
	"""
	try:
		# execute the actual run method.
		run(input, output, wildcards, threads, log)
		# finish all spawned shells.
		shell.join_all()
	except (Exception, BaseException) as ex:
		# this ensures that exception can be re-raised in the parent thread
		lineno, file = get_exception_origin(ex, linemaps)
		raise RuleException(format_error(ex, lineno, linemaps=linemaps, snakefile=file))
