
import os

from snakemake.logging import logger

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
	
	def _run(self, job):
		assert job.needrun
		if job.message:
			logger.info(job.message)
		else:
			# TODO show file types
			items = [job.input, job.output]
			if self.printreason:
				items.append(job.reason)
			desc = ["rule {}:".format(job.rule.name), *map(self.format_ruleitem, items)]
			if self.printshellcmds and job.shellcmd:
				desc.append(job.shellcmd)
			logger.info("\n".join(desc))
		
	@staticmethod
	def format_ruleitem(name, item):
		return "" if not item else "{}: {}".format(name, item)

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
		future = self._pool.submit(run_wrapper, *job.get_run_args())
		future.add_done_callback(partial(self._callback, job, callback, error_callback))
	
	def _callback(self, job, callback, error_callback, future):
		try:
			ex = future.exception()
			if ex:
				raise ex
			for f in job.expanded_output:
				f.created()
			callback()
		except (Exception, BaseException) as ex:
			print_exception(ex, self.workflow.linemaps)
			self.cleanup()
			error_callback()
			
class ClusterExecutor(DryrunExecutor):

	def __init__(self, workflow, cores, submitcmd, printreason=False, quiet=False, printshellcmds=False):
		if workflow.snakemakepath is None:
			raise ValueError("Cluster executor needs to know the path to the snakemake binary.")
		super().__init__(workflow, printreason=printreason, quiet=quiet, printshellcmds=printshellcmds)
		self.cores = cores
		self.submitcmd = submitcmd
		self.tmpdir = tempfile.mkdtemp(prefix="snakemake")
		self.startedjobs = 0
	
	def __del__(self):
		os.remove(self.tmpdir)
	
	def run(self, job, callback = None, error_callback = None):
		super()._run(job)
		workdir = os.getcwd()
		jobid = self.startedjobs
		
		jobscript = os.path.join(self.tmpdir, "{}.sh".format(jobid))
		jobfinished = os.path.join(self.tmpdir, "{}.jobfinished".format(jobid))
		jobfailed = os.path.join(self.tmpdir, "{}.jobfailed".format(jobid))
		cores = self._cores if self._cores else ""
		with open(jobscript, "w") as f:
			print(format_wildcards(textwrap.dedent("""
			                                       #!/bin/sh
			                                       #rule: {job}
			                                       #input: {job.input}
			                                       #output: {job.output}
			                                       {self.workflow.snakemakepath} --force -j{self._cores} --directory {workdir} --nocolor --quiet {job.output} && touch "{jobfinished}" || touch "{jobfailed}"
			                                       """)), file=f)
		os.chmod(fpath, stat.S_IEXEC)
		shell('{self.submitcmd} "{scriptpath}"')
		self.startedjobs += 1
		threading.Thread(target=self._wait_for_job, args=(job, jobscript, jobfinished, jobfailed)).start()
	
	def _wait_for_job(self, job, jobfinished, jobfailed, jobscript):
		while True:
			if os.path.exists(jobfinished):
				os.remove(jobfinished)
				os.remove(jobscript)
				job.finished()
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


def run_wrapper(run, rulename, ruledesc, input, output, wildcards, 
		threads, log, rowmaps, rulelineno, rulesnakefile):
	"""
	Wrapper around the run method that handles directory creation and
	output file deletion on error.
	
	Arguments
	run -- the run method
	input -- list of input files
	output -- list of output files
	wildcards -- so far processed wildcards
	"""
	logger.info(ruledesc)

			
	t0 = time.time()
	try:
		# execute the actual run method.
		run(input, output, wildcards, threads, log)
		# finish all spawned shells.
		shell.join_all()
		t1 = time.time()
		return t0, t1
	except (Exception, BaseException) as ex:
		# this ensures that exception can be re-raised in the parent thread
		lineno, file = get_exception_origin(ex, rowmaps)
		raise RuleException(format_error(ex, lineno, rowmaps=rowmaps, snakefile=file), )
