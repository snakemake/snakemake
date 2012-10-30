


class KnapsackJobScheduler:
	def __init__(self, jobs, workflow):
		""" Create a new instance of KnapsackJobScheduler. """
		self.workflow = workflow
		self._maxcores = workflow.cores if workflow.cores else multiprocessing.cpu_count()
		self._cores = self._maxcores
		self._pool = PoolExecutor(max_workers = self._cores)
		self._jobs = set()
		self.add_jobs(jobs)
		self._open_jobs = Event()
		self._open_jobs.set()
		self._errors = False

	def get_jobs(self):
		return self._jobs

	def add_jobs(self, jobs):
		for job in jobs:
			job.scheduler = self
			self._jobs.add(job)

	def remove_job(self, job):
		self._jobs.remove(job)

	def schedule(self):
		""" Schedule jobs that are ready, maximizing cpu usage. """
		while True:
			self._open_jobs.wait()
			self._open_jobs.clear()
			if self._errors:
				logger.warning("Will exit after finishing currently running jobs.")
				self._pool.shutdown()
				return False
			if not self._jobs:
				self._pool.shutdown()
				return True

			needrun, norun = [], set()
			for job in self._jobs:
				if job.depends:
					continue
				if job.needrun:
					if job.threads > self._maxcores:
						# reduce the number of threads so that it 
						# fits to available cores.
						if not job.dryrun:
							logger.warn(
								"Rule {} defines too many threads ({}), Scaling down to {}."
								.format(job.rule, job.threads, self._maxcores))
						job.threads = self._maxcores
					needrun.append(job)
				else: norun.add(job)

			
			run = self._knapsack(needrun)
			self._jobs -= run
			self._jobs -= norun
			self._cores -= sum(job.threads for job in run)
			for job in chain(run, norun):
				job.add_callback(self._finished)
				job.add_error_callback(self._error)
				job.run(self._run_job)
			
	
	def _run_job(self, job):
		future = self._pool.submit(run_wrapper, *job.get_run_args())
		future.add_done_callback(job.finished)
		
	def _finished(self, job):
		if job.needrun:
			self._cores += job.threads
		self._open_jobs.set()
	
	def _error(self):
		# clear jobs and stop the workflow
		self._errors = True
		self._jobs = set()
		self._open_jobs.set()
	
	def _knapsack(self, jobs):
		""" Solve 0-1 knapsack to maximize cpu utilization. """
		dimi, dimj = len(jobs) + 1, self._cores + 1
		K = [[0 for c in range(dimj)] for i in range(dimi)]
		for i in range(1, dimi):
			for j in range(1, dimj):
				t = jobs[i-1].threads
				if t > j:
					K[i][j] = K[i - 1][j]
				else:
					K[i][j] = max(K[i - 1][j], t + K[i - 1][j - t])
		
		solution = set()
		i = dimi - 1
		j = dimj - 1
		while i > 0:
			if K[i][j] != K[i-1][j]:
				job = jobs[i - 1]
				solution.add(job)
				j = j - job.threads
			i -= 1
		
		return solution

class ClusterJobScheduler:
	def __init__(self, jobs, workflow, submitcmd = "qsub"):
		self.workflow = workflow
		self._jobs = set()
		self.add_jobs(jobs)
		self._submitcmd = submitcmd
		self._open_jobs = Event()
		self._open_jobs.set()
		self._error = False
		self._cores = workflow.cores

	def add_jobs(self, jobs):
		for job in jobs:
			job.scheduler = self
			self._jobs.add(job)

	def schedule(self):
		while True:
			self._open_jobs.wait()
			self._open_jobs.clear()
			if self._error:
				logger.warning("Will exit after finishing currently running jobs.")
				return False
			if not self._jobs:
				return True
			needrun, norun = set(), set()
			for job in self._jobs:
				if job.depends:
					continue
				if job.needrun:
					needrun.add(job)
				else: norun.add(job)

			self._jobs -= needrun
			self._jobs -= norun
			for job in chain(needrun, norun):
				job.add_callback(self._finished)
				job.run(self._run_job)
	
	def _run_job(self, job):
		job.print_message()
		workdir = os.getcwd()
		prefix = ".snakemake.{}.".format(job.rule.name)
		jobid = "_".join(job.output).replace("/", "_")
		jobscript = "{}.{}.sh".format(prefix, jobid)
		jobfinished = "{}.{}.jobfinished".format(prefix, jobid)
		jobfailed = "{}.{}.jobfailed".format(prefix, jobid)
		cores = self._cores if self._cores else ""
		scriptpath = self.workflow.scriptpath
		if not scriptpath:
			scriptpath = "snakemake"
		shell("""
			echo '#!/bin/sh' > "{jobscript}"
			echo '#rule: {job}' >> "{jobscript}"
			echo '#input: {job.input}' >> "{jobscript}"
			echo '#output: {job.output}' >> "{jobscript}"
			echo '{scriptpath} --force -j{self._cores} --directory {workdir} --nocolor --quiet {job.output} && touch "{jobfinished}" || touch "{jobfailed}"' >> "{jobscript}"
			chmod +x "{jobscript}"
			{self._submitcmd} "{jobscript}"
		""")
		threading.Thread(target=self._wait_for_job, args=(job, jobfinished, jobfailed, jobscript)).start()

	def _finished(self, job):
		self._open_jobs.set()
		
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
				print_exception(ClusterJobException(job), self.workflow.rowmaps)
				self._error = True
				self._jobs = set()
				self._open_jobs.set()
				return
			time.sleep(1)

