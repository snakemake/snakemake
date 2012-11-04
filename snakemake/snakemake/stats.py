
import time

class Stats:
	def __init__(self):
		self.starttime = dict()
		self.endtime = dict()
	
	def report_job_start(self, job):
		self.starttime[job] = time.time()
	
	def report_job_end(self, job):
		self.endtime[job] = time.time()
	
	def rule_runtimes(self):
		runtimes = defaultdict(list)
		for job, t in self.starttime.items():
			runtimes[job.rule].append(self.endtime[job] - t)
		for rule, runtimes in runtimes:
			yield rule, sum(runtimes) / len(runtimes), min(runtimes), max(runtimes)
	
	def job_runtimes(self):
		for job, t in self.starttime.items():
			yield job, t, self.endtime[job]
