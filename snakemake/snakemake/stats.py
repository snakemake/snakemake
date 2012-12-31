
import time, csv
from collections import defaultdict

class Stats:
	def __init__(self):
		self.starttime = dict()
		self.endtime = dict()
	
	def report_job_start(self, job):
		self.starttime[job] = time.time()
	
	def report_job_end(self, job):
		self.endtime[job] = time.time()
	
	@property
	def rule_runtimes(self):
		runtimes = defaultdict(list)
		for job, t in self.starttime.items():
			runtimes[job.rule].append(self.endtime[job] - t)
		for rule, runtimes in runtimes.items():
			yield rule, sum(runtimes) / len(runtimes), min(runtimes), max(runtimes)
	
	@property
	def job_runtimes(self):
		for job, t in self.starttime.items():
			yield job, t, self.endtime[job]
	
	@property
	def overall_runtime(self):
		if self.starttime and self.endtime:
			return max(self.endtime.values()) - min(self.starttime.values())
		else:
			return 0

	def to_csv(self, path):
		with open(path, "w") as f:
			writer = csv.writer(f, delimiter="\t")
			writer.writerow("overall-runtime".split())
			writer.writerow([self.overall_runtime])
			writer.writerow("rule mean-runtime min-runtime max-runtime".split())
			for runtime in self.rule_runtimes:
				writer.writerow(runtime)
			writer.writerow(list())
			writer.writerow("file starttime endtime".split())
			for runtime in self.job_runtimes:
				writer.writerow(runtime)
