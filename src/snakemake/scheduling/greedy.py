from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Sequence, Union
from snakemake_interface_scheduler_plugins.base import SchedulerBase
from snakemake_interface_scheduler_plugins.settings import SchedulerSettingsBase
from snakemake_interface_scheduler_plugins.interfaces.jobs import JobSchedulerInterface
from snakemake_interface_common.io import AnnotatedStringInterface


@dataclass
class SchedulerSettings(SchedulerSettingsBase):
    greediness: Optional[float] = field(
        default=1.0,
        metadata={
            "help": "Set the greediness of scheduling. This value between 0 and 1 "
            "determines how careful jobs are selected for execution. The default "
            "value (1.0) provides the best speed and still acceptable scheduling "
            "quality.",
        },
    )
    omit_prioritize_by_temp_and_input: bool = field(
        default=False,
        metadata={
            "help": "If set, jobs with larger temporary or input files are "
            "not prioritized. The rationale of the prioritization is that temp files "
            "should be removed as soon as possible, and larger input files may take "
            "longer to process, so it is better to start them earlier.",
        },
    )

    def __post_init__(self):
        if self.greediness is None:
            self.greediness = 1.0
        if not (0 <= self.greediness <= 1.0):
            raise ValueError("greediness must be >=0 and <=1")


class Scheduler(SchedulerBase):

    def __post_init__(self):
        self._input_sizes = {}

    def select_jobs(
        self,
        selectable_jobs: Sequence[JobSchedulerInterface],
        remaining_jobs: Sequence[JobSchedulerInterface],
        available_resources: Mapping[str, Union[int, str]],
        input_sizes: Dict[AnnotatedStringInterface, int],
    ) -> Sequence[JobSchedulerInterface]:
        """
        Using the greedy heuristic from
        "A Greedy Algorithm for the General Multidimensional Knapsack
        Problem", Akcay, Li, Xu, Annals of Operations Research, 2012

        Args:
            jobs (list):    list of jobs
        """
        # each job is an item with one copy (0-1 MDKP)
        n = len(selectable_jobs)
        x = [0] * n  # selected jobs
        E = set(range(n))  # jobs still free to select
        u = [1] * n
        a = [
            self.job_weight(job, available_resources) for job in selectable_jobs
        ]  # resource usage of jobs

        c = [self.job_reward(job, input_sizes) for job in selectable_jobs]

        def calc_reward():
            return [c_j * y_j for c_j, y_j in zip(c, y)]

        b = [
            available_resources[name] for name in available_resources
        ]  # resource capacities

        while True:
            # Step 2: compute effective capacities
            y = [
                (
                    min(
                        (min(u[j], b_i // a_j_i) if a_j_i > 0 else u[j])
                        for b_i, a_j_i in zip(b, a[j])
                        if a_j_i
                    )
                    if j in E
                    else 0
                )
                for j in range(n)
            ]
            if not any(y):
                break
            y = [
                (max(1, int(self.settings.greediness * y_j)) if y_j > 0 else 0)
                for y_j in y
            ]

            # Step 3: compute rewards on cumulative sums
            reward = calc_reward()
            j_sel = max(E, key=reward.__getitem__)  # argmax

            # Step 4: batch increment
            y_sel = y[j_sel]

            # Step 5: update information
            x[j_sel] += y_sel
            b = [b_i - (a_j_i * y_sel) for b_i, a_j_i in zip(b, a[j_sel])]
            u[j_sel] -= y_sel
            if not u[j_sel] or self.settings.greediness == 1:
                E.remove(j_sel)
            if not E:
                break

        solution = [job for job, sel in zip(selectable_jobs, x) if sel]
        return solution

    def job_weight(self, job, available_resources: Mapping[str, Union[int, str]]):
        res = job.scheduler_resources
        return [res.get(name, 0) for name in available_resources]

    def job_reward(self, job, input_sizes: Dict[AnnotatedStringInterface, int]):
        if self.settings.omit_prioritize_by_temp_and_input:
            return job.priority
        else:
            # Usually, this should guide the scheduler to first schedule all jobs
            # that remove the largest temp file, then the second largest and so on.
            # Since the weight is summed up, it can in theory be that it sometimes
            # prefers a set of many jobs that all depend on smaller temp files though.
            # A real solution to the problem is therefore to use dummy jobs that
            # ensure selection of groups of jobs that together delete the same temp
            # file.
            return (
                job.priority,
                sum(input_sizes[f] or 0 for f in job.input if f.is_flagged("temp")),
                sum(input_sizes[f] or 0 for f in job.input),
            )
