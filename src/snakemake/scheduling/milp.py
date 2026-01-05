from dataclasses import dataclass, field
import math
import os
from functools import cached_property
from pathlib import Path
from typing import Collection, Dict, Iterator, List, Mapping, Optional, Sequence, Union
from snakemake_interface_scheduler_plugins.base import SchedulerBase
from snakemake_interface_scheduler_plugins.settings import SchedulerSettingsBase
from snakemake_interface_scheduler_plugins.interfaces.jobs import JobSchedulerInterface
from snakemake_interface_common.io import AnnotatedStringInterface


class LpSolverCollection(Collection[str]):
    """
    A lazy collection that avoids calling pulp.listSolvers if the default solver is selected
    """

    def __init__(self):
        """
        Initialize the collection and determine the default PuLP solver name if available.
        
        Sets self.default to the name of pulp.apis.LpSolverDefault when PuLP is importable and a default is defined; otherwise sets self.default to None. Import errors are suppressed.
        """
        default = None
        try:
            import pulp

            solver_default = pulp.apis.LpSolverDefault
            if solver_default is not None:
                default = solver_default.name
        except ImportError:
            pass

        self.default = default

    @cached_property
    def nondefault_solvers(self) -> List[str]:
        """
        Return a sorted list of available MILP solver names excluding the configured default.
        
        If the pulp package is not installed, returns an empty list.
        
        Returns:
            List[str]: Sorted solver names available on the system, excluding self.default.
        """
        try:
            import pulp

            return sorted(
                solver
                for solver in pulp.listSolvers(onlyAvailable=True)
                if solver != self.default
            )
        except ImportError:
            return []

    def __iter__(self) -> Iterator[str]:
        """
        Iterate over available solver names, yielding the default first if present.
        
        Returns:
            Iterator[str]: Solver names with the default solver (if set) yielded first, followed by non-default solvers.
        """
        if self.default is not None:
            yield self.default
        yield from self.nondefault_solvers

    def __contains__(self, x: object) -> bool:
        """
        Check whether the given solver identifier is available to pulp.
        
        Parameters:
            x (object): Solver name or solver specification passed to pulp.getSolver.
        
        Returns:
            `true` if pulp reports the solver is available, `false` otherwise (including when pulp is unavailable or the lookup fails).
        """
        try:
            import pulp

            return pulp.getSolver(x).available()
        except Exception:  # noqa: BLE001
            return False

    def __len__(self) -> int:
        """
        Compute the number of available MILP solvers, counting the default solver if present.
        
        Returns:
            total (int): Total number of solvers; includes 1 if a default solver exists plus the number of nondefault solvers.
        """
        return (1 if self.default is not None else 0) + len(self.nondefault_solvers)


lp_solvers = LpSolverCollection()


@dataclass
class SchedulerSettings(SchedulerSettingsBase):
    solver: Optional[str] = field(
        default=lp_solvers.default,
        metadata={
            "help": "Set MILP solver to use",
            "choices": lp_solvers,
        },
    )
    solver_path: Optional[Path] = field(
        default=None,
        metadata={"help": "Set the PATH to search for scheduler solver binaries."},
    )

    @property
    def lp_solver_available(self):
        """
        Indicates whether the configured MILP solver is available.
        
        Returns:
            True if the configured solver name is present in lp_solvers, False otherwise.
        """
        return self.solver in lp_solvers


class Scheduler(SchedulerBase):
    def __post_init__(self) -> None:
        """
        Initialize post-construction state for the scheduler.
        
        Set the internal `_technical_failure` flag to False to indicate no solver-related technical failure has occurred yet.
        """
        self._technical_failure = False

    def select_jobs(
        self,
        selectable_jobs: Sequence[JobSchedulerInterface],
        remaining_jobs: Sequence[JobSchedulerInterface],
        available_resources: Mapping[str, Union[int, str]],
        input_sizes: Dict[AnnotatedStringInterface, int],
    ) -> Optional[Sequence[JobSchedulerInterface]]:
        if self._technical_failure:
            # fallback early since we failed before already
            return None
        import pulp
        from pulp import lpSum
        from pulp import PulpSolverError

        scheduled_jobs = {
            job: pulp.LpVariable(
                f"job_{idx}", lowBound=0, upBound=1, cat=pulp.LpInteger
            )
            for idx, job in enumerate(selectable_jobs)
        }

        job_temp_files = {}
        for job in remaining_jobs:
            job_temp_files[job] = {
                infile for infile in job.input if infile.is_flagged("temp")
            }

        temp_files = {
            f for job in selectable_jobs for f in job.input if f.is_flagged("temp")
        }

        temp_sizes_gb = {f: input_sizes[f] / 1e9 for f in temp_files}

        temp_job_improvement = {
            temp_file: pulp.LpVariable(
                f"temp_file_{idx}", lowBound=0, upBound=1, cat="Continuous"
            )
            for idx, temp_file in enumerate(temp_files)
        }

        temp_file_deletable = {
            temp_file: pulp.LpVariable(
                f"deletable_{idx}",
                lowBound=0,
                upBound=1,
                cat=pulp.LpInteger,
            )
            for idx, temp_file in enumerate(temp_files)
        }
        prob = pulp.LpProblem("JobScheduler", pulp.LpMaximize)

        total_temp_size = max(
            sum([temp_sizes_gb[temp_file] for temp_file in temp_files]), 1
        )
        total_core_requirement = sum(
            [
                max(job.scheduler_resources.get("_cores", 1), 1)
                for job in selectable_jobs
            ]
        )
        # Objective function
        # Job priority > Core load
        # Core load > temp file removal
        # Instant removal > temp size
        prob += (
            2
            * total_core_requirement
            * 2
            * total_temp_size
            * lpSum([job.priority * scheduled_jobs[job] for job in selectable_jobs])
            + 2
            * total_temp_size
            * lpSum(
                [
                    max(job.scheduler_resources.get("_cores", 1), 1)
                    * scheduled_jobs[job]
                    for job in selectable_jobs
                ]
            )
            + total_temp_size
            * lpSum(
                [
                    temp_file_deletable[temp_file] * temp_sizes_gb[temp_file]
                    for temp_file in temp_files
                ]
            )
            + lpSum(
                [
                    temp_job_improvement[temp_file] * temp_sizes_gb[temp_file]
                    for temp_file in temp_files
                ]
            )
        )

        # Constraints:
        for name in available_resources:
            prob += (
                lpSum(
                    [
                        scheduled_jobs[job] * job.scheduler_resources.get(name, 0)
                        for job in selectable_jobs
                    ]
                )
                <= available_resources[name]
            )

        # Choose jobs that lead to "fastest" (minimum steps) removal of existing temp file
        for temp_file in temp_files:
            prob += temp_job_improvement[temp_file] <= lpSum(
                [
                    scheduled_jobs[job]
                    for job in selectable_jobs
                    if temp_file in job_temp_files[job]
                ]
            ) / lpSum([1 for job in remaining_jobs if temp_file in job_temp_files[job]])

            prob += temp_file_deletable[temp_file] <= temp_job_improvement[temp_file]

        try:
            status = self._solve_ilp(prob, time_limit=10)
        except PulpSolverError as e:
            self._technical_failure = True
            self.logger.warning(
                "Failed to solve scheduling problem with ILP solver, falling back to "
                "greedy scheduler. You likely have to fix your ILP solver "
                f"installation. Error message: {e}"
            )
            return None
        if pulp.LpStatus[status] != "Optimal":
            if pulp.LpStatus[status] == "Not Solved":
                self.logger.warning(
                    "Failed to solve scheduling problem with ILP solver in time (10s), "
                    "falling back to greedy scheduler."
                )
            elif pulp.LpStatus[status] == "Infeasible":
                self.logger.warning(
                    "Failed to solve scheduling problem with ILP solver, falling back "
                    "to greedy scheduler."
                )
            return None

        selected_jobs = [
            job
            for job, variable in scheduled_jobs.items()
            if math.isclose(variable.value(), 1.0)
        ]

        if not selected_jobs:
            # No selected jobs. This could be due to insufficient resources or a failure in the ILP solver
            # Hence, we silently fall back to the greedy solver to make sure that we don't miss anything.
            return None

        return selected_jobs

    def _solve_ilp(self, prob, threads=2, time_limit=10):
        """
        Solve the given ILP problem using the configured PuLP solver and return the solver status.
        
        Parameters:
            prob: The PuLP LpProblem to solve.
            threads (int): Number of threads to request from the solver.
            time_limit (int): Time limit (in seconds) to pass to the solver.
        
        Notes:
            If SchedulerSettings.solver_path is set, this function temporarily prepends that path to the process PATH so the solver binary can be found and restores the original PATH before solving.
        
        Returns:
            status (str): Solver result status such as "Optimal", "Infeasible", or "Not Solved".
        """
        import pulp

        old_path = os.environ["PATH"]
        if self.settings.solver_path is not None:
            # Temporarily prepend the given snakemake env to the path, such that the solver can be found in any case.
            # This is needed for cluster envs, where the cluster job might have a different environment but
            # still needs access to the solver binary.
            os.environ["PATH"] = "{}:{}".format(
                self.settings.solver_path,
                os.environ["PATH"],
            )
        try:
            solver = pulp.getSolver(self.settings.solver)
        finally:
            os.environ["PATH"] = old_path
        solver.optionsDict["threads"] = threads
        solver.timeLimit = time_limit
        solver.msg = False  # Suppress solver output
        return prob.solve(solver)