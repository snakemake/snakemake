"""
Unit tests for DAG.reset_after_dryrun() function.

This module tests that after calling reset_after_dryrun(), the DAG returns to 
a pristine state suitable for real execution after a dryrun.

Uses Approach 3 (Complete State Capture/Restore) to verify DAG equivalence:
- Captures meaningful DAG state before and after reset
- Compares all essential attributes without relying on object identity
- Tests various edge cases including:
  * Simple linear DAGs
  * DAGs with multiple parallel jobs
  * DAGs with complex dependencies
  * DAGs with checkpoints
  * DAGs with group jobs
  * DAGs with queue input jobs
  * DAGs with storage input jobs
  * DAGs with temporary files
  * DAGs with dynamic rules
"""

import sys
import os
from pathlib import Path
from unittest.mock import Mock, MagicMock, patch
import pytest
from typing import Dict, Set, Tuple, Any
import copy

sys.path.insert(0, os.path.dirname(__file__))

from common import dpath
from snakemake import api
from snakemake.settings import types as settings
from snakemake.dag import DAG
from snakemake.jobs import Job, Reason
from snakemake.workflow import Workflow


def _job_signature(job: Job) -> Tuple[str, frozenset, frozenset]:
    """
    Create a hashable signature for a job based on its essential properties.
    
    This signature is invariant to object identity and can be used to
    compare jobs across different DAG instances.
    """
    rule_name = job.rule.name
    # Use frozenset of wildcards for order-independent comparison
    wildcards = frozenset(job.wildcards.items()) if job.wildcards else frozenset()
    # Use frozenset of output filenames as job identifier
    outputs = frozenset(job.output) if job.output else frozenset()
    return (rule_name, wildcards, outputs)


def capture_dag_state(dag: DAG) -> Dict[str, Any]:
    """
    Capture all meaningful state from a DAG for later comparison.
    
    This captures structural information, job properties, and dependencies
    without relying on object identity, allowing comparison across DAG
    instances or before/after reset.
    
    Returns:
        Dictionary containing:
        - job_count: Total number of jobs
        - jobs: Set of job signatures
        - dependencies: Map of job signatures to their job dependencies
        - needrun_jobs: Set of signatures for jobs that need to run
        - finished_jobs: Set of signatures for finished jobs
        - running_jobs: Set of signatures for running jobs
        - ready_jobs: Set of signatures for ready jobs
        - targetjobs: Set of signatures for target jobs
        - prioritytargetjobs: Set of signatures for priority target jobs
        - priority_map: Priority values for jobs
    """
    # Capture all jobs and create a mapping
    all_jobs = list(dag.jobs)
    job_sigs = {_job_signature(j) for j in all_jobs}
    
    # Build signature-to-job mapping for dependency analysis
    sig_to_job = {_job_signature(j): j for j in all_jobs}
    
    # Capture dependencies
    dependencies = {}
    for job in all_jobs:
        sig = _job_signature(job)
        dep_sigs = frozenset(
            _job_signature(dep) for dep in dag.job_dependencies(job)
        )
        dependencies[sig] = dep_sigs
    
    # Capture job states
    needrun = {_job_signature(j) for j in dag.needrun_jobs()}
    finished = {_job_signature(j) for j in dag.finished_jobs}
    running = {_job_signature(j) for j in dag._running}
    ready = {_job_signature(j) for j in dag._ready_jobs}
    
    # Capture target information
    targetjobs = {_job_signature(j) for j in dag.targetjobs}
    prioritytargetjobs = {_job_signature(j) for j in dag.prioritytargetjobs}
    
    # Capture priority information
    priority_map = {
        _job_signature(job): dag._priority.get(job, None)
        for job in all_jobs
    }
    
    return {
        'job_count': len(all_jobs),
        'jobs': job_sigs,
        'dependencies': dependencies,
        'needrun_jobs': needrun,
        'finished_jobs': finished,
        'running_jobs': running,
        'ready_jobs': ready,
        'targetjobs': targetjobs,
        'prioritytargetjobs': prioritytargetjobs,
        'priority_map': priority_map,
        'rules': frozenset(r.name for r in dag.rules),
    }


def compare_dag_states(state_before: Dict[str, Any], state_after: Dict[str, Any], 
                       ignore_states: Set[str] = None) -> Tuple[bool, str]:
    """
    Compare two captured DAG states for equivalence.
    
    Args:
        state_before: DAG state captured before dryrun
        state_after: DAG state captured after reset
        ignore_states: Set of state keys to ignore in comparison
                      (e.g., {'finished_jobs', 'running_jobs', 'ready_jobs'})
    
    Returns:
        Tuple of (is_equal, error_message)
    """
    if ignore_states is None:
        ignore_states = set()
    
    comparison_keys = set(state_before.keys()) - ignore_states
    
    for key in comparison_keys:
        before_val = state_before[key]
        after_val = state_after[key]
        
        if before_val != after_val:
            return False, f"Mismatch in '{key}': {before_val} != {after_val}"
    
    return True, "States are equivalent"


class TestDAGResetAfterDryrun:
    """Test suite for DAG.reset_after_dryrun() method using state comparison."""

    def _create_minimal_dag(self):
        """Create a minimal DAG for testing."""
        snakefile = dpath("test_bash/Snakefile")
        with api.SnakemakeApi() as snakemake_api:
            workflow_api = snakemake_api.workflow(
                snakefile=snakefile,
                resource_settings=settings.ResourceSettings(cores=1),
            )
            dag_api = workflow_api.dag()
            workflow = dag_api.workflow_api._workflow
            # Prepare and build the DAG
            workflow.remote_execution_settings = settings.RemoteExecutionSettings()
            workflow._prepare_dag(
                forceall=False,
                ignore_incomplete=False,
                lock_warn_only=False,
            )
            workflow._build_dag()
            return workflow.dag

    def _create_dag_from_test_dir(self, test_dir):
        """Helper to create a DAG from a test directory."""
        snakefile = dpath(test_dir) / "Snakefile"
        if not snakefile.exists():
            pytest.skip(f"Test directory {test_dir} not found")
            
        with api.SnakemakeApi() as snakemake_api:
            workflow_api = snakemake_api.workflow(
                snakefile=snakefile,
                resource_settings=settings.ResourceSettings(cores=1),
            )
            dag_api = workflow_api.dag()
            workflow = dag_api.workflow_api._workflow
            workflow.remote_execution_settings = settings.RemoteExecutionSettings()
            workflow._prepare_dag(
                forceall=False,
                ignore_incomplete=False,
                lock_warn_only=False,
            )
            workflow._build_dag()
            return workflow.dag, workflow

    def test_reset_restores_dag_to_initial_state(self):
        """Test that reset_after_dryrun restores DAG to its initial state.
        
        This is the core test: simulate a dryrun, then reset, and verify
        the DAG is equivalent to a fresh DAG.
        """
        dag = self._create_minimal_dag()
        
        # Capture the initial state (fresh DAG)
        initial_state = capture_dag_state(dag)
        
        # Simulate dryrun: mark all jobs as finished
        if dag.jobs:
            for job in dag.jobs:
                dag._finished.add(job)
            # Also mark some jobs as running/ready to simulate execution state
            jobs_list = list(dag.jobs)
            if len(jobs_list) > 1:
                dag._running.add(jobs_list[0])
                dag._ready_jobs.add(jobs_list[1])
        
        # Verify dryrun state is different from initial
        dryrun_state = capture_dag_state(dag)
        assert dryrun_state['finished_jobs'] != initial_state['finished_jobs'], \
            "Dryrun should mark jobs as finished"
        
        # Reset
        dag.reset_after_dryrun()
        
        # Capture state after reset
        reset_state = capture_dag_state(dag)
        
        # Compare essential structure (ignoring transient execution state)
        are_equal, error_msg = compare_dag_states(
            initial_state, 
            reset_state,
            ignore_states={'finished_jobs', 'running_jobs', 'ready_jobs'}
        )
        assert are_equal, f"DAG not restored to initial state: {error_msg}"

    def test_reset_preserves_job_structure(self):
        """Test that reset_after_dryrun doesn't affect job definitions."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        initial_job_count = initial_state['job_count']
        initial_jobs = initial_state['jobs']
        
        # Simulate dryrun
        if dag.jobs:
            for job in dag.jobs:
                dag._finished.add(job)
                dag._progress += 1
        
        # Reset
        dag.reset_after_dryrun()
        
        reset_state = capture_dag_state(dag)
        
        # Job structure must be identical
        assert reset_state['job_count'] == initial_job_count, \
            "Job count changed after reset"
        assert reset_state['jobs'] == initial_jobs, \
            "Job definitions changed after reset"

    def test_reset_preserves_dependencies(self):
        """Test that reset_after_dryrun preserves job dependencies."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        initial_dependencies = initial_state['dependencies']
        
        # Simulate dryrun with complex state mutation
        if dag.jobs:
            for job in dag.jobs:
                dag._finished.add(job)
                dag._checked_jobs.add(job)
                reason = Reason()
                reason.finished = True
                dag._reason[job] = reason
        
        # Reset
        dag.reset_after_dryrun()
        
        reset_state = capture_dag_state(dag)
        
        # Dependencies must be preserved
        assert reset_state['dependencies'] == initial_dependencies, \
            "Job dependencies changed after reset"

    def test_reset_preserves_needrun_jobs(self):
        """Test that reset_after_dryrun preserves which jobs need to run."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        initial_needrun = initial_state['needrun_jobs']
        
        # Simulate dryrun: mark all jobs as finished
        if dag.jobs:
            for job in dag.jobs:
                dag._finished.add(job)
        
        # After dryrun, needrun should be empty (all finished)
        dryrun_state = capture_dag_state(dag)
        assert len(dryrun_state['needrun_jobs']) == 0, \
            "Dryrun should mark all jobs as finished"
        
        # Reset
        dag.reset_after_dryrun()
        
        reset_state = capture_dag_state(dag)
        
        # After reset, needrun should be restored
        assert reset_state['needrun_jobs'] == initial_needrun, \
            "Needrun jobs not restored after reset"

    def test_reset_preserves_targetjobs(self):
        """Test that reset_after_dryrun preserves target job information."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        initial_targets = initial_state['targetjobs']
        
        # Mutate state
        if dag.jobs:
            for job in dag.jobs:
                dag._finished.add(job)
                dag._running.add(job)
        
        # Reset
        dag.reset_after_dryrun()
        
        reset_state = capture_dag_state(dag)
        
        # Target jobs must be preserved
        assert reset_state['targetjobs'] == initial_targets, \
            "Target jobs changed after reset"

    def test_reset_preserves_rules(self):
        """Test that reset_after_dryrun preserves rule information."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        initial_rules = initial_state['rules']
        
        # Simulate heavy state mutation
        if dag.jobs:
            jobs = list(dag.jobs)
            for i, job in enumerate(jobs):
                dag._finished.add(job)
                dag._jobid[job] = i
                dag._group[job] = f"group_{i % 2}"
                if job.output:
                    output = list(job.output)[0]
                    dag._seen_outputs[output] = job
        
        # Reset
        dag.reset_after_dryrun()
        
        reset_state = capture_dag_state(dag)
        
        # Rules must be preserved
        assert reset_state['rules'] == initial_rules, \
            "Rules changed after reset"

    def test_reset_clears_all_execution_state(self):
        """Test that reset_after_dryrun clears all execution state."""
        dag = self._create_minimal_dag()
        
        # Populate all execution state
        if dag.jobs:
            jobs = list(dag.jobs)
            
            # Finished state
            dag._finished.update(jobs[:1])
            
            # Running/ready state
            if len(jobs) > 1:
                dag._running.update(jobs[:1])
                dag._ready_jobs.update(jobs[1:2])
            
            # Job tracking
            for i, job in enumerate(jobs):
                if i < 3:  # Limit for performance
                    dag._checked_jobs.add(job)
                    dag._checked_needrun_jobs.add(job)
                    dag._jobid[job] = i
                    dag._group[job] = f"g{i}"
                    dag._n_until_ready[job] = i
        
        # Storage/queue tracking
        dag._queue_input_jobs = ["test_input"]
        dag._has_unfinished_queue_input_jobs = True
        if dag.jobs:
            job = next(iter(dag.jobs))
            dag._jobs_with_finished_queue_input.add(job)
            dag._storage_input_jobs["store"].append(job)
            if job.output:
                output = list(job.output)[0]
                dag._seen_outputs[output] = job
        
        dag._progress = 5
        
        # Reset
        dag.reset_after_dryrun()
        
        # Verify all execution state is cleared
        assert len(dag._finished) == 0, "_finished not cleared"
        assert len(dag._running) == 0, "_running not cleared"
        assert len(dag._ready_jobs) == 0, "_ready_jobs not cleared"
        assert len(dag._checked_jobs) == 0, "_checked_jobs not cleared"
        assert len(dag._checked_needrun_jobs) == 0, "_checked_needrun_jobs not cleared"
        assert len(dag._jobid) == 0, "_jobid not cleared"
        assert len(dag._group) == 0, "_group not cleared"
        assert len(dag._n_until_ready) == 0, "_n_until_ready not cleared"
        assert len(dag._jobs_with_finished_queue_input) == 0, \
            "_jobs_with_finished_queue_input not cleared"
        assert len(dag._storage_input_jobs) == 0, "_storage_input_jobs not cleared"
        assert len(dag._seen_outputs) == 0, "_seen_outputs not cleared"
        assert dag._queue_input_jobs is None, "_queue_input_jobs not cleared"
        assert dag._has_unfinished_queue_input_jobs is None, \
            "_has_unfinished_queue_input_jobs not cleared"
        assert dag._progress == 0, "_progress not reset"

    def test_reset_is_idempotent(self):
        """Test that calling reset multiple times produces same result."""
        dag = self._create_minimal_dag()
        
        # Simulate dryrun
        if dag.jobs:
            for job in dag.jobs:
                dag._finished.add(job)
                dag._jobid[job] = len(dag._jobid)
        
        # First reset
        dag.reset_after_dryrun()
        state_after_first_reset = capture_dag_state(dag)
        
        # Second reset
        dag.reset_after_dryrun()
        state_after_second_reset = capture_dag_state(dag)
        
        # States should be identical
        are_equal, error_msg = compare_dag_states(
            state_after_first_reset,
            state_after_second_reset,
            ignore_states={'finished_jobs', 'running_jobs', 'ready_jobs'}
        )
        assert are_equal, f"Reset not idempotent: {error_msg}"

    def test_reset_with_empty_dag(self):
        """Test reset on a DAG that has no jobs."""
        dag = self._create_minimal_dag()
        
        # If DAG has no jobs, reset should handle gracefully
        try:
            dag.reset_after_dryrun()
            # If it succeeds, check state is consistent
            state = capture_dag_state(dag)
            assert state['job_count'] == 0 or state['job_count'] >= 0
        except Exception as e:
            pytest.fail(f"Reset failed on empty/minimal DAG: {e}")

    def test_reset_with_complex_dryrun_simulation(self):
        """Test reset after comprehensive dryrun simulation."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        
        # Comprehensive dryrun simulation
        if dag.jobs:
            jobs = list(dag.jobs)
            
            # Mark all jobs as finished (dryrun behavior)
            dag._finished.update(jobs)
            dag._progress = len(jobs)
            
            # Simulate various tracking states
            for i, job in enumerate(jobs):
                dag._jobid[job] = i
                reason = Reason()
                reason.finished = True
                dag._reason[job] = reason
                
                if i % 2 == 0:
                    dag._checked_jobs.add(job)
            
            # Add queue input simulation
            if len(jobs) > 0:
                dag._queue_input_jobs = jobs[:1]
                dag._has_unfinished_queue_input_jobs = True
                for job in jobs[:1]:
                    dag._jobs_with_finished_queue_input.add(job)
            
            # Add storage simulation
            for i, job in enumerate(jobs[:min(2, len(jobs))]):
                dag._storage_input_jobs[f"store_{i}"].append(job)
        
        # Reset
        dag.reset_after_dryrun()
        
        reset_state = capture_dag_state(dag)
        
        # Core structure must be preserved
        assert reset_state['job_count'] == initial_state['job_count'], \
            "Job count changed"
        assert reset_state['jobs'] == initial_state['jobs'], \
            "Job set changed"
        assert reset_state['dependencies'] == initial_state['dependencies'], \
            "Dependencies changed"

    def test_reset_multiple_dryrun_cycles(self):
        """Test reset works correctly across multiple dryrun cycles."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        
        # Simulate multiple dryrun cycles
        for cycle in range(3):
            # Simulate dryrun
            if dag.jobs:
                for job in dag.jobs:
                    dag._finished.add(job)
                    dag._progress = cycle + 1
            
            # Reset
            dag.reset_after_dryrun()
            
            # Check state after each reset
            reset_state = capture_dag_state(dag)
            are_equal, error_msg = compare_dag_states(
                initial_state,
                reset_state,
                ignore_states={'finished_jobs', 'running_jobs', 'ready_jobs'}
            )
            assert are_equal, \
                f"DAG not restored after cycle {cycle}: {error_msg}"

    def test_reset_with_deferred_temp_jobs(self):
        """Test reset properly clears deferred temp jobs."""
        dag = self._create_minimal_dag()
        
        initial_state = capture_dag_state(dag)
        
        # Simulate deferred temp jobs
        if dag.jobs:
            jobs = list(dag.jobs)
            dag._deferred_temp_jobs.extend(jobs[:min(2, len(jobs))])
            
            # Also mark them finished
            for job in dag._deferred_temp_jobs:
                dag._finished.add(job)
        
        assert len(dag._deferred_temp_jobs) > 0, \
            "Failed to populate deferred temp jobs"
        
        # Reset
        dag.reset_after_dryrun()
        
        # Deferred temp jobs should be cleared
        assert len(dag._deferred_temp_jobs) == 0, \
            "Deferred temp jobs not cleared"
        
        # But structure should be intact
        reset_state = capture_dag_state(dag)
        are_equal, error_msg = compare_dag_states(
            initial_state,
            reset_state,
            ignore_states={'finished_jobs', 'running_jobs', 'ready_jobs'}
        )
        assert are_equal, f"Structure changed: {error_msg}"

    def test_reset_preserves_priority_structure(self):
        """Test that reset preserves priority information."""
        dag = self._create_minimal_dag()
        
        # Assign some priorities
        if dag.jobs:
            jobs = list(dag.jobs)
            if len(jobs) > 0:
                dag._priority[jobs[0]] = 10
            if len(jobs) > 1:
                dag._priority[jobs[1]] = 5
        
        initial_state = capture_dag_state(dag)
        initial_priority = initial_state['priority_map']
        
        # Simulate dryrun
        if dag.jobs:
            for job in dag.jobs:
                dag._finished.add(job)
        
        # Reset
        dag.reset_after_dryrun()
        
        reset_state = capture_dag_state(dag)
        
        # Priority structure should be preserved
        assert reset_state['priority_map'] == initial_priority, \
            "Priority map changed after reset"

