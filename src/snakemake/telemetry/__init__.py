"""PSB (Planetary Scale Benchmark) telemetry for Snakemake.

Public API for telemetry collection, initialization, and data submission.

Protocol specification:
https://github.com/btraven00/psb/blob/main/docs/spec.md
"""

from ._client import (
    DEFAULT_TOKEN,
    BenchmarkTelemetryClient,
    add_psb_record,
    flush_psb,
    get_psb_client,
    get_session_url,
    init_psb,
    submit_benchmark_records,
)
from ._environment import collect_environment
from ._metrics import (
    capture_system_snapshot,
    compute_iowait_pct,
    read_iowait_ticks,
    read_loadavg,
    read_mem_swap,
)
from ._parsing import GENERIC_INTERPRETERS, get_file_type, parse_shell_tool

__all__ = [
    # Client
    "BenchmarkTelemetryClient",
    "init_psb",
    "get_psb_client",
    "add_psb_record",
    "flush_psb",
    "get_session_url",
    "submit_benchmark_records",
    "DEFAULT_TOKEN",
    # Environment
    "collect_environment",
    # Metrics
    "capture_system_snapshot",
    "read_loadavg",
    "read_mem_swap",
    "read_iowait_ticks",
    "compute_iowait_pct",
    # Parsing
    "parse_shell_tool",
    "get_file_type",
    "GENERIC_INTERPRETERS",
]
