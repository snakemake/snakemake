"""Telemetry client for Snakemake.

Thread-safe singleton that accumulates benchmark records during a workflow run
and flushes them as a single JSONL POST to the PSB collector on shutdown.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import threading
import urllib.error
import urllib.request
import uuid
from typing import Any

from snakemake.logging import logger

DEFAULT_TOKEN = "dev-secret"


class BenchmarkTelemetryClient:
    """Accumulates benchmark records and flushes as a single JSONL POST."""

    def __init__(self, endpoint: str, token: str, session_id: str | None = None):
        self.endpoint = endpoint.rstrip("/")
        self.token = token
        self.session_id = session_id or str(uuid.uuid4())
        self.session_start = datetime.datetime.now(datetime.timezone.utc).isoformat()
        self._cookie = str(uuid.uuid4())
        self._env: dict[str, str] = {}
        self._session_meta: dict[str, str] = {}
        self._buffer: list[dict[str, Any]] = []

    def make_record_id(self, rule_name: str, wildcards: dict | None = None) -> str:
        """Deterministic record ID from session start + rule name + wildcards."""
        key = f"{self.session_start}:{rule_name}"
        if wildcards:
            wc_str = ",".join(f"{k}={v}" for k, v in sorted(wildcards.items()))
            key = f"{key}:{wc_str}"
        h = hashlib.sha256(key.encode()).hexdigest()
        return h[:16]

    def set_session_meta(
        self, workflow_url: str = "", workflow_version: str = ""
    ) -> None:
        """Set session-level metadata (workflow URL and version)."""
        self._session_meta = {}
        if workflow_url:
            self._session_meta["workflow_url"] = workflow_url
        if workflow_version:
            self._session_meta["workflow_version"] = workflow_version

    def set_environment(self, **kwargs: str) -> None:
        self._env = dict(kwargs)

    def add_record(
        self,
        *,
        record_id: str,
        tool: str,
        runtime_sec: float,
        max_rss_mb: float,
        cpu_percent: float,
        **optional: Any,
    ) -> None:
        if not tool or runtime_sec <= 0 or not record_id:
            return  # silently skip invalid
        record: dict[str, Any] = {
            "record_id": record_id,
            "tool": tool,
            "runtime_sec": runtime_sec,
            "max_rss_mb": max_rss_mb,
            "cpu_percent": cpu_percent,
        }
        for key in (
            "command",
            "params",
            "inputs",
            "outputs",
            "threads",
            "exit_code",
            "load_avg",
            "mem_avail_mb",
            "swap_used_mb",
            "io_wait_pct",
            "shell_block",
            "max_vms_mb",
            "max_uss_mb",
            "max_pss_mb",
            "io_in_mb",
            "io_out_mb",
            "cpu_time_sec",
            "resources",
            "tool_version",
            "category",
        ):
            if key in optional:
                record[key] = optional[key]
        self._buffer.append(record)

    def flush(self, timeout: float = 10.0) -> dict | None:
        if not self._buffer or not self._env:
            return None

        lines = []
        for record in self._buffer:
            line = {
                "session_id": self.session_id,
                **self._session_meta,
                **self._env,
                **record,
            }
            lines.append(json.dumps(line, separators=(",", ":")))

        body = "\n".join(lines).encode("utf-8")
        url = f"{self.endpoint}/v1/telemetry"
        req = urllib.request.Request(url, data=body, method="POST")
        req.add_header("Content-Type", "text/jsonl")
        req.add_header("User-Agent", "snakemake-psb/1.0")
        req.add_header("X-PSB-Token", self.token)
        req.add_header("X-PSB-Nonce", str(uuid.uuid4()))
        req.add_header("Cookie", f"_psb_session={self._cookie}")

        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                resp_body = resp.read().decode("utf-8")
                if resp.status == 201:
                    result = json.loads(resp_body)
                    self._buffer.clear()
                    return result
        except Exception:
            pass  # buffer preserved for retry
        return None


# Global client singleton
_client: BenchmarkTelemetryClient | None = None
_lock = threading.Lock()


def init_psb(
    endpoint: str,
    sm_version: str,
    token: str = DEFAULT_TOKEN,
    deploy_mode: str = "host",
    workflow_url: str = "",
    workflow_version: str = "",
) -> None:
    """Initialize the global PSB client. Call once at workflow start."""
    from ._environment import collect_environment

    global _client
    with _lock:
        _client = BenchmarkTelemetryClient(endpoint, token)
        _client.set_environment(**collect_environment(sm_version, deploy_mode))
        _client.set_session_meta(
            workflow_url=workflow_url, workflow_version=workflow_version
        )
    logger.info(f"Telemetry enabled (session {_client.session_id})")


def get_psb_client() -> BenchmarkTelemetryClient | None:
    """Return the global PSB client, or None if not initialized."""
    with _lock:
        return _client


def add_psb_record(**kwargs: Any) -> None:
    """Thread-safe append of a benchmark record. No-op if not initialized."""
    with _lock:
        if _client is not None:
            try:
                _client.add_record(**kwargs)
            except Exception as e:
                logger.debug(f"Telemetry: failed to add record: {e}")


def flush_psb() -> None:
    """Flush all buffered records. Fire-and-forget, never raises."""
    with _lock:
        client = _client
    if client is None:
        return
    try:
        result = client.flush(timeout=10.0)
        if result:
            logger.info(
                f"Telemetry: {result.get('accepted', 0)} accepted, "
                f"{result.get('duplicates', 0)} duplicates, "
                f"{result.get('rejected', 0)} rejected"
            )
        else:
            n = len(client._buffer)
            if n:
                logger.warning(f"Telemetry: flush failed, {n} records lost")
    except Exception as e:
        logger.warning(f"Telemetry flush error: {e}")


def get_session_url() -> str | None:
    """Return the full URL to view the current session, or None."""
    with _lock:
        if _client is not None:
            return f"{_client.endpoint}/session/{_client.session_id}"
    return None


def submit_benchmark_records(
    bench_records,
    job_rule,
    params,
    wildcards,
    threads,
    input_files,
    output_files,
    psb_snapshot,
):
    """Submit benchmark records to Telemetry.

    Handles annotation resolution, file metadata collection, and record submission.

    Args:
        bench_records: List of BenchmarkRecord objects
        job_rule: The rule being executed
        params: Rule parameters (may contain _psb_* annotations)
        wildcards: Wildcards dict
        threads: Number of threads
        input_files: List of input file paths
        output_files: List of output file paths
        psb_snapshot: System snapshot captured before execution
    """
    import json
    import os

    from ._parsing import GENERIC_INTERPRETERS, parse_shell_tool, get_file_type
    from ._metrics import read_iowait_ticks, compute_iowait_pct

    # Resolve _psb_* annotations
    # Precedence: rule-level param > workflow config > auto-detect
    psb_cfg = {}
    try:
        psb_cfg = job_rule.workflow.config.get("psb", {}) or {}
    except Exception:
        pass

    def _annot(name):
        """Get annotation: rule param > config default."""
        val = getattr(params, f"_psb_{name}", None)
        if val is not None:
            return val
        return psb_cfg.get(name)

    # Auto-detect tool/params from shell template
    auto_tool, auto_params = parse_shell_tool(job_rule.shellcmd)

    # _psb_tool replaces tool field
    psb_tool = _annot("tool")
    tool = psb_tool if psb_tool else auto_tool

    # _psb_primary_cmd replaces command field
    psb_cmd = _annot("primary_cmd")
    command = psb_cmd if psb_cmd else tool

    # _psb_params replaces params field
    psb_params = _annot("params")
    if psb_params is not None:
        params_str = psb_params
    elif psb_cmd and not psb_params:
        # When _psb_primary_cmd is set without _psb_params, clear params
        params_str = ""
    else:
        params_str = auto_params

    # tool_version and category from annotations
    tool_version = _annot("tool_version") or ""
    category = _annot("category") or ""

    # GENERIC_INTERPRETERS drop rule always applies
    psb_client = get_psb_client()
    if not tool or tool in GENERIC_INTERPRETERS or not psb_client:
        return

    # Compute iowait delta
    io_wait_pct = 0.0
    if psb_snapshot:
        iowait_end = read_iowait_ticks()
        io_wait_pct = compute_iowait_pct(psb_snapshot["iowait_start"], iowait_end)

    for br in bench_records:
        if br.running_time is None or br.running_time <= 0:
            continue
        wc_dict = dict(wildcards) if wildcards else None
        record_id = psb_client.make_record_id(job_rule.name, wc_dict)

        # Build per-file input entries
        inputs_list: list[dict[str, object]] = []
        # Ensure input_files is iterable as a list, not a string
        # Namedlist works correctly, but plain strings would iterate by character
        if isinstance(input_files, str):
            input_files = [input_files] if input_files.strip() else []
        if input_files:
            for inp in input_files:
                try:
                    size = os.path.getsize(inp)
                except OSError:
                    size = 0
                t = get_file_type(inp) or ""
                inputs_list.append({"type": t, "size": size})

        # Build per-file output entries
        outputs_list: list[dict[str, object]] = []
        # Ensure output_files is iterable as a list, not a string
        if isinstance(output_files, str):
            output_files = [output_files] if output_files.strip() else []
        if output_files:
            for out in output_files:
                try:
                    size = os.path.getsize(out)
                except OSError:
                    size = 0
                t = get_file_type(out) or ""
                outputs_list.append({"type": t, "size": size})

        # Resources as JSON string
        res_str = ""
        try:
            res_str = json.dumps(br.parse_resources())
        except Exception:
            pass

        add_psb_record(
            record_id=record_id,
            tool=tool,
            command=command,
            params=params_str,
            shell_block=job_rule.shellcmd,
            runtime_sec=br.running_time,
            threads=threads or 0,
            max_rss_mb=br.max_rss or 0.0,
            cpu_percent=br.cpu_usage or 0.0,
            max_vms_mb=br.max_vms or 0.0,
            max_uss_mb=br.max_uss or 0.0,
            max_pss_mb=br.max_pss or 0.0,
            io_in_mb=br.io_in or 0.0,
            io_out_mb=br.io_out or 0.0,
            cpu_time_sec=br.cpu_time or 0.0,
            resources=res_str,
            tool_version=tool_version,
            category=category,
            inputs=inputs_list,
            outputs=outputs_list,
            exit_code=0,
            load_avg=psb_snapshot["load_avg"] if psb_snapshot else 0.0,
            mem_avail_mb=psb_snapshot["mem_avail_mb"] if psb_snapshot else 0,
            swap_used_mb=psb_snapshot["swap_used_mb"] if psb_snapshot else 0,
            io_wait_pct=io_wait_pct,
        )
