"""Simple behavior tests for PSB telemetry."""

import subprocess as sp
import sys
import pytest
from pathlib import Path

from .common import run, dpath


class TestTelemetryCLI:
    """Test CLI parameter handling for --share-benchmark."""

    def test_share_benchmark_requires_collector(self, tmp_path):
        """Test that --share-benchmark without --share-benchmark-collector fails."""
        # Create minimal Snakefile
        snakefile = tmp_path / "Snakefile"
        snakefile.write_text("""
rule all:
    output: "test.txt"
    shell: "echo hi > {output}"
""")

        result = sp.run(
            [sys.executable, "-m", "snakemake", "--share-benchmark", "-c1"],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )

        # Should fail with error message
        assert result.returncode != 0
        assert (
            "share-benchmark-collector" in result.stderr.lower()
            or "share-benchmark-collector" in result.stdout.lower()
        )

    def test_share_benchmark_with_collector_accepted(self, tmp_path):
        """Test that both flags together are accepted."""
        snakefile = tmp_path / "Snakefile"
        snakefile.write_text("""
rule all:
    output: "test.txt"
    shell: "echo hi > {output}"
""")

        # Use fake collector URL - workflow should start (may fail on flush, but that's ok)
        result = sp.run(
            [
                sys.executable,
                "-m",
                "snakemake",
                "--share-benchmark",
                "--share-benchmark-collector",
                "http://localhost:9999",
                "-c1",
                "-n",
            ],  # -n for dry-run
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )

        # Should not fail due to missing collector URL
        assert "--share-benchmark-collector" not in result.stderr


class TestAnnotationResolution:
    """Test PSB annotation resolution: per-rule > config > auto-detect."""

    def test_auto_detect_extracts_tool_from_shell(self):
        """Test that tool name is auto-detected from shell command."""
        from snakemake.telemetry import parse_shell_tool

        tool, params = parse_shell_tool("samtools sort -@ 4 input.bam")
        assert tool == "samtools"
        assert "sort" in params

    def test_auto_detect_skips_generic_interpreters(self):
        """Test that python/bash are skipped, revealing actual tool."""
        from snakemake.telemetry import parse_shell_tool

        tool, params = parse_shell_tool("python scripts/process.py --input data")
        assert tool == "scripts/process.py"
        assert tool != "python"


class TestClientInit:
    """Test basic PSB client initialization."""

    def test_client_can_be_initialized(self):
        """Test that PSB client initializes without errors."""
        from snakemake.telemetry import BenchmarkTelemetryClient

        client = BenchmarkTelemetryClient(
            endpoint="http://test:8080", token="test-token"
        )

        assert client.endpoint == "http://test:8080"
        assert client.token == "test-token"
        assert client.session_id is not None

    def test_global_init_creates_client(self):
        """Test that init_psb creates a global client."""
        from snakemake.telemetry import init_psb, get_psb_client

        init_psb(endpoint="http://test:8080", sm_version="8.0.0")

        client = get_psb_client()
        assert client is not None
        assert client.endpoint == "http://test:8080"

    def test_client_accumulates_records(self):
        """Test that client can accumulate multiple records."""
        from snakemake.telemetry import BenchmarkTelemetryClient

        client = BenchmarkTelemetryClient("http://test", "token")
        client.set_environment(sm_version="8.0.0", deploy_mode="host")

        client.add_record(
            record_id="r1",
            tool="tool1",
            runtime_sec=10.0,
            max_rss_mb=100.0,
            cpu_percent=50.0,
        )

        client.add_record(
            record_id="r2",
            tool="tool2",
            runtime_sec=5.0,
            max_rss_mb=50.0,
            cpu_percent=30.0,
        )

        # Buffer should contain both records
        assert len(client._buffer) == 2


class TestPSBCollectorIntegration:
    """Integration test with a real (fake) HTTP server."""

    def test_client_sends_valid_jsonl_to_collector(self, tmp_path):
        """Behavior: client should send valid JSONL conforming to PSB spec."""
        import threading
        import http.server
        import socketserver
        from urllib.parse import urlparse

        from snakemake.telemetry import init_psb, add_psb_record, flush_psb

        import json

        # Storage for received data
        received_records = []

        class FakePSBHandler(http.server.BaseHTTPRequestHandler):
            def do_POST(self):
                if self.path == "/v1/telemetry":
                    content_length = int(self.headers["Content-Length"])
                    body = self.rfile.read(content_length).decode("utf-8")

                    # Verify it's valid JSONL
                    lines = body.strip().split("\n")
                    for line in lines:
                        record = json.loads(line)  # Will raise if invalid JSON
                        received_records.append(record)

                    # Send success response
                    response = json.dumps(
                        {"accepted": len(lines), "duplicates": 0, "rejected": 0}
                    ).encode("utf-8")

                    self.send_response(201)
                    self.send_header("Content-Type", "application/json")
                    self.send_header("Content-Length", str(len(response)))
                    self.end_headers()
                    self.wfile.write(response)
                else:
                    self.send_response(404)
                    self.end_headers()

            def log_message(self, format, *args):
                pass  # Suppress logs

        # Start fake server on random port
        with socketserver.TCPServer(("127.0.0.1", 0), FakePSBHandler) as httpd:
            port = httpd.server_address[1]
            server_thread = threading.Thread(target=httpd.serve_forever)
            server_thread.daemon = True
            server_thread.start()

            try:
                # Initialize PSB client pointing to fake server
                init_psb(
                    endpoint=f"http://127.0.0.1:{port}",
                    sm_version="9.0.0",
                    deploy_mode="host",
                    workflow_url="https://github.com/test/repo",
                    workflow_version="v1.2.3",
                )

                # Add some records
                add_psb_record(
                    record_id="test-rec-1",
                    tool="samtools",
                    command="samtools",
                    params="sort -@ 4",
                    runtime_sec=10.5,
                    max_rss_mb=256.0,
                    cpu_percent=85.0,
                    threads=4,
                    inputs=[{"type": ".bam", "size": 1024000}],
                    outputs=[{"type": ".sorted.bam", "size": 1024000}],
                )

                add_psb_record(
                    record_id="test-rec-2",
                    tool="bwa-mem2",
                    command="bwa-mem2",
                    params="mem -t 8",
                    runtime_sec=45.2,
                    max_rss_mb=2048.0,
                    cpu_percent=750.0,
                    threads=8,
                )

                # Flush and verify
                flush_psb()

                # Give server a moment to process
                import time

                time.sleep(0.1)

                # Verify received records
                assert len(received_records) == 2

                # Verify required fields in spec
                for record in received_records:
                    # Session metadata
                    assert "session_id" in record
                    assert "sm_version" in record
                    assert record["sm_version"] == "9.0.0"
                    assert "deploy_mode" in record
                    assert "workflow_url" in record
                    assert record["workflow_url"] == "https://github.com/test/repo"
                    assert "workflow_version" in record
                    assert record["workflow_version"] == "v1.2.3"

                    # Environment metadata
                    assert "host_hash" in record
                    assert "cpu_model" in record
                    assert "cpu_cores" in record
                    assert "os" in record

                    # Record-level required fields
                    assert "record_id" in record
                    assert "tool" in record
                    assert "runtime_sec" in record
                    assert "max_rss_mb" in record
                    assert "cpu_percent" in record

                    # Verify types
                    assert isinstance(record["runtime_sec"], (int, float))
                    assert isinstance(record["max_rss_mb"], (int, float))
                    assert isinstance(record["cpu_percent"], (int, float))

                # Verify specific tool names made it through
                tools = [r["tool"] for r in received_records]
                assert "samtools" in tools
                assert "bwa-mem2" in tools

                # Verify inputs/outputs format
                rec1 = [r for r in received_records if r["tool"] == "samtools"][0]
                assert "inputs" in rec1
                assert isinstance(rec1["inputs"], list)
                assert len(rec1["inputs"]) == 1
                assert rec1["inputs"][0]["type"] == ".bam"
                assert rec1["inputs"][0]["size"] == 1024000

            finally:
                httpd.shutdown()
                # Cleanup
                import snakemake.telemetry._client as client_module

                client_module._client = None
