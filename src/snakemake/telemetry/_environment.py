"""System environment introspection for telemetry.

Collects hardware and OS information: CPU model, features, cores, caches,
frequency, OS version, etc.
"""

from __future__ import annotations

import hashlib
import os
import platform
import socket
import subprocess
import sys


def host_hash() -> str:
    """Return a deterministic hash of the hostname."""
    return hashlib.sha256(socket.gethostname().encode()).hexdigest()


def cpu_model() -> str:
    """Return CPU model name."""
    if sys.platform == "linux":
        try:
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
        except OSError:
            pass
    elif sys.platform == "darwin":
        try:
            out = subprocess.check_output(
                ["sysctl", "-n", "machdep.cpu.brand_string"], text=True, timeout=2
            )
            return out.strip()
        except (subprocess.SubprocessError, FileNotFoundError):
            pass
    return platform.processor() or "unknown"


# Curated CPU feature bitmask. bit positions MUST match the Go server
# (internal/cpufeatures/cpufeatures.go).
_CPU_FEATURE_REGISTRY: list[tuple[str, int]] = [
    # x86 SIMD
    ("sse2", 1 << 0),
    ("pni", 1 << 1),  # Linux name for SSE3
    ("sse3", 1 << 1),
    ("ssse3", 1 << 2),
    ("sse4_1", 1 << 3),
    ("sse4_2", 1 << 4),
    ("avx", 1 << 5),
    ("avx2", 1 << 6),
    ("fma", 1 << 7),
    ("avx512f", 1 << 8),
    ("avx512bw", 1 << 9),
    ("avx512vl", 1 << 10),
    ("avx512dq", 1 << 11),
    ("avx512_vnni", 1 << 12),
    ("avx512vnni", 1 << 12),
    ("f16c", 1 << 13),
    # Bit manipulation
    ("popcnt", 1 << 14),
    ("bmi1", 1 << 15),
    ("bmi2", 1 << 16),
    ("abm", 1 << 17),
    ("lzcnt", 1 << 17),
    # Crypto / hashing
    ("aes", 1 << 18),
    ("aesni", 1 << 18),
    ("sha_ni", 1 << 19),
    ("sha1", 1 << 19),
    ("sha2", 1 << 19),
    ("pclmulqdq", 1 << 20),
    ("pclmul", 1 << 20),
    # Misc x86
    ("rdrand", 1 << 21),
    ("rtm", 1 << 22),
    ("hle", 1 << 22),
    # ARM
    ("neon", 1 << 23),
    ("asimd", 1 << 23),
    ("sve", 1 << 24),
    ("sve2", 1 << 25),
    ("crc32", 1 << 26),
    # AMD
    ("xop", 1 << 27),
]


def raw_cpu_flags() -> list[str]:
    """Return raw CPU flag names from the OS."""
    if sys.platform == "linux":
        try:
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("flags"):
                        return line.split(":", 1)[1].strip().split()
        except OSError:
            pass
    elif sys.platform == "darwin":
        try:
            out = subprocess.check_output(
                ["sysctl", "-n", "machdep.cpu.features"], text=True, timeout=2
            )
            return [f.lower() for f in out.strip().split()]
        except (subprocess.SubprocessError, FileNotFoundError):
            pass
    return []


def encode_cpu_features(flags: list[str]) -> int:
    """Encode a list of CPU flag names into a bitmask."""
    lookup: dict[str, int] = {}
    for name, bit in _CPU_FEATURE_REGISTRY:
        lookup[name] = bit
    mask = 0
    for f in flags:
        f = f.lower().strip()
        if f in lookup:
            mask |= lookup[f]
    return mask


def cpu_cores() -> int:
    """Return the number of physical CPU cores (0 if unknown)."""
    try:
        if sys.platform == "linux":
            # Count unique physical core IDs
            cores: set[tuple[str, str]] = set()
            phys_id = core_id = ""
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("physical id"):
                        phys_id = line.split(":", 1)[1].strip()
                    elif line.startswith("core id"):
                        core_id = line.split(":", 1)[1].strip()
                        cores.add((phys_id, core_id))
            return len(cores) if cores else os.cpu_count() or 0
        elif sys.platform == "darwin":
            out = subprocess.check_output(
                ["sysctl", "-n", "hw.physicalcpu"], text=True, timeout=2
            )
            return int(out.strip())
    except (OSError, subprocess.SubprocessError, ValueError):
        pass
    return os.cpu_count() or 0


def cache_sizes() -> tuple[int, int]:
    """Return (l2_cache_kb, l3_cache_kb). Returns 0 for unknown values."""
    l2 = l3 = 0
    try:
        if sys.platform == "linux":
            import glob as _glob

            for idx_dir in sorted(
                _glob.glob("/sys/devices/system/cpu/cpu0/cache/index*")
            ):
                try:
                    with open(os.path.join(idx_dir, "level")) as f:
                        level = int(f.read().strip())
                    with open(os.path.join(idx_dir, "size")) as f:
                        size_str = f.read().strip()  # e.g. "256K" or "8192K"
                        size_kb = int(size_str.rstrip("K"))
                    if level == 2:
                        l2 = size_kb
                    elif level == 3:
                        l3 = size_kb
                except (OSError, ValueError):
                    continue
        elif sys.platform == "darwin":
            for key, setter in [
                ("hw.l2cachesize", "l2"),
                ("hw.l3cachesize", "l3"),
            ]:
                try:
                    out = subprocess.check_output(
                        ["sysctl", "-n", key], text=True, timeout=2
                    )
                    val = int(out.strip()) // 1024  # bytes -> KB
                    if setter == "l2":
                        l2 = val
                    else:
                        l3 = val
                except (subprocess.SubprocessError, ValueError, FileNotFoundError):
                    pass
    except Exception:
        pass
    return l2, l3


def cpu_freq_mhz() -> int:
    """Return max CPU frequency in MHz (0 if unknown)."""
    try:
        if sys.platform == "linux":
            with open("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq") as f:
                return int(f.read().strip()) // 1000  # kHz -> MHz
        elif sys.platform == "darwin":
            out = subprocess.check_output(
                ["sysctl", "-n", "hw.cpufrequency_max"], text=True, timeout=2
            )
            return int(out.strip()) // 1_000_000  # Hz -> MHz
    except (OSError, subprocess.SubprocessError, ValueError, FileNotFoundError):
        pass
    return 0


def platform_os() -> str:
    """Return normalized OS name."""
    mapping = {
        "linux": "linux",
        "darwin": "darwin",
        "freebsd": "freebsd",
        "win32": "windows",
        "cygwin": "windows",
    }
    return mapping.get(sys.platform, sys.platform)


def collect_environment(sm_version: str, deploy_mode: str = "host") -> dict:
    """Collect all environment metadata for telemetry session.

    Args:
        sm_version: Snakemake version string
        deploy_mode: Deployment method (e.g., "host", "conda", "conda+apptainer")

    Returns:
        Dictionary with CPU, OS, and version metadata
    """
    raw_flags = raw_cpu_flags()
    l2, l3 = cache_sizes()
    return {
        "host_hash": host_hash(),
        "cpu_model": cpu_model(),
        "cpu_features": encode_cpu_features(raw_flags),
        "cpu_cores": cpu_cores(),
        "l2_cache_kb": l2,
        "l3_cache_kb": l3,
        "cpu_freq_mhz": cpu_freq_mhz(),
        "os": platform_os(),
        "kernel_version": platform.release(),
        "kernel_string": platform.platform(),
        "sm_version": sm_version,
        "deploy_mode": deploy_mode,
    }
