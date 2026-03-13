"""Parsing utilities for telemetry.

Extracts tool names and parameters from shell commands, and determines
file types from file paths.
"""

from __future__ import annotations

import os

# Generic interpreters that should be skipped when detecting tool names
GENERIC_INTERPRETERS = frozenset(
    {"python", "python3", "bash", "sh", "Rscript", "perl", "java", "ruby"}
)


def parse_shell_tool(shellcmd: str) -> tuple[str | None, str]:
    """Extract the tool and params from a shell command template.

    Skips comment lines, blank lines, and leading generic interpreter tokens.
    Returns (tool, params) where tool is the first non-interpreter token
    and params is everything after it.

    Args:
        shellcmd: Shell command string, potentially multi-line

    Returns:
        Tuple of (tool_name, params_string). Returns (None, "") if no tool found.

    Examples:
        "samtools sort -@ 4 {input}"       -> ("samtools", "sort -@ 4 {input}")
        "python scripts/foo.py --arg"      -> ("scripts/foo.py", "--arg")
        "# comment\\ngzip -9 {input}"       -> ("gzip", "-9 {input}")
        "python3"                           -> (None, "")
        ""                                  -> (None, "")
    """
    for line in shellcmd.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.split()
        # Skip leading interpreter tokens
        idx = 0
        while idx < len(tokens) and tokens[idx] in GENERIC_INTERPRETERS:
            idx += 1
        if idx >= len(tokens):
            return None, ""
        tool = tokens[idx]
        params = " ".join(tokens[idx + 1 :]) if idx + 1 < len(tokens) else ""
        return tool, params
    return None, ""


def get_file_type(filepath: str) -> str | None:
    """Extract file extension, handling compound extensions like .fastq.gz.

    Args:
        filepath: Path to file (can be relative or absolute). Must be a single
            file path. If a space-separated list of paths is passed, only the
            first path will be processed.

    Returns:
        File extension including the dot, or None if no extension found.
        Compound extensions (e.g., .fastq.gz) are preserved.

    Examples:
        "data.txt" -> ".txt"
        "reads.fastq.gz" -> ".fastq.gz"
        "archive.tar.bz2" -> ".tar.bz2"
        "Makefile" -> None
        "" -> None
    """
    if not filepath:
        return None

    # Handle the case where a space-separated list of files is passed
    # (e.g., from accidentally converting a list to string).
    # We take the first file only as a best-effort approach.
    stripped = filepath.strip()
    if " " in stripped and "." in stripped.split()[0]:
        # Looks like space-separated files (e.g., "file1.txt file2.txt")
        # Take the first file only
        stripped = stripped.split()[0]
        filepath = stripped

    name = os.path.basename(filepath)
    parts = name.split(".")

    # Handle compound extensions (e.g., .fastq.gz, .tar.bz2, .tar.gz)
    # Recognized compression suffixes
    COMPRESSION_SUFFIXES = frozenset({"gz", "bz2", "xz", "zst", "zip", "rar", "7z"})

    if len(parts) >= 3 and parts[-1] in COMPRESSION_SUFFIXES:
        # It's a compound extension like file.tar.gz or reads.fastq.gz
        ext = "." + ".".join(parts[-2:])
    elif len(parts) >= 2:
        # Simple extension like .txt, .fa, .bam
        ext = "." + parts[-1]
    else:
        # No extension found (e.g., "Makefile", "README")
        return None

    return ext
