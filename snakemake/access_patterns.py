from dataclasses import dataclass
from typing import Optional
from snakemake.io import flag
from snakemake_interface_common.exceptions import WorkflowError

VALID_PATTERNS = ["head", "tail", "random", "sequential"]


@dataclass
class AccessPattern:
    pattern: str
    lines: Optional[int] = None
    bytes: Optional[int] = None


def access_pattern(
    path, pattern: str, lines: Optional[int] = None, bytes: Optional[int] = None
):
    if pattern not in VALID_PATTERNS:
        raise WorkflowError(
            f"Invalid access pattern: {pattern}. Valid patterns are: {', '.join(VALID_PATTERNS)}"
        )
    return flag(path, "access_pattern", AccessPattern(pattern, lines, bytes))
