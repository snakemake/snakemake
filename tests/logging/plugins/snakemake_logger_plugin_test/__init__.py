"""Test logger plugin.

* Outputs records in JSONL format.
* Writes to a file or stdout based on CLI options.
* Other LogHandlerBase properties set by CLI options (but formatter not used if set).
* First line is special "logger_info" event with information on whether default handler/filter have
  been added.
"""

from typing import Any, Optional
from dataclasses import dataclass, field
import logging
import json
import sys

from snakemake_interface_logger_plugins.settings import LogHandlerSettingsBase
from snakemake_interface_logger_plugins.base import LogHandlerBase


DEFAULT_FORMATTER = logging.Formatter()


@dataclass
class LogHandlerSettings(LogHandlerSettingsBase):
    outfile: Optional[str] = field(
        default=None,
        metadata={
            "help": "Path to output file",
        },
    )
    has_formatter: bool = field(
        default=False,
        metadata={
            "help": "Value of has_formatter()",
        },
    )
    has_filter: bool = field(
        default=False,
        metadata={
            "help": "Value of has_filter()",
        },
    )
    needs_rulegraph: bool = field(
        default=False,
        metadata={
            "help": "Value of needs_rulegraph()",
        },
    )


class LogHandler(LogHandlerBase):
    settings: LogHandlerSettings
    handler = logging.Handler

    def __post_init__(self) -> None:

        # These should be set by the logging manager when the handler is configured, which happens
        # after construction.
        assert self.formatter is None
        assert not self.filters

        if self.settings.outfile is None:
            # Output to stderr
            self.stream = sys.stderr
            self.needs_close = False
        else:
            self.stream = open(self.settings.outfile, "w")
            self.needs_close = True

        self.first_record = True

    def emit(self, record: logging.LogRecord) -> None:
        # Emit info about logger first
        if self.first_record:
            self._emit(
                dict(
                    event="logger_info",
                    formatter_set=self.formatter is not None,
                    filter_added=bool(self.filters),
                )
            )
            self.first_record = False

        event = getattr(record, "event", None)
        self._emit(
            dict(
                event=None if event is None else str(event),
                msg=DEFAULT_FORMATTER.format(record),
                level=record.levelname,
            )
        )

    def _emit(self, data: dict[str, Any]) -> None:
        json.dump(data, self.stream)
        self.stream.write("\n")

    def close(self):
        super().close()
        # Close file handle
        if self.needs_close:
            self.stream.close()

    @property
    def writes_to_stream(self) -> bool:
        return self.settings.outfile is None

    @property
    def writes_to_file(self) -> bool:
        return self.settings.outfile is not None

    @property
    def has_filter(self) -> bool:
        return self.settings.has_filter

    @property
    def has_formatter(self) -> bool:
        return self.settings.has_formatter

    @property
    def needs_rulegraph(self) -> bool:
        return self.settings.needs_rulegraph

    @property
    def baseFilename(self) -> str | None:
        return self.settings.outfile
