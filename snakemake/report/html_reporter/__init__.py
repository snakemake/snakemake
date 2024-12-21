from dataclasses import dataclass, field
import datetime
import json
import os
from pathlib import Path
import textwrap
import time
from typing import Optional
from zipfile import ZIP_DEFLATED, ZipFile

from jinja2 import Environment, PackageLoader
from snakemake.report.common import data_uri_from_file
from snakemake.report.html_reporter import data

from snakemake.report.rulegraph_spec import rulegraph_spec
from snakemake.report.html_reporter.common import get_resource_as_string, get_result_uri

from snakemake_interface_common.exceptions import WorkflowError
from snakemake_interface_report_plugins.reporter import ReporterBase
from snakemake_interface_report_plugins.settings import ReportSettingsBase


@dataclass
class ReportSettings(ReportSettingsBase):
    path: Optional[Path] = field(
        default=None,
        metadata={
            "help": "Path to the report file (either .html or .zip). Use zip if your "
            "report contains large results or directories with htmlindex as results.",
            "env_var": False,
            "required": True,
        },
    )
    stylesheet_path: Optional[Path] = field(
        default=None,
        metadata={
            "help": "Path to a custom stylesheet for the report.",
            "env_var": False,
            "required": False,
        },
    )


class Reporter(ReporterBase):
    def __post_init__(self):
        self.mode_embedded = True
        if self.settings.path.suffix == ".zip":
            self.mode_embedded = False
        elif self.settings.path.suffix != ".html":
            raise WorkflowError("Report file does not end with .html or .zip")

        self.custom_stylesheet = None
        if self.settings.stylesheet_path is not None:
            try:
                with open(self.settings.stylesheet_path) as s:
                    self.custom_stylesheet = s.read()
            except Exception as e:
                raise WorkflowError("Unable to read custom report stylesheet.", e)

        self.env = Environment(
            loader=PackageLoader("snakemake", "report/html_reporter/template"),
            trim_blocks=True,
            lstrip_blocks=True,
        )
        self.env.filters["get_resource_as_string"] = get_resource_as_string

    def render(self):
        # prepare runtimes
        runtimes = [
            {"rule": rec.rule, "runtime": rec.endtime - rec.starttime}
            for rec in self.jobs
        ]

        def get_datetime(rectime):
            try:
                return datetime.datetime.fromtimestamp(rectime).isoformat()
            except OSError:
                return None

        # prepare end times
        timeline = [
            {
                "rule": rec.rule,
                "starttime": get_datetime(rec.starttime),
                "endtime": get_datetime(rec.endtime),
            }
            for rec in self.jobs
        ]

        # rulegraph
        rulegraph, _, _ = rulegraph_spec(self.dag)

        # record time
        now = f"{datetime.datetime.now().ctime()} {time.tzname[0]}"

        try:
            from pygments.formatters import HtmlFormatter
        except ImportError:
            raise WorkflowError(
                "Python package pygments must be installed to create reports."
            )

        categories = data.render_categories(self.results)
        rendered_results = data.render_results(self.results, self.mode_embedded)
        rulegraph = data.render_rulegraph(
            rulegraph["nodes"], rulegraph["links"], rulegraph["links_direct"]
        )
        rules = data.render_rules(self.rules)
        runtimes = data.render_runtimes(runtimes)
        timeline = data.render_timeline(timeline)
        packages = data.get_packages()

        template = self.env.get_template("index.html.jinja2")

        rendered = template.render(
            results=rendered_results,
            categories=categories,
            rulegraph=rulegraph,
            rules=rules,
            workflow_desc=json.dumps(self.workflow_description),
            runtimes=runtimes,
            timeline=timeline,
            packages=packages,
            pygments_css=HtmlFormatter(style="stata-dark").get_style_defs(".source"),
            custom_stylesheet=self.custom_stylesheet,
            logo=data_uri_from_file(Path(__file__).parent / "template" / "logo.svg"),
            now=now,
        )

        # TODO look into supporting .WARC format, also see (https://webrecorder.io)

        if not self.mode_embedded:
            with ZipFile(
                self.settings.path, compression=ZIP_DEFLATED, mode="w"
            ) as zipout:
                folder = Path(Path(self.settings.path).stem)
                # store results in data folder
                for subcats in self.results.values():
                    for catresults in subcats.values():
                        for result in catresults:
                            if self.mode_embedded and result.aux_files:
                                raise WorkflowError(
                                    "Directory marked for inclusion in report. "
                                    "This is unsupported when requesting a pure HTML report. "
                                    "Please use store as zip instead (--report report.zip)."
                                )
                            # write raw data
                            zipout.write(
                                result.path,
                                str(
                                    folder.joinpath(
                                        get_result_uri(result, self.mode_embedded)
                                    )
                                ),
                            )
                            # write aux files
                            parent = folder.joinpath(
                                get_result_uri(result, self.mode_embedded)
                            ).parent
                            for aux_path in result.aux_files:
                                zipout.write(
                                    aux_path,
                                    str(
                                        parent.joinpath(
                                            os.path.relpath(
                                                aux_path, os.path.dirname(result.path)
                                            )
                                        )
                                    ),
                                )

                # write report html
                zipout.writestr(str(folder.joinpath("report.html")), rendered)
        else:
            with open(self.settings.path, "w", encoding="utf-8") as htmlout:
                htmlout.write(rendered)
