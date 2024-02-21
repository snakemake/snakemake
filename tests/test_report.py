__authors__ = ["Tobias Marschall", "Marcel Martin", "Johannes Köster"]
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import sys


sys.path.insert(0, os.path.dirname(__file__))

from .common import *
def test_report():
    run(
        dpath("test_report"),
        report="report.html",
        report_stylesheet="custom-stylesheet.css",
        check_md5=False,
    )