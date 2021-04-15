# -*- coding: UTF-8 -*-

from __future__ import print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import sys
import versioneer


if sys.version_info < (3, 5):
    print("At least Python 3.5 is required for Snakemake.\n", file=sys.stderr)
    exit(1)


try:
    from setuptools import setup
except ImportError:
    print("Please install setuptools before installing snakemake.", file=sys.stderr)
    exit(1)


setup(
    name="snakemake",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Johannes Köster",
    author_email="johannes.koester@tu-dortmund.de",
    description="Snakemake is a workflow management system that aims to reduce the complexity "
    "of creating workflows by providing a fast and comfortable execution environment, "
    "together with a clean and modern specification language in python style. "
    "Snakemake workflows are essentially Python scripts extended by declarative "
    "code to define rules. Rules describe how to create output files from input files.",
    zip_safe=False,
    license="MIT",
    url="https://snakemake.readthedocs.io",
    packages=[
        "snakemake",
        "snakemake.remote",
        "snakemake.report",
        "snakemake.caching",
        "snakemake.deployment",
        "snakemake.linting",
        "snakemake.executors",
        "snakemake.unit_tests",
        "snakemake.unit_tests.templates"
    ],
    entry_points={
        "console_scripts": [
            "snakemake = snakemake:main",
            "snakemake-bash-completion = snakemake:bash_completion",
        ]
    },
    package_data={"": ["*.css", "*.sh", "*.html", "*.jinja2"]},
    install_requires=[
        "wrapt",
        "requests",
        "ratelimiter",
        "pyyaml",
        "configargparse",
        "appdirs",
        "datrie",
        "jsonschema",
        "docutils",
        "gitpython",
        "psutil",
        "nbformat",
        "toposort",
        "pulp >=2.0",
        "smart_open",
        "filelock",
        "stopit",
    ],
    extras_require={
        "reports": ["jinja2", "networkx", "pygments", "pygraphviz"],
        "messaging": ["slacker"],
        "google-cloud": [
            "oauth2client",
            "google-crc32c",
            "google-api-python-client",
            "google-cloud-storage",
        ],
        "pep": [
            "peppy",
            "eido",
        ]
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
