# -*- coding: UTF-8 -*-

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"

import sys
import versioneer


if sys.version_info < (3, 5):
    print("At least Python 3.5 is required.\n", file=sys.stderr)
    exit(1)


try:
    from setuptools import setup
except ImportError:
    print("Please install setuptools before installing snakemake.",
          file=sys.stderr)
    exit(1)


setup(
    name='snakemake',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author='Johannes Köster',
    author_email='johannes.koester@tu-dortmund.de',
    description=
    'Snakemake is a workflow management system that aims to reduce the complexity '
    'of creating workflows by providing a fast and comfortable execution environment, '
    'together with a clean and modern specification language in python style. '
    'Snakemake workflows are essentially Python scripts extended by declarative '
    'code to define rules. Rules describe how to create output files from input files.',
    zip_safe=False,
    license='MIT',
    url='http://snakemake.bitbucket.io',
    packages=['snakemake', 'snakemake.remote', 'snakemake.report'],
    entry_points={
        "console_scripts":
        ["snakemake = snakemake:main",
         "snakemake-bash-completion = snakemake:bash_completion"]
    },
    package_data={'': ['*.css', '*.sh', '*.html']},
    install_requires=['wrapt', 'requests', 'ratelimiter', 'pyyaml',
                      'configargparse', 'appdirs', 'datrie', 'jsonschema',
                      'docutils', 'gitpython'],
    extras_require={"reports": ['jinja2', 'networkx']},
    classifiers=
    ["Development Status :: 5 - Production/Stable", "Environment :: Console",
     "Intended Audience :: Science/Research",
     "License :: OSI Approved :: MIT License", "Natural Language :: English",
     "Programming Language :: Python :: 3.5",
     "Topic :: Scientific/Engineering :: Bio-Informatics"])
