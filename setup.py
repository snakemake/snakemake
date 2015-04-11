# coding: utf-8

import sys

if sys.version_info < (3,2):
    print("At least Python 3.2 is required.\n", file=sys.stderr)
    exit(1)


try:
    from setuptools import setup
except ImportError:
    print("Please install setuptools before installing snakemake.", file=sys.stderr)
    exit(1)


# load version info
exec(open("snakemake/version.py").read())


setup(
    name='snakemake',
    version=__version__,
    author='Johannes Köster',
    author_email='johannes.koester@tu-dortmund.de',
    description='Build systems like make are frequently used to create complicated workflows, e.g. in bioinformatics. This project aims to reduce the complexity of creating workflows by providing a clean and modern domain specific language (DSL) in python style, together with a fast and comfortable execution environment.',
    zip_safe=False,
    license='MIT',
    url='http://snakemake.googlecode.com',
    packages=['snakemake'],
    entry_points={
        "console_scripts": ["snakemake = snakemake:main", "snakemake-bash-completion = snakemake:bash_completion"]
    },
    package_data={'': ['*.css', '*.sh', '*.html']},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
