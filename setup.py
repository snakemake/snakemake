__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from setuptools.command.test import test as TestCommand
import sys


if sys.version_info < (3, 3):
    print("At least Python 3.3 is required.\n", file=sys.stderr)
    exit(1)


try:
    from setuptools import setup
except ImportError:
    print("Please install setuptools before installing snakemake.",
          file=sys.stderr)
    exit(1)


# load version info
exec(open("snakemake/version.py").read())


class NoseTestCommand(TestCommand):
    def run_tests(self):
        # Run nose ensuring that argv simulates running nosetests directly
        import nose
        nose.run_exit(argv=['nosetests'])


setup(
    name='snakemake',
    version=__version__,
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
    url='http://snakemake.bitbucket.org',
    packages=['snakemake', 'snakemake.remote'],
    entry_points={
        "console_scripts":
        ["snakemake = snakemake:main",
         "snakemake-bash-completion = snakemake:bash_completion"]
    },
    package_data={'': ['*.css', '*.sh', '*.html']},
    tests_require=['rpy2', 'httpretty==0.8.10', 'docutils', 'nose>=1.3', 'boto>=2.38.0', 'filechunkio>=1.6', 
                     'moto>=0.4.14', 'ftputil>=3.2', 'pysftp>=0.2.8', 'requests>=2.8.1', 'dropbox>=5.2'],
    cmdclass={'test': NoseTestCommand},
    classifiers=
    ["Development Status :: 5 - Production/Stable", "Environment :: Console",
     "Intended Audience :: Science/Research",
     "License :: OSI Approved :: MIT License", "Natural Language :: English",
     "Programming Language :: Python :: 3",
     "Topic :: Scientific/Engineering :: Bio-Informatics"])
