.. _project_info-contributing:

============
Contributing
============

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

You can contribute in many ways:


----------------------
Types of Contributions
----------------------


Report Bugs
===========

Report bugs at https://github.com/snakemake/snakemake/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.


Fix Bugs
========

Look through the Github issues for bugs.
If you want to start working on a bug then please write short message on the issue tracker to prevent duplicate work.


Implement Features
==================

Look through the Github issues for features.
If you want to start working on an issue then please write short message on the issue tracker to prevent duplicate work.

Contributing a new cluster or cloud execution backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execution backends are added by implementing a so-called ``Executor``.
All executors are located in `snakemake/executors.py <https://github.com/snakemake/snakemake/tree/master/snakemake/executors.py>`_.
In order to implement a new executor, you have to inherit from the class ``ClusterExecutor``.
Below you find a skeleton

.. code-block:: python

    class SkeletonExecutor(ClusterExecutor):
        def __init__(self, workflow, dag, cores,
                 jobname="snakejob.{name}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 cluster_config=None,
                 local_input=None,
                 restart_times=None,
                 exec_job=None,
                 assume_shared_fs=True,
                 max_status_checks_per_second=1):

            # overwrite the command to execute a single snakemake job if necessary
            # exec_job = "..."

            super().__init__(workflow, dag, None,
                             jobname=jobname,
                             printreason=printreason,
                             quiet=quiet,
                             printshellcmds=printshellcmds,
                             latency_wait=latency_wait,
                             cluster_config=cluster_config,
                             local_input=local_input,
                             restart_times=restart_times,
                             exec_job=exec_job,
                             assume_shared_fs=False,
                             max_status_checks_per_second=10)

            # add additional attributes

        def shutdown(self):
            # perform additional steps on shutdown if necessary
            super().shutdown()

        def cancel(self):
            for job in self.active_jobs:
                # cancel active jobs here
            self.shutdown()
        
        def run_jobs(self, jobs, callback=None, submit_callback=None, error_callback=None):
            """Run a list of jobs that is ready at a given point in time.

            By default, this method just runs each job individually.
            This behavior is inherited and therefore this method can be removed from the skeleton if the
            default behavior is intended.
            This method can be overwritten to submit many jobs in a more efficient way than one-by-one.

            Note that in any case, for each job, the callback functions have to be called individually!
            """
            for job in jobs:
                self.run(
                    job,
                    callback=callback,
                    submit_callback=submit_callback,
                    error_callback=error_callback,
                )

        def run(self, job,
                callback=None,
                submit_callback=None,
                error_callback=None):
            """Run an individual job or a job group.
            """

            super()._run(job)
            # obtain job execution command
            exec_job = self.format_job(
                self.exec_job, job, _quote_all=True,
                use_threads="--force-use-threads" if not job.is_group() else "")

            # submit job here, and obtain job ids from the backend

            # register job as active, using your own namedtuple.
            # The namedtuple must at least contain the attributes
            # job, jobid, callback, error_callback.
            self.active_jobs.append(MyJob(
                job, jobid, callback, error_callback))

        def _wait_for_jobs(self):
            # busy wait on job completion
            # This is only needed if your backend does not allow to use callbacks
            # for obtaining job status.
            while True:
                # always use self.lock to avoid race conditions
                with self.lock:
                    if not self.wait:
                        return
                    active_jobs = self.active_jobs
                    self.active_jobs = list()
                    still_running = list()
                for j in active_jobs:
                    # use self.status_rate_limiter to avoid too many API calls.
                    with self.status_rate_limiter:

                        # Retrieve status of job j from your backend via j.jobid
                        # Handle completion and errors, calling either j.callback(j.job)
                        # or j.error_callback(j.job)
                        # In case of error, add job j to still_running.
                with self.lock:
                    self.active_jobs.extend(still_running)
                sleep()


Write Documentation
===================

Snakemake could always use more documentation, whether as part of the official vcfpy docs, in docstrings, or even on the web in blog posts, articles, and such.

Snakemake uses `Sphinx <https://sphinx-doc.org>`_ for the user manual (that you are currently reading).
See `project_info-doc_guidelines` on how the documentation reStructuredText is used.


Submit Feedback
===============

The best way to send feedback is to file an issue at https://github.com/snakemake/snakemake/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :)

-----------------------
Pull Request Guidelines
-----------------------

To update the documentation, fix bugs or add new features you need to create a Pull Request
. A PR is a change you make to your local copy of the code for us to review and potentially integrate into the code base.

To create a Pull Request you need to do these steps:

1. Create a Github account.
2. Fork the repository.
3. Clone your fork locally.
4. Go to the created snakemake folder with :code:`cd snakemake`.
5. Create a new branch with :code:`git checkout -b <descriptive_branch_name>`.
6. Make your changes to the code or documentation.
7. Run :code:`git add .` to add all the changed files to the commit (to see what files will be added you can run :code:`git add . --dry-run`).
8. To commit the added files use :code:`git commit`. (This will open a command line editor to write a commit message. These should have a descriptive 80 line header, followed by an empty line, and then a description of what you did and why. To use your command line text editor of choice use (for example) :code:`export GIT_EDITOR=vim` before running :code:`git commit`).
9. Now you can push your changes to your Github copy of Snakemake by running :code:`git push origin <descriptive_branch_name>`.
10. If you now go to the webpage for your Github copy of Snakemake you should see a link in the sidebar called "Create Pull Request".
11. Now you need to choose your PR from the menu and click the "Create pull request" button. Be sure to change the pull request target branch to <descriptive_branch_name>!

If you want to create more pull requests, first run :code:`git checkout master` and then start at step 5. with a new branch name.

Feel free to ask questions about this if you want to contribute to Snakemake :)

------------------
Testing Guidelines
------------------

To ensure that you do not introduce bugs into Snakemake, you should test your code thouroughly.

To have integration tests run automatically when commiting code changes to Github, you need to sign up on wercker.com and register a user.

The easiest way to run your development version of Snakemake is perhaps to go to the folder containing your local copy of Snakemake and call:

.. code-block:: console

    $ conda env create -f environment.yml -n snakemake-testing
    $ conda activate snakemake-testing
    $ pip install -e .

This will make your development version of Snakemake the one called when running snakemake. You do not need to run this command after each time you make code changes.

From the base snakemake folder you call :code:`nosetests` to run all the tests, or choose one specific test. For this to work, Nose (the testing framework we use) can be installed to the conda environment using pip:

.. code-block:: console

   $ pip install nose
   $ nosetests
   $ nosetests tests.tests:test_log_input

If you introduce a new feature you should add a new test to the tests directory. See the folder for examples.

.. project_info-doc_guidelines:

------------------------
Documentation Guidelines
------------------------

For the documentation, please adhere to the following guidelines:

- Put each sentence on its own line, this makes tracking changes through Git SCM easier.
- Provide hyperlink targets, at least for the first two section levels.
  For this, use the format ``<document_part>-<section_name>``, e.g., ``project_info-doc_guidelines``.
- Use the section structure from below.

::

    .. document_part-heading_1:

    =========
    Heading 1
    =========


    .. document_part-heading_2:

    ---------
    Heading 2
    ---------


    .. document_part-heading_3:

    Heading 3
    =========


    .. document_part-heading_4:

    Heading 4
    ---------


    .. document_part-heading_5:

    Heading 5
    ~~~~~~~~~


    .. document_part-heading_6:

    Heading 6
    :::::::::

.. _doc_setup:

-------------------
Documentation Setup
-------------------

For building the documentation, you have to install the Sphinx.
If you have already installed Conda, all you need to do is to create a
Snakemake development environment via

.. code-block:: console

    $ git clone git@github.com:snakemake/snakemake.git
    $ cd snakemake
    $ conda env create -f environment.yml -n snakemake

Then, the docs can be built with

.. code-block:: console

    $ conda activate snakemake
    $ cd docs
    $ make html
    $ make clean && make html  # force rebuild

Alternatively, you can use virtualenv.
The following assumes you have a working Python 3 setup.

.. code-block:: console

    $ git clone git@github.org:snakemake/snakemake.git
    $ cd snakemake/docs
    $ virtualenv -p python3 .venv
    $ source .venv/bin/activate
    $ pip install --upgrade -r requirements.txt

Afterwards, the docs can be built with

.. code-block:: console

    $ source .venv/bin/activate
    $ make html  # rebuild for changed files only
    $ make clean && make html  # force rebuild
