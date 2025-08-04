.. _snakefiles-testing:

===================================
Automatically generating unit tests
===================================

Snakemake can automatically generate unit tests from a workflow that has already been successfully executed.
By running

.. code-block:: bash

    snakemake --generate-unit-tests

Snakemake is instructed to take one representative job for each rule and copy its input files to a hidden folder ``.tests/unit``,
along with generating test cases for Pytest_.
Pytest_ tests can can be run as:

.. code-block:: bash

    pytest .tests/unit/

or, optionally, if you want to use a local conda cache and disable pytest caching:

.. code-block:: bash

    pytest -p no:cacheprovider .tests/unit/ --conda-prefix /path/to/cache/conda/

Each auto-generated unit test is stored in a file ``.tests/unit/test_<rulename>.py``, and executes just the one representative job of the respective rule.
After successful execution of the job, it will compare the obtained results with those that have been present when running ``snakemake --generate-unit-tests``.
By default, the comparison happens byte by byte (using ``cmp/zcmp/bzcmp``). This behavior can be overwritten by modifying the test file.

NOTE: Importantly, such unit tests shall not be generated from big data, as they should usually be finished in a few seconds.
Furthermore, it makes sense to store the generated unit tests in version control (e.g. git), such that huge files are also not recommended.
Instead, we suggest to first execute the workflow that shall be tested with some kind of small dummy datasets while keeping all temp files (``--notemp``),
and then use the results thereof to generate the unit tests.
The small dummy datasets can in addition be used to generate an integration test, that could e.g. be stored under ``.tests/integration``, next to the unit tests.

.. _Pytest: https://pytest.org
