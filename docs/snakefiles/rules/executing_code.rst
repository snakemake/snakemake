.. _snakefiles-plain-python-rules:
==================
Executing Code
==================

Plain python rules
------------------

Instead of a shell command, a rule can run some python code to generate the output.
It is highly advisable to limit such code to a few lines.
Otherwise, use Snakemake's :ref:`script support <snakefiles-external_scripts>`.

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile",
        output:
            "path/to/outputfile",
            somename="path/to/another/outputfile",
        run:
            for f in input:
                ...
                with open(output[0], "w") as out:
                    out.write(...)
            with open(output.somename, "w") as out:
                out.write(...)

As can be seen, instead of accessing input and output as a whole, we can also access by index (``output[0]``) or by keyword (``output.somename``).
Note that, when adding keywords or names for input or output files, their order won't be preserved when accessing them as a whole via e.g. ``{output}`` in a shell command.

Shell commands like above can also be invoked inside a python based rule, via the function ``shell`` that takes a string with the command and allows the same formatting like in the rule above, e.g.:

.. code-block:: python

    shell("somecommand {output.somename}")

Further, this combination of python and shell commands allows us to iterate over the output of the shell command, e.g.:

.. code-block:: python

    for line in shell("somecommand {output.somename}", iterable=True):
        ... # do something in python

.. _snakefiles-external_scripts:

External scripts
----------------

A rule can also point to an external script instead of a shell command or inline Python code, e.g.

Python
~~~~~~

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        script:
            "scripts/script.py"

.. note::

    It is possible to refer to wildcards and params in the script path, e.g. by specifying ``"scripts/{params.scriptname}.py"`` or ``"scripts/{wildcards.scriptname}.py"``.

The script path is always relative to the Snakefile containing the directive (in contrast to the input and output file paths, which are relative to the working directory).
It is recommended to put all scripts into a subfolder ``scripts`` as above.
Inside the script, you have access to an object ``snakemake`` that provides access to the same objects that are available in the ``run`` and ``shell`` directives (input, output, params, wildcards, log, threads, resources, config), e.g. you can use ``snakemake.input[0]`` to access the first input file of above rule.
It is also possible to explicitly import the snakemake object in the script like ``from snakemake.script import snakemake`` to enable code completion, linting and type checking your python code in IDEs.

An example external Python script could look like this:

.. code-block:: python

    def do_something(data_path, out_path, threads, myparam):
        # python code

    do_something(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config["myparam"])

or using the explicit import:

.. code-block:: python

    from snakemake.script import snakemake

    def do_something(data_path, out_path, threads, myparam):
        # python code

    do_something(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config["myparam"])

You can use the Python debugger from within the script if you invoke Snakemake with ``--debug``.

Xonsh_
""""""

.. _Xonsh: https://xon.sh

Because Xonsh is a superset of Python, you can use a Xonsh script as you would a Python script, but with all the additional shell primitives that Xonsh provides.

For example, with this rule:

.. code-block:: python

    rule get_variants_in_genes:
        input:
            vcf="input.vcf",
            gene_locations="genes.bed",
        output:
            "output.tsv"
        conda:
            "envs/variant_calling.yaml"
        log:
            "logs/get_variants_in_genes.log"
        script:
            "scripts/get_variants_in_genes.xsh"

the Xonsh script might look like this:

.. code-block::

    $XONSH_TRACEBACK_LOGFILE = snakemake.log[0]

    annotations = ", ".join(
        f'ANN["{field}"]' for field in ["Consequence", "SYMBOL", "Feature", "BIOTYPE"]
    )

    bcftools view -R @(snakemake.input.gene_locations) @(snakemake.input.vcf) \
    | vembrane table --output @(snakemake.output[0]) @(f'CHROM, POS, ID, {annotations}')


Hy_
"""

.. _Hy: https://hylang.org/

Hy allows you to interact with Python using a Lisp-like syntax.

For example, with this rule:

.. code-block:: python

    rule get_sum_of_odd_numbers:
        input:
            "list_of_numbers.txt"
        output:
            results_file="sum_of_odd_numbers.txt"
        conda:
            "envs/hy.yaml"
        script:
            "scripts/sum_odd_numbers.hy"

the Hy script might look like this:

.. code-block:: hy

    (require hyrule [-> ->>])

    (defn is-odd? [n] (!= (% n 2) 0))

    (setv result
          (->> (get snakemake.input 0)
               open
               .readlines
               (map int)
               (filter is-odd?)
               sum))

    (with [f
           (-> (get snakemake.output "results_file")
               (open "w"))]
      (print result :file f))



R and R Markdown
~~~~~~~~~~~~~~~~

Apart from Python scripts, this mechanism also allows you to integrate R_ and R Markdown_ scripts with Snakemake, e.g.

.. _R: https://www.r-project.org
.. _Markdown: https://rmarkdown.rstudio.com

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        script:
            "scripts/script.R"

In the R script, an S4 object named ``snakemake`` analogous to the Python case above is available and allows access to input and output files and other parameters. Here the syntax follows that of S4 classes with attributes that are R lists, e.g. we can access the first input file with ``snakemake@input[[1]]`` (note that the first file does not have index ``0`` here, because R starts counting from ``1``). Named input and output files can be accessed in the same way, by just providing the name instead of an index, e.g. ``snakemake@input[["myfile"]]``.

An equivalent script (:ref:`to the Python one above <Python>`) written in R would look like this:

.. code-block:: r

    do_something <- function(data_path, out_path, threads, myparam) {
        # R code
    }

    do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config[["myparam"]])


To debug R scripts, you can save the workspace with ``save.image()``, and invoke R after Snakemake has terminated. Then you can use the usual R debugging facilities while having access to the ``snakemake`` variable.
It is best practice to wrap the actual code into a separate function. This increases the portability if the code shall be invoked outside of Snakemake or from a different rule.
A convenience method, ``snakemake@source()``, acts as a wrapper for the normal R ``source()`` function, and can be used to source files relative to the original script directory.

An R Markdown file can be integrated in the same way as R and Python scripts, but only a single output (html) file can be used:

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/report.html",
        script:
            "path/to/report.Rmd"

In the R Markdown file you can insert output from a R command, and access variables stored in the S4 object named ``snakemake``

.. code-block:: R

    ---
    title: "Test Report"
    author:
        - "Your Name"
    date: "`r format(Sys.time(), '%d %B, %Y')`"
    params:
       rmd: "report.Rmd"
    output:
      html_document:
      highlight: tango
      number_sections: no
      theme: default
      toc: yes
      toc_depth: 3
      toc_float:
        collapsed: no
        smooth_scroll: yes
    ---

    ## R Markdown

    This is an R Markdown document.

    Test include from snakemake `r snakemake@input`.

    ## Source
    <a download="report.Rmd" href="`r base64enc::dataURI(file = params$rmd, mime = 'text/rmd', encoding = 'base64')`">R Markdown source file (to produce this document)</a>

A link to the R Markdown document with the snakemake object can be inserted. Therefore a variable called ``rmd`` needs to be added to the ``params`` section in the header of the ``report.Rmd`` file. The generated R Markdown file with snakemake object will be saved in the file specified in this ``rmd`` variable. This file can be embedded into the HTML document using base64 encoding and a link can be inserted as shown in the example above.
Also other input and output files can be embedded in this way to make a portable report. Note that the above method with a data URI only works for small files. An experimental technology to embed larger files is using Javascript Blob `object <https://developer.mozilla.org/en-US/docs/Web/API/Blob>`_.

Julia_
~~~~~~

.. _Julia: https://julialang.org

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        script:
            "path/to/script.jl"

In the Julia_ script, a ``snakemake`` object is available, which can be accessed similar to the :ref:`Python case <Python>`, with the only difference that you have to index from 1 instead of 0.

Rust_
~~~~~

.. _Rust: https://www.rust-lang.org/

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile",
            named_input="path/to/named/inputfile",
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        params:
            seed=4
        conda:
            "rust.yaml"
        log:
            stdout="path/to/stdout.log",
            stderr="path/to/stderr.log",
        script:
            "path/to/script.rs"

The ability to execute Rust scripts is facilitated by |rust-script|_.
As such, the script must be a valid ``rust-script`` script and ``rust-script``
(plus OpenSSL and a C compiler toolchain, provided by Conda packages ``openssl``, ``c-compiler``, ``pkg-config``)
must be available in the environment the rule is run in.
The minimum required ``rust-script`` version is 0.35.0, so in the example above, the contents of ``rust.yaml`` might look like this:

.. code-block:: yaml

    channels:
      - conda-forge
      - bioconda
    dependencies:
      - rust-script>=0.35.0
      - openssl
      - c-compiler
      - pkg-config

Some example scripts can be found in the
`tests directory <https://github.com/snakemake/snakemake/tree/main/tests/test_script/scripts>`_.

In the Rust script, a ``snakemake`` instance is available, which is automatically generated from the python snakemake object using |json_typegen|_.
It usually looks like this:

.. code-block:: rust

    pub struct Snakemake {
        input: Input,
        output: Output,
        params: Params,
        wildcards: Wildcards,
        threads: u64,
        log: Log,
        resources: Resources,
        config: Config,
        rulename: String,
        bench_iteration: Option<usize>,
        scriptdir: String,
    }

Any named parameter is translated to a corresponding ``field_name: Type``, such that ``params.seed`` from the example above can be accessed just like in python, i.e.:

.. code-block:: rust

    let seed = snakemake.params.seed;
    assert_eq!(seed, 4);

Positional arguments for ``input``, ``output``, ``log`` and ``wildcards`` can be accessed by index and iterated over:

.. code-block:: rust

    let input = &snakemake.input;

    // Input implements Index<usize>
    let inputfile = input[0];
    assert_eq!(inputfile, "path/to/inputfile");

    // Input implements IntoIterator
    //
    // prints
    // > 'path/to/inputfile'
    // > 'path/to/other/inputfile'
    for f in input {
        println!("> '{}'", &f);
    }


It is also possible to redirect ``stdout`` and ``stderr``:

.. code-block:: rust

    println!("This will NOT be written to path/to/stdout.log");
    // redirect stdout to "path/to/stdout.log"
    let _stdout_redirect = snakemake.redirect_stdout(snakemake.log.stdout)?;
    println!("This will be written to path/to/stdout.log");

    // redirect stderr to "path/to/stderr.log"
    let _stderr_redirect = snakemake.redirect_stderr(snakemake.log.stderr)?;
    eprintln!("This will be written to path/to/stderr.log");
    drop(_stderr_redirect);
    eprintln!("This will NOT be written to path/to/stderr.log");

Redirection of stdout/stderr is only "active" as long as the returned ``Redirect`` instance is alive; in order to stop redirecting, drop the respective instance.

In order to work, rust-script support for snakemake has some dependencies enabled by default:

#. ``anyhow=1``, for its ``Result`` type
#. ``gag=1``, to enable stdout/stderr redirects
#. ``json_typegen=0.6``, for generating rust structs from a json representation of the snakemake object
#. ``lazy_static=1.4``, to make a ``snakemake`` instance easily accessible
#. ``serde=1.0``, explicit dependency of ``json_typegen``
#. ``serde_derive=1.0``, explicit dependency of ``json_typegen``
#. ``serde_json=1.0``, explicit dependency of ``json_typegen``

If your script uses any of these packages, you do not need to ``use`` them in your script. Trying to ``use`` them will cause a compilation error.

.. |rust-script| replace:: ``rust-script``
.. _rust-script: https://rust-script.org/
.. |json_typegen| replace:: ``json_typegen``
.. _json_typegen: https://github.com/evestera/json_typegen


Bash
~~~~

Bash scripts work much the same as the other script languages above, but with some important differences. Access to the
rule's directives is provided through the use of `associative arrays <arrays_>`_ - **requiring Bash version 4.0 or greater**.
One "limitation" of associative arrays is they cannot be nested. As such, the following rule directives are found in a separate
variable, named as ``snakemake_<directive>``:

* ``input``
* ``output``
* ``log``
* ``wildcards``
* ``resources``
* ``params``
* ``config``

Access to the ``input`` directive is facilitated through the bash associative array named ``snakemake_input``. The
remaining directives can be found in the variable ``snakemake``.

.. note::

    As arrays cannot be nested in Bash, use of python's ``dict`` in directives is not supported. So, adding a ``params`` key of ``data={"foo": "bar"}`` will not be reflected - ``${snakemake_params[data]}`` actually only returns ``"foo"``.

Bash Example 1
""""""""""""""

.. code-block:: python

    rule align:
        input:
            "{sample}.fq",
            reference="ref.fa",
        output:
            "{sample}.sam"
        params:
            opts="-a -x map-ont",
        threads: 4
        log:
            "align/{sample}.log"
        conda:
            "envs/align.yaml"
        script:
            "scripts/align.sh"



``align.sh``

.. code-block:: bash

    #!/usr/bin/env bash

    echo "Aligning sample ${snakemake_wildcards[sample]} with minimap2" 2> "${snakemake_log[0]}"

    minimap2 ${snakemake_params[opts]} -t ${snakemake[threads]} "${snakemake_input[reference]}" \
        "${snakemake_input[0]}" > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"


If you don't add a shebang, the default ``#!/usr/bin/env bash`` will be inserted for you. A tutorial on how to use
associative arrays can be found `here <https://www.xmodulo.com/key-value-dictionary-bash.html>`__.

You may also have noticed the mixed use of double-quotes when accessing some variables. It is generally good practice in
Bash to double-quote variables for which you want to `prevent word splitting <split_>`_; generally, you will want to
double-quote any variable that could contain a file name. However, `in some cases <exception_>`_, word splitting *is* desired,
such as ``${snakemake_params[opts]}`` in the above example.

Bash Example 2
""""""""""""""

.. code-block:: python

    rule align:
        input:
            reads=["{sample}_R1.fq", "{sample}_R2.fq]"],
            reference="ref.fa",
        output:
            "{sample}.sam"
        params:
            opts="-M",
        threads: 4
        log:
            "align/{sample}.log"
        conda:
            "envs/align.yaml"
        script:
            "scripts/align.sh"


In this example, the ``input`` variable ``reads``, which is a python list, actually gets stored as a space-separated string
in Bash because, you guessed it, you can't nest arrays in Bash! So in order to access the individual members, we turn the
string into an array; allowing us to access individual elements of the list/array. See `this stackoverflow question <so_>`_ for other solutions.

``align.sh``

.. code-block:: bash

    #!/usr/bin/env bash

    exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

    reads=(${snakemake_input[reads]})  # don't double-quote this - we want word splitting

    r1="${reads[0]}"
    r2="${reads[1]}"

    bwa index "${snakemake_input[reference]}"
    bwa mem ${snakemake_params[opts]} -t ${snakemake[threads]} \
        "${snakemake_input[reference]}" "$r1" "$r2" > "${snakemake_output[0]}"

If, in the above example, the fastq reads were not in a named variable, but were instead just a list, they would be available
as ``"${snakemake_input[0]}"`` and ``"${snakemake_input[1]}"``.

.. _arrays: https://www.gnu.org/software/bash/manual/html_node/Arrays.html#Arrays
.. _split: https://github.com/koalaman/shellcheck/wiki/SC2046
.. _exception: https://github.com/koalaman/shellcheck/wiki/SC2046#exceptions
.. _so: https://stackoverflow.com/q/1469849/5299417

----

For technical reasons, scripts are executed in ``.snakemake/scripts``. The original script directory is available as ``scriptdir`` in the ``snakemake`` object.

.. _snakefiles_notebook-integration:

Jupyter notebook integration
----------------------------

Instead of plain scripts (see above), one can integrate Jupyter_ Notebooks.
This enables the interactive development of data analysis components (e.g. for plotting).
Integration works as follows (note the use of `notebook:` instead of `script:`):

.. _Jupyter: https://jupyter.org/

.. code-block:: python

    rule hello:
        output:
            "test.txt"
        log:
            # optional path to the processed notebook
            notebook="logs/notebooks/processed_notebook.ipynb"
        notebook:
            "notebooks/hello.py.ipynb"

.. note::

    Consider Jupyter notebook integration as a way to get the best of both worlds.
    A modular, readable workflow definition with Snakemake, and the ability to quickly explore and plot data with Jupyter.
    The benefit will be maximal when integrating many small notebooks that each do a particular job, hence allowing to get away from large monolithic, and therefore unreadable notebooks.

It is recommended to prefix the ``.ipynb`` suffix with either ``.py`` or ``.r`` to indicate the notebook language.
In the notebook, a snakemake object is available, which can be accessed in the same way as the with :ref:`script integration <snakefiles-external_scripts>`.
In other words, you have access to input files via ``snakemake.input`` (in the Python case) and ``snakemake@input`` (in the R case) etc..
Optionally it is possible to automatically store the processed notebook.
This can be achieved by adding a named logfile ``notebook=...`` to the ``log`` directive.

.. note::

    It is possible to refer to wildcards and params in the notebook path, e.g. by specifying ``"notebook/{params.name}.py"`` or ``"notebook/{wildcards.name}.py"``.

Normally, notebooks are executed headlessly (without a Jupyter interface being presented to you).
This is achieved with Papermill_ if that is installed in your software environment,
or `nbconvert`_ otherwise.
The latter will be installed automatically along with Jupyter, but will not output
an executed (logfile) notebook until the entire execution is complete, and won't output
a notebook if execution encounters an error.

.. _Papermill: https://github.com/nteract/papermill

.. _nbconvert: https://nbconvert.readthedocs.io/en/latest/

In order to simplify the coding of notebooks given the automatically inserted ``snakemake`` object, Snakemake provides an interactive edit mode for notebook rules.
Let us assume you have written above rule, but the notebook does not yet exist.
By running

.. code-block:: console

    snakemake --cores 1 --edit-notebook test.txt

you instruct Snakemake to allow interactive editing of the notebook needed to create the file ``test.txt``.
Snakemake will run all dependencies of the notebook rule, such that all input files are present.
Then, it will start a jupyter notebook server with an empty draft of the notebook, in which you can interactively program everything needed for this particular step.
Once done, you should save the notebook from the jupyter web interface, go to the jupyter dashboard and hit the ``Quit`` button on the top right in order to shut down the jupyter server.
Snakemake will detect that the server is closed and automatically store the drafted notebook into the path given in the rule (here ``hello.py.ipynb``).
If the notebook already exists, above procedure can be used to easily modify it.
Note that Snakemake requires local execution for the notebook edit mode.
On a cluster or the cloud, you can generate all dependencies of the notebook rule via

.. code-block:: console

    snakemake --cluster ... --jobs 100 --until test.txt

Then, the notebook rule can easily be executed locally.
An demo of the entire interactive editing process can be found by clicking below:

.. image:: images/snakemake-notebook-demo.gif
    :scale: 20%
    :alt: Notebook integration demo
    :align: center

Finally, it is advisable to combine the ``notebook`` directive with the ``conda`` directive (see :ref:`integrated_package_management`) in order to define a software stack to use.
At least, this software stack should contain jupyter and the language to use (e.g. Python or R).
For the above case, this means

.. code-block:: python

    rule hello:
        output:
            "test.txt"
        conda:
            "envs/hello.yaml"
        notebook:
            "notebooks/hello.py.ipynb"

with

.. code-block:: yaml

    channels:
      - conda-forge
    dependencies:
      - python =3.8
      - jupyter =1.0
      - jupyterlab_code_formatter =1.4

The last dependency is advisable in order to enable autoformatting of notebook cells when editing.
When using other languages than Python in the notebook, one needs to additionally add the respective kernel, e.g. ``r-irkernel`` for R support.

When using an IDE with built-in Jupyter support, an alternative to ``--edit-notebook`` is ``--draft-notebook``.
Instead of firing up a notebook server, ``--draft-notebook`` just creates a skeleton notebook for editing within the IDE.
In addition, it prints instructions for configuring the IDE's notebook environment to use the interpreter from the
Conda environment defined in the corresponding rule.
For example, running

.. code-block:: console

    snakemake --cores 1 --draft-notebook test.txt --software-deployment-method conda

or the short form

.. code-block:: console

    snakemake -c 1 --draft-notebook test.txt --sdm conda

will generate skeleton code in ``notebooks/hello.py.ipynb`` and additionally print instructions on how to open and execute the notebook in VSCode.



.. _shell_settings:

Shell settings
--------------

By default, Snakemake uses the ``bash`` shell.
This can be overridden in two ways.
First, by globally setting the shell executable (e.g. to zsh) via

.. code-block:: python

    shell.executable("/bin/zsh")

Note that this is usually not recommended, as it requires others who want to use the workflow to have that shell installed.
Second, by setting the shell executable via the :ref:`resources directive <snakefiles-resources>` of a rule, e.g.

.. code-block:: python

    rule a:
        input: ...
        output: ...
        resources:
            shell_exec="zsh"
        shell:
            "echo 'hello world' > {output}"

This can be particularly important in case you use a :ref:`container image <apptainer>` for the rule which does not contain bash, e.g.

.. code-block:: python

    rule a:
        output:
            "test.out"
        resources:
            shell_exec="sh"
        # image does not have bash, hence this would fail if shell_exec is not set to sh
        container: "docker://busybox:1.33"
        shell:
            "echo 'hello world' > {output}"

Shell behavior
^^^^^^^^^^^^^^

In case of bash shell, Snakemake always uses the so-called `strict mode <http://redsymbol.net/articles/unofficial-bash-strict-mode/>`_.
For individual rules, you can deactivate aspects of the strict mode by unsetting them at the beginning of the shell command.
Further, it is possible to set global prefixes and suffixes for all shell commands via

.. code-block:: python

    shell.prefix("some prefix command;")
    shell.suffix("; some suffix command")

anywhere in your snakefile (preferably at the beginning for clarity).
This can sometimes be useful for debugging, but is not recommended for production workflows and releases because it might hamper reproducibility and readability.



.. _snakefiles_mpi_support:

MPI support
-----------

Some highly parallel programs or scripts implement the `message passing interface (MPI)) <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_, which enables a program to span work across multiple compute nodes on a compute cluster (where a node is an individual machine, which will usually have multiple CPUs nowadays).
Let us assume, we have such a program that can parallelize using the MPI and its name is ``calc-pi-mpi``.
To run such a program, the user will usually launch it via the `mpirun <https://docs.open-mpi.org/en/v5.0.x/launching-apps/quickstart.html>`_ command from Open MPI.
But because it can make sense to use another MPI launch command in some circumstances, we recommend the following pattern for a snakemake rule using an MPI-parallelized tool:

.. code-block:: python

  rule calc_pi:
    output:
        "pi.calc",
    log:
        "logs/calc_pi.log",
    resources:
    resources:
        tasks=10,
        mpi="mpirun",
    shell:
        "{resources.mpi} -n {resources.tasks} calc-pi-mpi 10 > {output} 2> {log}"

Here, you provide the MPI wrapper command used to launch the program under ``resources: mpi=``.
This enables users to override this command if their execution environment requires this, via providing the ``mpi`` resource in :ref:`executing-profiles` or via the command line option |set-resources|_.
While ``mpirun`` `should work in most compute environments, including cluster systems like Slurm, LSF or PBS <https://docs.open-mpi.org/en/v5.0.x/launching-apps/quickstart.html>`_, the exact MPI wrapper command to launch programs may differ on your system.
To find out if and which command your execution environment provides, you will have to consult local documentation, check out if any known mpi wrapper commands are available or ask your system's administrators.
A good reference point for getting mpirun to work on your execution environment is the `documentation of the mpirun prerequisites <https://docs.open-mpi.org/en/v5.0.x/launching-apps/prerequisites.html>`_.

.. |set-resources| replace:: ``--set-resources``
.. _set-resources: https://snakemake.readthedocs.io/en/stable/executing/cli.html#snakemake.cli-get_argument_parser-execution

While a number of cluster scheduling systems are able to figure the ``tasks`` resource out for you, other execution environments will require setting ``-n`` manually.
This includes running a snakemake workflow with an MPI program `on a single host <https://docs.open-mpi.org/en/v5.0.x/launching-apps/quickstart.html#launching-on-a-single-host>`_ or `in a non-scheduled environment via ssh <https://docs.open-mpi.org/en/v5.0.x/launching-apps/quickstart.html#launching-in-a-non-scheduled-environments-via-ssh>`_, but will also include snakemake remote execution plugins for cluster systems that don't integrate handling this for you.
It is thus good practice to provide this explicitly.
To understand how a remote execution plugin for a particular cluster scheduling system supports MPI job execution, please consult the `documentation for the respective plugin <https://snakemake.github.io/snakemake-plugin-catalog/>`_.

In addition to overriding the MPI wrapper command, you can also provide extra parameters to the MPI wrapper command with the above construct, should your execution environment require it, for example:

.. code-block:: console

  $ snakemake --set-resources calc_pi:mpi="mpirun -arch x86" ...

or:

.. code-block:: console

  $ snakemake --set-resources calc_pi:mpi="srun --hint nomultithread" ...
