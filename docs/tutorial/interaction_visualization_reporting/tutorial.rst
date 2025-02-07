.. _interaction_visualization_reporting_tutorial:

==================================================================
Tutorial: Interaction, Visualization, and Reporting with Snakemake
==================================================================

.. _Snakemake: https://snakemake.github.io
.. _Jupyter: https://jupyter.org
.. _Datavzrd: https://datavzrd.github.io
.. _Altair: https://altair-viz.github.io

Via its Jupyter_ :ref:`notebook <snakefiles_notebook-integration>` and :ref:`script <snakefiles-external_scripts>` integration, Snakemake_ offers powerful support for interactive data exploration.
This is particularly useful for the last mile of data analysis, where results obtained from established tools and libraries are summarized and visualized for e.g. presentation in a scientific publication.
Moreover, via its :ref:`reporting <snakefiles-reports>` capabilities, Snakemake can be used to generate reproducible reports that contain the results of a data analysis pipeline, the code that was used to generate them, and the environment in which the pipeline was executed.
The latter integrates well with Datavzrd_, a tool for generating interactive views of tabular data.

In this short tutorial, we bring all these capabilities together in order to demonstrate how Snakemake enables last mile data analysis without losing reproducibility, adaptability, and transparency.
The tutorial assumes that you already know how to write Snakemake workflows in general, at least by having conducted the entire :ref:`general tutorial <tutorial>`.
For simplicity, this tutorial does not use any :ref:`wildcards <snakefiles-wildcards>`, which render both the notebook integration and Snakemake's last mile capabilities even more powerful.

Setup
-----

Install Snakemake into a new Conda environment as instructed in the :ref:`installation guide <getting_started-installation>`.
In the following, we assume that you conduct all commands within this environment.
If you already have Snakemake installed, make sure that you use the latest stable version.

Next, create an empty working directory at a reasonable place in your file system.
In the following, we assume that you have an open terminal inside of that directory.
If you use visual studio code (recommended), open a new instance inside of that directory and open a terminal via the terminal menu.

Step 1: Obtain example data
---------------------------

First, we need some example data.
We will conduct this step as our first interactive exploration via Snakemake's Jupyter integration.
We create a new directory ``workflow`` and inside of that an empty ``Snakefile``.
In the ``Snakefile``, define the following rule:

.. code-block:: python

    
    rule get_data:
        output:
            "resources/data/cars.tsv"
        conda:
            "envs/download.yaml"
        log:
            notebook="logs/get_data.ipynb"
        notebook:
            "notebooks/get_data.py.ipynb"

The rule shall create an output TSV file ``resources/data/cars.tsv`` that will contain our example dataset.
It defines a Conda environment with the required software in ``envs/download.yaml``.
Finally, it specifies that the output shall be created via the Jupyter_ notebook ``notebooks/get_data.py.ipynb``.
Snakemake's notebook support allows to write a notebook with the per-cell output as a log file, which we utilize here by specifying ``notebook="logs/get_data.ipynb"`` in the ``log`` directive.
After saving the new rule into the ``Snakefile``, we store the following in the Conda environment file ``envs/download.yaml``:

.. code-block:: yaml

    channels:
      - conda-forge
      - nodefaults
    dependencies:
      - python =3.11
      - polars =1.1
      - vega_datasets =0.9
      - ipykernel =6.29
      - notebook =7.2

Since we use a notebook for the step, we not only need ``python``, ``polars`` and ``vega_datasets`` for data processing, but also ``ipykernel`` (for Jupyter's Python support) and ``notebook`` for being able to run a Jupyter notebook.

Now, instead of creating the notebook ourselves, we let Snakemake do the job by executing the following in a terminal:

.. code-block:: console

    $ snakemake --sdm conda --cores 1 --edit-notebook resources/data/cars.tsv

This tells Snakemake to create a skeleton notebook and start a Jupyter server.
The output of the server provides three options to open the notebook in the browser.
Use Ctrl-Click on one of the options to open the notebook server in your browser.
In the presented interface, select the temporary notebook file and start editing.
We aim for the following content:

.. code-block:: python
    
    import polars as pl
    from vega_datasets import data

    cars = pl.from_pandas(data.cars()).with_columns(
        pl.col("Year").dt.year()
    ).select(
        pl.col("*").name.map(lambda name: name.lower().replace("_", " "))
    )

    cars.write_csv(snakemake.output[0], separator="\t")

This code snippet loads the cars dataset from the ``vega_datasets`` package, converts the ``Year`` column (which actually contains dates) to an integer representing just the year, normalizes column names, and writes the resulting table to the output file.
After saving the notebook, stop the Jupyter server by selecting shut down in the ``File`` menu.
Snakemake then cleans up your notebook and stores it in the desired place.

Future reruns of this rule can just treat the notebook as an ordinary script.
Whenever you want to modify the notebook, you can do so in interactive mode again using the ``--edit-notebook`` option in combination with the output file path.

Step 2: Create a plot with R
----------------------------

We now add a rule that creates a plot from the given dataset using the `R tidyverse <https://www.tidyverse.org>`__.
Add the following rule to the ``Snakefile``:

.. code-block:: python

    
    rule plot_with_r:
        input:
            "resources/data/cars.tsv"
        output:
            "results/plots/horsepower_vs_mpg.ggplot.svg",
        log:
            notebook="logs/plot_horsepower_vs_mpg.r.ipynb"
        conda:
            "envs/rstats.yaml"
        notebook:
            "notebooks/plot_horsepower_vs_mpg.r.ipynb"

Analogously to before, we specify a Conda environment in ``envs/rstats.yaml`` with the following content:

.. code-block:: yaml

    channels:
      - conda-forge
      - nodefaults
    dependencies:
      - r-base =4.4
      - r-readr =2.1.5
      - r-dplyr =1.1
      - r-ggplot2 =3.5
      - notebook =7.2
      - r-irkernel =1.3

Here, we require ``r-irkernel`` instead of the ``ipykernel`` for the Python case.
Now recall the terminal command from the previous step and run the equivalent for the new rule:

.. code-block:: console

    $ snakemake --sdm conda --cores 1 --edit-notebook results/plots/horsepower_vs_mpg.ggplot.svg

Again open the notebook in your browser and start editing.
We want to create a plot of the horsepower vs. miles per gallon from the cars dataset.
The plot shall finally be saves as an SVG file:

.. code-block:: r
    
    library(readr)
    library(ggplot2)

    cars <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)
    svg(snakemake@output[[1]])
    ggplot(cars, aes(`miles per gallon`, horsepower)) + geom_point() + theme_classic(16)
    dev.off()

After saving the notebook, shut down the Jupyter server as before (the ``File`` menu) and let Snakemake clean up and store the notebook automatically.

Step 3: Create a plot with Python
---------------------------------

Now, for illustration purposes, we want to create the same plot with Python, using the plotting library Altair_.
Add the following rule to the ``Snakefile``:

.. code-block:: python

    
    rule plot_with_python:
        input:
            "resources/data/cars.tsv"
        output:
            "results/plots/horsepower_vs_mpg.altair.html",
        log:
            notebook="logs/plot_horsepower_vs_mpg.py.ipynb"
        conda:
            "envs/pystats.yaml"
        notebook:
            "notebooks/plot_horsepower_vs_mpg.py.ipynb"

The corresponding Conda environment in ``envs/pystats.yaml`` is:

.. code-block:: yaml

    channels:
      - conda-forge
      - nodefaults
    dependencies:
      - python =3.11
      - polars =1.1
      - altair =5.3
      - altair_saver =0.5
      - vl-convert-python =1.5
      - vegafusion =1.6
      - vegafusion-python-embed =1.6
      - notebook =7.2
      - ipykernel =6.29

In addition to the packages already used for the download step, we now require ``altair``, ``altair_saver`` and ``vl-convert-python`` for output format support in Altair_. In addition, adding the two ``vegafusion`` packages, enables support for efficient plotting that involves a lot of datapoints.
While we don't need that in this example, it is a good practice to include them in the environment file in order to be prepared for such cases.

Adapt the edit notebook command from above to edit this notebook interactively.
The content of the notebook shall be:

.. code-block:: python
    
    import altair as alt
    import polars as pl
    alt.data_transformers.enable("vegafusion")

    data = pl.read_csv(snakemake.input[0], separator="\t")

    chart = alt.Chart(data).mark_point(tooltip=True).encode(
        alt.X("miles per gallon"),
        alt.Y("horsepower"),
        alt.Color("origin").scale(scheme="accent"),
    ).interactive()

    chart.save(snakemake.output[0])

Here, in addition to the plot before, we color the points by the origin of the car.
Moreover, we define the chart to be interactive and offer tooltips at each point.
After running and saving the notebook, followed by shutting down the Jupyter server like before, open the generated output file in your browser and explore the interactivity (zoom with the mouse wheel, pan by click-hold-move, and hover the points for tooltips).

Step 4: Create an interactive table view with Datavzrd
------------------------------------------------------

While plots are a central part of data exploration and enable to reveal e.g. relationships between variables, providing transparent access to the underlying dataset is crucial for generating trust and easing communication of the results.
Datavzrd_ is a tool that allows to rapidly generate interactive views of tabular data without requiring a web server to be set up and maintained.
It is available as a `Snakemake wrapper <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/datavzrd.html>`__ and can be included here as follows:

.. code-block:: python

    
    rule view_with_datavzrd:
        input:
            config=workflow.source_path("resources/datavzrd/cars.yaml"),
            table="resources/data/cars.tsv",
        output:
            report(
                directory("results/tables/cars"),
                htmlindex="index.html",
                caption="report/cars.rst",
                category="Tables",
                labels={"table": "cars"},
            ),
        log:
            "logs/datavzrd.log",
        wrapper:
            "v4.7.2/utils/datavzrd"

Note that the wrapper already suggests to annotate the output for inclusion into the Snakemake report.
While we ignore that for now, we will revisit it at a later step.
Snakemake wrappers are a powerful way to integrate external tools into a Snakemake workflow.
They come with predefined Conda environments (so that we don't need to specify any software here) and usually wrap around a library or command line tool using e.g. Python or R scripts.
The Datavzrd wrapper expects one or multiple tables as input, as well as a configuration file that defines the appearance and behavior of the table view(s).
In this case, this config file shall be stored in ``resources/datavzrd/cars.yaml`` and contain the following:

.. code-block:: yaml

    __use_yte__: true

    datasets:
      cars:
        path: ?input.table
        separator: "\t"

    views:
      cars:
        dataset: cars
        render-table:
          columns:
            name:
              link-to-url:
                Wikipedia:
                url: "https://en.wikipedia.org/wiki/{value}"
            miles per gallon:
              plot:
                ticks:
                  scale: linear
            cylinders:
              plot:
                heatmap:
                  scale: linear
                  domain:
                    - 0
                    - 16
                  range:
                    # white to blue
                    - "#ffffff"
                    - "#6baed6"
            displacement:
              plot:
                ticks:
                  scale: linear
              display-mode: detail
            horsepower:
              plot:
                ticks:
                  scale: linear
              display-mode: detail
            weight in lbs:
              plot:
                ticks:
                  scale: linear
            acceleration:
              plot:
                ticks:
                  scale: linear
            year:
              plot:
                ticks:
                  scale: linear
            origin:
              plot:
                heatmap:
                  scale: ordinal
                  color-scheme: category10

While we refer to the Datavzrd_ documentation for the details, let us highlight a few important aspects here:
The configuration file defines a dataset named ``cars`` that is loaded from the input table file.
Since the path of the table file shall rather not be hardcoded into the configuration we use the `YTE template engine <https://yte-template-engine.github.io>`__ to render the path dynamically into the config file, thereby accessing ``input.table`` as provided by the Snakemake job.
In this case, we define one view on the dataset (in principle there can be multiple views on one or multiple tables encoded into one Datavzrd report).
For each column of the table, we specify how it shall be visualized (as a link, as tick or heatmap plot).
Beyond these, Datavzrd offers a lot more possibilities, including the ability to publish `spells <https://datavzrd.github.io/docs/spells.html>`__ with custom visualizations for common column types.

Step 5: Adding default targets
------------------------------

While we have so far generated each output file manually via the ``--edit-notebook`` option, it is time to define the default targets of our workflow.
Add a rule ``all`` to the top of the ``Snakefile``:

.. code-block:: python

    
    rule all:
        input:
            "results/plots/horsepower_vs_mpg.ggplot.svg",
            "results/plots/horsepower_vs_mpg.altair.html",
            "results/tables/cars",

This rule just defines input files.
Since Snakemake wants to run the first rule in the workflow by default, all the inputs have to be created by combinations of other rules, such that the desired files are generated.

Run the workflow with

.. code-block:: console

    $ snakemake --sdm conda --cores 1

Since the plots are already present, Snakemake will just run the Datavzrd rule.
Afterwards explore the Datavzrd output by opening ``results/tables/cars/index.html`` in your browser.

Step 6: Reporting
-----------------

At the end of the last mile, the data analysis results are usually communicated, e.g. in a scienficic manuscript or by first sending them to collaborators.
Just communicating the individual outputs disconnects them from the code that generated them, which is a problem for transparency.
To overcome this caveat, Snakemake offers the generation of automatic reports, which connect code and parameters with results.
Since these reports are standalone, server free HTML files, they can be easily added as supplementary files to a manuscript or shared with collaborators.
The Datavrzrd rule already marks the output for inclusion in the report.
It does so by specifying a category, a caption, and labels for showing the output inside of the report menu structure.
Let us now add equivalent annotations to the two plot rules by editing them into the following form:

.. code-block:: python

    rule plot_with_r:
        input:
            "resources/data/cars.tsv"
        output:
            report(
                "results/plots/horsepower_vs_mpg.ggplot.svg",
                category="Plots",
                labels={"plot": "horsepower_vs_mpg", "approach": "ggplot"},
                caption="report/horsepower_vs_mpg.rst",
            ),
        log:
            notebook="logs/plot_horsepower_vs_mpg.r.ipynb"
        conda:
            "envs/rstats.yaml"
        notebook:
            "notebooks/plot_horsepower_vs_mpg.r.ipynb"


    rule plot_with_python:
        input:
            "resources/data/cars.tsv"
        output:
            report(
                "results/plots/horsepower_vs_mpg.altair.html",
                category="Plots",
                labels={"plot": "horsepower_vs_mpg", "approach": "altair"},
                caption="report/horsepower_vs_mpg.rst",
            ),
        log:
            notebook="logs/plot_horsepower_vs_mpg.py.ipynb"
        conda:
            "envs/pystats.yaml"
        notebook:
            "notebooks/plot_horsepower_vs_mpg.py.ipynb"

The caption files shall contain a human readable description, which we will keep short in this example case.
Create the file ``workflow/report/cars.rst`` with the following content:

.. code-block:: rst

    The cars dataset as provided by the vega project: https://github.com/vega/vega-datasets.

Create the file ``workflow/report/horsepower_vs_mpg.r.rst`` with the following content:

.. code-block:: rst

    A plot of the horsepower vs. miles per gallon from the cars dataset, created with ggplot2.

Create the file ``workflow/report/horsepower_vs_mpg.py.rst`` with the following content:

.. code-block:: rst

    A plot of the horsepower vs. miles per gallon from the cars dataset, created with Altair.

Finally, add a global report directive to the top of the ``Snakefile``

.. code-block:: python

    report: "workflow/report/workflow.rst"

and create the file ``workflow/report/workflow.rst`` with the following content:

.. code-block:: rst

    Workflow illustrating Snakemake's capabilities for last mile data analysis with Jupyter, Datavzrd, and R/Python plotting libraries.

Now create the report by running

.. code-block:: console

    $ snakemake --sdm conda --report report.zip

This will create a report in the form of a zip file.
Unzip it with

.. code-block:: console

    $ unzip report.zip

and open the file ``report/report.html`` file in your browser.
You will see that the report brings all desired outputs together in a structured way, including the captions and the global description.
It not only allows to view the results, but also to explore the code of each rule and all involved parameters and software tools.
This way, it generates transparency without requiring people to manually inspect the workflow codebase.
In many ways, it can be seen as a self-contained next generation supplementary file of a scientific manuscript.