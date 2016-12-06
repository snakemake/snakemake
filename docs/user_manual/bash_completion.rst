.. _user_manual-bash_completion:

===============
Bash completion
===============

Snakemake supports bash completion for filenames, rulenames and arguments.
To enable it globally, just append

.. code-block:: bash

    `snakemake --bash-completion`

including the accents to your ``.bashrc``.
This only works if the ``snakemake`` command is in your path.
