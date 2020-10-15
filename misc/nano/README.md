A nano syntax highlighting definition for Snakemake.

To use this file in nano, copy the `syntax/snakemake.nanorc` file
to your `$HOME` directory and add a line to your ~/.nanorc saying:

    include $HOME/snakemake.nanorc


NB. Line 12 of the syntax file contains a regular expression for
identifying a shebang (#!) header line. This command is not
supported in some versions of nano, so is commented out. If you
wish to enable it, simply remove the # at the beginning of the
line. If you are uncertain if your version of nano supports this
command, you may remove it and try. You will see an error if it
does not. Leaving the line commented out will result in the header
being treated as a regular comment for highlighting purposes. 
