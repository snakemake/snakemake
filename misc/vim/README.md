A vim syntax highlighting definition for Snakemake.
You can copy the `snakemake.vim` file to `$HOME/.vim/syntax` directory and add

    au BufNewFile,BufRead Snakefile set syntax=snakemake
    au BufNewFile,BufRead *.rules set syntax=snakemake
    au BufNewFile,BufRead *.snakefile set syntax=snakemake
    au BufNewFile,BufRead *.snake set syntax=snakemake

to your `$HOME/.vimrc` file. Highlighting can be forced in a vim session with `:set syntax=snakemake`.

To install via Vundle use:

    Plugin 'https://github.com/snakemake/snakemake.git', {'rtp': 'misc/vim/'}

To install via [vim-plug]( https://github.com/junegunn/vim-plug):

    Plug 'snakemake/snakemake', {'rtp': 'misc/vim'}

