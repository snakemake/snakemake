" Vim ftdetect file
" Language: Snakemake (extended from python.vim)
" Maintainer: Jay Hesselberth (jay.hesselberth@gmail.com)
" Last Change: 2020 Oct 6
"
" Usage
"
" copy to $HOME/.vim/ftdetect directory
au BufNewFile,BufRead Snakefile set filetype=snakemake 
au BufNewFile,BufRead *.rules set filetype=snakemake 
au BufNewFile,BufRead *.snakefile set filetype=snakemake 
au BufNewFile,BufRead *.snake set filetype=snakemake 
