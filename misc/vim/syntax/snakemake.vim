" Vim syntax file
" Language: Snakemake (extended from python.vim)
" Maintainer: Jay Hesselberth (jay.hesselberth@gmail.com)
" Last Change: 2025 Mar 6
"
" Usage
"
" copy to $HOME/.vim/syntax directory and
" copy to ftdetect/snakemake.vim to $HOME/.vim/ftdetect directory
"
" force coloring in a vim session with:
"
" :set syntax=snakemake
"
if exists("b:current_syntax")
    finish
endif

" load settings from system python.vim (7.4)
source $VIMRUNTIME/syntax/python.vim
source $VIMRUNTIME/indent/python.vim

" Snakemake rules, as of version 8.29
"
"
" rule         = "rule" (identifier | "") ":" ruleparams
" include      = "include:" stringliteral
" workdir      = "workdir:" stringliteral
" module       = "module" identifier ":" moduleparams
" configfile   = "configfile" ":" stringliteral
" userule      = "use" "rule" (identifier | "*") "from" identifier ["as" identifier] ["with" ":" norunparams]
" ni           = NEWLINE INDENT
" norunparams  = [ni input] [ni output] [ni params] [ni message] [ni threads] [ni resources] [ni log] [ni conda] [ni container] [ni benchmark] [ni cache] [ni priority]
" ruleparams   = norunparams [ni (run | shell | script | notebook)] NEWLINE snakemake
" input        = "input" ":" parameter_list
" output       = "output" ":" parameter_list
" params       = "params" ":" parameter_list
" log          = "log" ":" parameter_list
" benchmark    = "benchmark" ":" statement
" cache        = "cache" ":" bool
" message      = "message" ":" stringliteral
" threads      = "threads" ":" integer
" priority     = "priority" ":" integer
" resources    = "resources" ":" parameter_list
" version      = "version" ":" statement
" conda        = "conda" ":" stringliteral
" container    = "container" ":" stringliteral
" run          = "run" ":" ni statement
" shell        = "shell" ":" stringliteral
" script       = "script" ":" stringliteral
" notebook     = "notebook" ":" stringliteral
" moduleparams = [ni snakefile] [ni metawrapper] [ni config] [ni skipval]
" snakefile    = "snakefile" ":" stringliteral
" metawrapper  = "meta_wrapper" ":" stringliteral
" config       = "config" ":" stringliteral
" skipval      = "skip_validation" ":" stringliteral


" general directives (e.g. input)
syn keyword pythonStatement 
      \ benchmark
      \ conda
      \ configfile
      \ container
      \ default_target
      \ envmodules
      \ group
      \ include
      \ input
      \ localrule
      \ localrules
      \ log
      \ message
      \ notebook
      \ onerror
      \ onstart
      \ onsuccess
      \ output
      \ params
      \ priority
      \ resources
      \ ruleorder
      \ run
      \ scattergather
      \ script
      \ shadow
      \ shell
      \ singularity
      \ snakefile
      \ template_engine
      \ threads
      \ version
      \ wildcard_constraints
      \ wildcards
      \ workdir
      \ wrapper

" directives with a label (e.g. rule)
syn keyword pythonStatement 
      \ checkpoint
      \ rule
      \ subworkflow
      \ use
      \ nextgroup=pythonFunction skipwhite

" common snakemake objects
syn keyword pythonBuiltinObj 
      \ Paramspace
      \ checkpoints
      \ config
      \ gather
      \ rules
      \ scatter
      \ workflow

" snakemake functions
syn keyword pythonBuiltinFunc 
      \ ancient
      \ before_update
      \ directory
      \ expand
      \ exists
      \ multiext
      \ pipe
      \ protected
      \ read_job_properties
      \ service
      \ temp
      \ touch
      \ unpack
      \ update

" similar to special def and class treatment from python.vim, except
" parenthetical part of def and class
syn match pythonFunction
      \ "\%(\%(rule\s\|subworkflow\s\|checkpoint\s\)\s*\)\@<=\h\w*" contained

syn sync match pythonSync grouphere NONE "^\s*\%(rule\|subworkflow\|checkpoint\)\s\+\h\w*\s*"

let b:current_syntax = "snakemake"

" vim:set sw=2 sts=2 ts=8 noet:
