setlocal foldmethod=expr
setlocal foldexpr=GetSnakemakeFold(v:lnum)

function! GetSnakemakeFold(lnum)
    " fold preamble
    if a:lnum == 1
        return '>1'
    endif

    let thisline = getline(a:lnum)

    " blank lines end folds
    if thisline =~? '\v^\s*$'
        return '-1'
    " start fold on top level rules or python objects
    elseif thisline =~? '\v^(rule|def|checkpoint|class)'
        return ">1"
    elseif thisline =~? '\v^\S'
        if PreviousLineIndented(a:lnum) && NextRuleIndented(a:lnum)
            return ">1"
        endif
    endif

    return "="

endfunction

function! NextRuleIndented(lnum)
    let numlines = line('$')
    let current = a:lnum + 1

    while current <= numlines
        let thisline = getline(current)
        if thisline =~? '\v^(rule|def|checkpoint|class)'
            return 0
        elseif thisline =~? '\v^\s+(rule|checkpoint)'
            return 1
        endif

        let current += 1
    endwhile

    return 0
endfunction

function! PreviousLineIndented(lnum)
    let current = a:lnum - 1

    while current >= 1
        let thisline = getline(current)
        if thisline =~? '\v^\S'
            return 0
        elseif thisline =~? '\v^\s+\S'
            return 1
        endif

        let current -= 1
    endwhile

    return 0
endfunction
