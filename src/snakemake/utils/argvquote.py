# compat: python 3.7 (script support)
# Must be kept compatible to Python 3.7 because it is used in Snakemake's
# Python script support.
# Only modules from python standard library or other Snakemake modules that are
# Python 3.7 compatible should be imported here (except for methods
# that are only called by Snakemake itself)!


def argvquote(arg, force=True):
    """Returns an argument quoted in such a way that CommandLineToArgvW
    on Windows will return the argument string unchanged.
    This is the same thing Popen does when supplied with a list of arguments.
    Arguments in a command line should be separated by spaces; this
    function does not add these spaces. This implementation follows the
    suggestions outlined here:
    https://blogs.msdn.microsoft.com/twistylittlepassagesallalike/2011/04/23/everyone-quotes-command-line-arguments-the-wrong-way/
    """
    if not force and len(arg) != 0 and not any([c in arg for c in ' \t\n\v"']):
        return arg
    else:
        n_backslashes = 0
        cmdline = '"'
        for c in arg:
            if c == "\\":
                # first count the number of current backslashes
                n_backslashes += 1
                continue
            if c == '"':
                # Escape all backslashes and the following double quotation mark
                cmdline += (n_backslashes * 2 + 1) * "\\"
            else:
                # backslashes are not special here
                cmdline += n_backslashes * "\\"
            n_backslashes = 0
            cmdline += c
        # Escape all backslashes, but let the terminating
        # double quotation mark we add below be interpreted
        # as a metacharacter
        cmdline += +n_backslashes * 2 * "\\" + '"'
        return cmdline


def cmd_exe_quote(arg):
    """Quotes an argument in a cmd.exe compliant way."""
    arg = argvquote(arg)
    cmd_exe_metachars = '^()%!"<>&|'
    for char in cmd_exe_metachars:
        arg.replace(char, "^" + char)
    return arg
