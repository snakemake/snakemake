from snakemake.io import not_iterable


def format_cli_arg(flag, value, quote=True, skip=False):
    if not skip and value:
        if isinstance(value, bool):
            value = ""
        else:
            value = format_cli_pos_arg(value, quote=quote)
        return f"{flag} {value}"
    return ""


def format_cli_pos_arg(value, quote=True):
    if not_iterable(value):
        return repr(value)
    elif isinstance(value, dict):
        return join_cli_args(repr(f"{key}={val}") for key, val in value)
    else:
        return join_cli_args(repr(v) for v in value)


def join_cli_args(args):
    return " ".join(arg for arg in args if arg)
