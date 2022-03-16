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
    else:
        return " ".join(repr(v) for v in value)
