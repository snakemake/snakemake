from snakemake.io import flag, is_flagged


def update(value):
    """
    A flag for an input or output file that shall be removed after usage.
    """
    return flag(value, "update")


def register_in_globals(_globals):
    _globals.update(
        {
            "update": update,
        }
    )