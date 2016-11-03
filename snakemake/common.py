DYNAMIC_FILL = "__snakemake_dynamic__"


class Mode:
    """
    Enum for execution mode of Snakemake.
    This handles the behavior of e.g. the logger.
    """
    default = 0
    subprocess = 1
    cluster = 2
