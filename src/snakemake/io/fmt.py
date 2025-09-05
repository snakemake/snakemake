from snakemake.common.tbdstring import TBDString
from snakemake.io import get_flag_value, is_flagged


def fmt_iofile(f, as_input: bool = False, as_output: bool = False):
    as_io = as_input or as_output

    if f.is_storage:
        if as_input:
            if f.storage_object.retrieve:
                storage_phrase = "retrieve from storage"
            else:
                storage_phrase = "only provide storage uri"
        elif as_output:
            storage_phrase = "send to storage"
        else:
            storage_phrase = ""
        f_str = f.storage_object.print_query
    else:
        f_str = f
        storage_phrase = ""

    def annotate(f_str, label=""):
        sep = ", " if label and storage_phrase else ""
        ann = f" ({label}{sep}{storage_phrase})" if label or storage_phrase else ""
        return f"{f_str}{ann}"

    if is_flagged(f, "sourcecache_entry"):
        orig_path_or_uri = get_flag_value(f, "sourcecache_entry")
        return annotate(orig_path_or_uri, "cached")
    if as_io:
        if is_flagged(f, "pipe"):
            return annotate(f_str, "pipe")
        elif is_flagged(f, "service"):
            return annotate(f_str, "service")
        elif is_flagged(f, "nodelocal"):
            return annotate(f_str, "nodelocal")
        elif is_flagged(f, "update"):
            return annotate(f_str, "update")
        elif is_flagged(f, "before_update"):
            return annotate(f_str, "before update")
        elif is_flagged(f, "access_pattern"):
            pattern = f"access: {get_flag_value(f, 'access_pattern')}"
            return annotate(f_str, pattern)
        elif is_flagged(f, "checkpoint_target"):
            return TBDString()
        else:
            return annotate(f_str)
    return f
