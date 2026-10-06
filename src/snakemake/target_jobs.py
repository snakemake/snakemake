from snakemake_interface_executor_plugins.utils import TargetSpec

from snakemake.common.misc import parse_key_value_arg


def parse_target_jobs_cli_args(target_jobs_args):
    errmsg = "Invalid target wildcards definition: entries have to be defined as WILDCARD=VALUE pairs"
    if target_jobs_args is not None:
        target_jobs = list()
        for entry in target_jobs_args:
            rulename, wildcards = entry.split(":", 1)
            if wildcards:

                def parse_wildcard(entry):
                    key, value = parse_key_value_arg(entry, errmsg, strip_quotes=False)
                    # The encoder only wraps a value in double quotes when the
                    # value contains a comma, so a matched pair of surrounding
                    # quotes around a comma-containing value was added by
                    # encoding and must be removed. Quotes are otherwise part
                    # of the wildcard value itself.
                    if (
                        len(value) >= 2
                        and value[0] == value[-1] == '"'
                        and "," in value
                    ):
                        value = value[1:-1]
                    return key, value

                wildcards = dict(
                    parse_wildcard(entry)
                    for entry in _split_at_unquoted_commas(wildcards)
                )
                target_jobs.append(TargetSpec(rulename, wildcards))
            else:
                target_jobs.append(TargetSpec(rulename, dict()))
        return target_jobs


def _split_at_unquoted_commas(arg):
    items = []
    buf = []
    quoted = False
    prev = ""
    equals_in_item = 0
    for pos, char in enumerate(arg):
        # The encoder only wraps values containing commas in double quotes, so
        # a quote is structural only when it opens a value (right after the
        # entry's first "=") or closes one (followed by "," or the end of the
        # argument). Quotes elsewhere are literal.
        if char == '"':
            if not quoted:
                if prev == "=" and equals_in_item == 1:
                    quoted = True
            elif pos + 1 == len(arg) or arg[pos + 1] == ",":
                quoted = False
        elif char == "=":
            equals_in_item += 1
        elif char == "," and not quoted:
            items.append("".join(buf))
            buf = []
            equals_in_item = 0
            prev = char
            continue
        buf.append(char)
        prev = char
    items.append("".join(buf))
    return items
