import json
from snakemake_interface_common.exceptions import WorkflowError


def render_configfiles(configfiles):
    return json.dumps(
        [
            {
                "path": configfile.path,
                "code": render_code(configfile),
            }
            for configfile in configfiles.items()
        ]
    )


def render_code(configfile):
    try:
        from pygments.lexers import get_lexer_by_name
        from pygments.formatters import HtmlFormatter
        from pygments import highlight
    except ImportError:
        raise WorkflowError(
            "Python package pygments must be installed to create reports."
        )

    file_ext = configfile.path.suffix
    if file_ext in (".yml", ".yaml"):
        language = "yaml"
    elif file_ext == ".json":
        language = "json"
    else:
        raise ValueError(
            "Config file extension {} is not supported - must be YAML or JSON".format(
                file_ext
            )
        )

    source = configfile.source
    lexer = get_lexer_by_name(language)
    return highlight(
        source,
        lexer,
        HtmlFormatter(linenos=True, cssclass="source", wrapcode=True),
    )
