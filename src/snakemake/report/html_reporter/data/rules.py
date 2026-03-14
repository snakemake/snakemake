import json
from snakemake_interface_common.exceptions import WorkflowError


def render_rules(rules):
    return json.dumps(
        {
            rulename: {
                "input": rule.input,
                "output": rule.output,
                "conda_env": rule.conda_env,
                "container_img_url": rule.container_img_url,
                "code": render_code(rule),
                "n_jobs": rule.n_jobs,
            }
            for rulename, rule in rules.items()
        }
    )


def render_code(rule):
    try:
        from pygments.lexers import get_lexer_by_name
        from pygments.formatters import HtmlFormatter
        from pygments import highlight
        import pygments.util
    except ImportError:
        raise WorkflowError(
            "Python package pygments must be installed to create reports."
        )

    try:
        lexer = get_lexer_by_name(rule.language)

        highlighted = highlight(
            rule.source,
            lexer,
            HtmlFormatter(linenos=True, cssclass="source", wrapcode=True),
        )

        return highlighted
    except pygments.util.ClassNotFound:
        return [f'<pre class="source"><code>{rule.source}</code></pre>']
