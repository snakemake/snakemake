__author__ = "Thomas Mulvaney"
__copyright__ = "Copyright 2026, Thomas Mulvaney"
__email__ = "mulvaney@mailbox.org"
__license__ = "MIT"
from collections import defaultdict
import textwrap
import html


def resolve_functions(input_files):
    """Iterate over all input files and replace input functions with a
    fixed string.
    """
    files = []
    for f in input_files:
        if callable(f):
            files.append("<input function>")
        else:
            files.append(repr(f).strip("'"))
    return files


def enumerate_input_files(job):
    """Given a job, enumerates all outputs as tuples in the form:

        (id, job_file, rule_file)

    where job_file is a real path and rule_file is the input string as defined in the rule
    """
    for i, (job_file, rule_file) in enumerate(
        zip(job._input, resolve_functions(job.rule._input))
    ):
        yield (i, job_file, rule_file)


def enumerate_output_files(job):
    """Given a job, enumerates all inputs as tuples:

         (id, job_file, rule_file)

    where job_file is a real path and rule_file is the output as defined in the rule
    """
    for i, (job_file, rule_file) in enumerate(zip(job._output, job.rule._output)):
        yield (i, job_file, rule_file)


def rule_file_edges(dag, node_to_id):
    """Assuming nodes are rules with input and output files, produces
    edges connecting input and output files."""
    edges = defaultdict(int)

    for job in dag.jobs:
        for i, job_input, rule_input in enumerate_input_files(job):
            for dep in dag._dependencies[job]:
                for j, dep_output, dep_rule_output in enumerate_output_files(dep):
                    if job_input == dep_output:
                        edges[
                            (
                                f"{node_to_id[dep.rule]}:out{j}",
                                f"{node_to_id[job.rule]}:in{i}",
                            )
                        ] += 1
    return edges


def job_file_edges(dag, node_to_id):
    """Assuming nodes are jobs with input and output files, produces
    edges connecting input and output files."""
    edges = defaultdict(int)

    for job in dag.jobs:
        for i, job_input, _ in enumerate_input_files(job):
            for dep in dag._dependencies[job]:
                for j, dep_output, _ in enumerate_output_files(dep):
                    if job_input == dep_output:
                        edges[
                            (f"{node_to_id[dep]}:out{j}", f"{node_to_id[job]}:in{i}")
                        ] += 1
    return edges


def job_edges(dag, node_to_id) -> dict[(str, str), int]:
    """Assuming nodes represent jobs, produces edges connecting these jobs"""
    edges = defaultdict(int)
    for job in dag.jobs:
        for dep in dag._dependencies[job]:
            edges[(f"{node_to_id[dep]}", f"{node_to_id[job]}")] += 1
    return edges


def rule_edges(dag, node_to_id) -> dict[(str, str), int]:
    """Assuming nodes represent rules, produces edges connecting these rules"""
    edges = defaultdict(int)
    for job in dag.jobs:
        for dep in dag._dependencies[job]:
            edges[(f"{node_to_id[dep.rule]}", f"{node_to_id[job.rule]}")] += 1
    return edges


def rule_to_nodes(dag):
    """Given a DAG create a dictionary of nodes.

    The dictionary maps the rule to an a numerical ID used by graphviz.
    """
    node_cnt = 0
    nodes = {}
    for job in dag.jobs:
        if nodes.get(job.rule) is None:
            nodes[job.rule] = node_cnt
            node_cnt += 1
    return nodes


def job_to_nodes(dag):
    """Given a DAG create a dictionary of nodes.

    The dictionary maps the job node to an a numerical ID used by graphviz.
    """
    node_cnt = 0
    nodes = {}
    for job in dag.jobs:
        if nodes.get(job) is None:
            nodes[job] = node_cnt
            node_cnt += 1
    return nodes


def rule_to_color(dag):
    # Sorting the rules by name before assigning colors gives deterministic output
    rules = sorted(set(job.rule for job in dag.jobs), key=lambda x: x.name)

    def hsv_to_htmlhexrgb(h, s, v):
        """Convert hsv colors to hex-encoded rgb colors usable by html."""
        import colorsys

        hex_r, hex_g, hex_b = (round(255 * x) for x in colorsys.hsv_to_rgb(h, s, v))
        return "#{hex_r:0>2X}{hex_g:0>2X}{hex_b:0>2X}".format(
            hex_r=hex_r, hex_g=hex_g, hex_b=hex_b
        )

    huefactor = 2 / (3 * len(rules))
    return {
        rule: hsv_to_htmlhexrgb(i * huefactor, 0.6, 0.85)
        for i, rule in enumerate(rules)
    }


def job_rule_to_color(dag):
    """A reasonable default for coloring jobs, is to color them according to rule"""
    r2c = rule_to_color(dag)
    return {j: r2c[j.rule] for j in dag.jobs}


def render_graph(nodes, edges, edge_labels=False, rankdir="LR"):
    edges = render_edges(edges, edge_labels)
    return textwrap.dedent("""\
        digraph snakemake_dag {{
            graph[bgcolor=white, margin=0, rankdir={rankdir}, ranksep=0.5];
            node[shape=box, style=rounded, fontname=sans, \
            fontsize=10, penwidth=2];
            edge[penwidth=2, color=grey];
        {items}
        }}\
        """).format(rankdir=rankdir, items="\n".join(nodes + list(edges)))

def render_edges(edges, edge_labels=False):
    if edge_labels:
        return [f"{source} -> {dest} [ label={edges[edge]} ]" for source, dest in edges]
    else:
        return [f"{source} -> {dest}" for source, dest in edges]


def render_node(
    node, node_to_id, node_to_color, dag, wildcards=False, inputs=False, outputs=False
) -> str:
    """Render the contents of a node. Takes as input the `node` which
    is either a Snakemake Job or Rule, a dictionary mapping nodes to
    integer ids `node_to_id` and a dictionary mapping nodes to colors
    `node_to_color`.

    The node can be rendered with additional fields, by setting the following
    boolean flags to True.

    * wildcards - show new wildcards introduced by the rule.
    * inputs    - show file inputs of the job or rule
    * outputs   - show file outputs of the job or rule

    """
    input_files = resolve_functions(node._input)
    input_file_ids = {input_file: id for id, input_file in enumerate(input_files)}
    output_files = resolve_functions(node._output)
    output_file_ids = {output_file: id for id, output_file in enumerate(output_files)}
    color = node_to_color[node]

    def render_file_row(value, port):
        return f'<tr><td port="{port}" align="left"><font face="monospace">{value}</font></td></tr>'

    def render_param_row(value):
        return f'<tr><td align="left"><font face="monospace">{value}</font></td></tr>'

    def render_header(value):
        return f'<tr><td align="left"><b>{value}</b></td></tr>'

    def render_wildcard_rows(node):
        # Note this is only compatible with nodes that are jobs, and
        # throws an exception on rules.
        html = ""
        try:
            new_wildcards = dag.new_wildcards(node)
            if len(new_wildcards) > 0:
                html += render_header("wildcards")
                for name, value in new_wildcards:
                    html += f'<tr><td align="left"><font face="monospace"><b>{name}:</b>{value}</font></td></tr>'
        except Exception:
            pass
        return html

    color = node_to_color[node]

    html_node = (
        f'<table border="2" color="{color}" cellspacing="0" cellborder="0">'
        f'<tr><td bgcolor="{color}">'
        f'<b><font point-size="10">{node.name}</font></b>'
        "</td></tr>"
    )

    if inputs:
        html_node += render_header("input")
        for f in sorted(input_files):
            html_node += render_file_row(html.escape(f), f"in{input_file_ids[f]}")

    if wildcards:
        html_node += render_wildcard_rows(node)

    if outputs:
        html_node += render_header("output")
        for f in sorted(output_files):
            html_node += render_file_row(html.escape(f), f"out{output_file_ids[f]}")

    html_node += "</table>"
    return f"{node_to_id[node]} [ shape=none, margin=0, label=<{html_node}>]"


def rule_graph(dag, edge_labels=True, file_edges=False, inputs=False, outputs=False):
    """Print a rule graph.

    Options to be turned into CLI args:
        * edge_labels, --edge-labels, Show the number of inputs a given edge represents
        * file_edges, --file-edges,   Show edges between files in the graph (Shows inputs and output files).
        * inputs,     --show-inputs   Show input files
        * outputs,    --show-oututs   Show output files
    """
    nodes = rule_to_nodes(dag)
    node_to_color = rule_to_color(dag)
    if file_edges:
        inputs = True
        outputs = True
        edges = rule_file_edges(dag, nodes)
    else:
        edges = rule_edges(dag, nodes)
    render_nodes = [
        render_node(
            node,
            nodes,
            node_to_color,
            dag,
            wildcards=False,
            inputs=inputs,
            outputs=outputs,
        )
        for node in nodes
    ]
    print(render_graph(render_nodes, edges, edge_labels))


def job_graph(dag, inputs=False, outputs=False,file_edges=False, wildcards=False):
    """Print a job graph.

    Options to be turned into CLI args:
        * file_edges, --file-edges,   Show edges between files in the graph (shows input and output files).
        * inputs,     --show-inputs   Show input files
        * outputs,    --show-oututs   Show output files
        * wildcards,  --show-wildcards   Show wildcards.
    """

    nodes = job_to_nodes(dag)
    node_to_color = job_rule_to_color(dag)
    if file_edges:
        inputs = True
        outputs = True
        edges = job_file_edges(dag, nodes)
    else:
        edges = job_edges(dag, nodes)
    render_nodes = [
        render_node(
            node, nodes, node_to_color, dag, wildcards=False, inputs=inputs, outputs=outputs
        )
        for node in nodes
    ]
    print(render_graph(render_nodes, edges,edge_labels=False))
