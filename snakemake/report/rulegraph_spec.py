from functools import partial
from collections import defaultdict


def rulegraph_spec(dag):
    # Graphviz necessary to generate the layout
    try:
        import networkx as nx
        from networkx.drawing.nx_agraph import graphviz_layout
        from networkx.readwrite import json_graph
    except ImportError as e:
        raise WorkflowError(
            "Python packages networkx and pygraphviz must be "
            "installed to create reports.",
            e,
        )

    # Using sets make sure you have a single representative
    graph = defaultdict(set)
    for job in dag.jobs:
        graph[job.rule].update(dep.rule for dep in dag.dependencies[job])

    # Generate networkx directed graph, add nodes and edges
    g = nx.DiGraph()
    g.add_nodes_from([rule.name for rule in graph.keys()])
    # For each job, add one edge per dependency
    for rule, deps in graph.items():
        target = rule.name
        for dep in deps:
            source = dep.name
            g.add_edge(source, target)

    # Generate layout
    pos = graphviz_layout(g, prog="dot", args="-Grankdir=BT")
    xmax = max(x for x, y in pos.values()) + 100  # add offset to account for labels
    ymax = max(y for x, y in pos.values())

    # Encode nodes for VEGA template
    def encode_node(node):
        x, y = pos[node]
        return {"rule": node, "fx": x, "fy": y}

    nodes = list(map(encode_node, g.nodes))
    nodeids = {node["rule"]: i for i, node in enumerate(nodes)}

    # Encode links and links_direct for VEGA template
    def encode_link(edge):
        return {
            "target": nodeids[edge[0]],
            "source": nodeids[edge[1]],
            "targetrule": edge[0],
            "sourcerule": edge[1],
            "value": 1,
        }

    links = list(map(encode_link, g.edges))
    # I’m not sure why links_direct is necessary here, keeping it for compatibility
    links_direct = list(map(encode_link, g.edges))

    spec = (
        {
            "nodes": nodes,
            "links": links,
            "links_direct": links_direct,
        },
        xmax,
        ymax,
    )

    return spec
