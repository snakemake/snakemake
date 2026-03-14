from functools import partial


def rulegraph_spec(dag):
    # get toposorting, and keep only one job of each rule per level
    representatives = dict()

    toposorted = [
        get_representatives(level, representatives) for level in dag.toposorted()
    ]

    jobs = [job for level in toposorted for job in level]

    nodes = [
        {"rule": job.rule.name, "fx": 10, "fy": i * 50} for i, job in enumerate(jobs)
    ]
    idx = {job: i for i, job in enumerate(jobs)}

    xmax = 100
    ymax = max(node["fy"] for node in nodes)

    _get_links = partial(get_links, jobs, dag, idx, nodes, representatives)

    return (
        {
            "nodes": nodes,
            "links": list(_get_links(direct=False)),
            "links_direct": list(_get_links(direct=True)),
        },
        xmax,
        ymax,
    )


def get_representatives(level: list, representatives: dict):
    unique = dict()
    for job in level:
        if job.rule.name in unique:
            representatives[job] = unique[job.rule.name]
        else:
            representatives[job] = job
            unique[job.rule.name] = job
    return sorted(unique.values(), key=lambda job: job.rule.name)


def get_links(jobs, dag, idx, nodes, representatives, direct: bool):
    for u in jobs:
        for v in dag.dependencies[u]:
            target = idx[u]
            source = idx[representatives[v]]
            if target - source == 1:
                if not direct:
                    continue
            else:
                if direct:
                    continue

            yield {
                "target": target,
                "source": source,
                "sourcerule": nodes[source]["rule"],
                "targetrule": nodes[target]["rule"],
                "value": 1,
            }
