def compare_version_geq(m:str,n:str) -> bool:
    """
    a poor man implementation of version comparison. 
    External implementations are not used, since Snakemake inside containers
     external packages do not work. e.g.: "tests/tests.py::test_singularity" 
     would fail

    Only tested with version number of conda,singularity and apptainer
    """
    if m==n:
        return True
    for q,r in zip(m.split("."),n.split(".")):
        
        if q.isdigit():
            q=int(q)
        elif "-rc" in q:
            q=int(q.split("-rc")[0])-1

        if r.isdigit():
            r=int(r)
        elif "-rc" in r:
            r=int(r.split("-rc")[0])-1

        if q>r:
            return True
        if q<r:
            return False

    return False