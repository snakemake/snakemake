storage:
    provider="http",
    keep_local=True

storage http_ret:
    provider="http",
    keep_local=True,
    retrieve=False,


rule retrieve_default:
    input: storage.http("https://github.com/snakemake/snakemake/blob/main/images/logo.png")
    output: "retrieve_default.flag"
    shell: "touch {output}"

rule retrieve_true:
    input: storage.http("https://github.com/snakemake/snakemake/blob/main/images/logo.png", retrieve=True)
    output: "retrieve_true.flag"
    shell: "touch {output}"

rule retrieve_false:
    input: storage.http("https://github.com/snakemake/snakemake/blob/main/images/logo.png", retrieve=False)
    output: "retrieve_false.flag"
    shell: "touch {output}"

rule retrieve_false_directive:
    input: storage.http_ret("https://github.com/snakemake/snakemake/blob/main/images/logo.png")
    output: "retrieve_false_directive.flag"
    shell: "touch {output}"
