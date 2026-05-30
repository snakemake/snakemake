rule all:
    input:
        "res.txt",


a = {"b": "!!!!"}


rule t1:
    input:
        file="../a.txt",
    output:
        res="tmp.txt",
    run:
        b = a["b"]
        shell("echo {b} >> {output.res}")


rule t2:
    input:
        file=rules.t1.output.res,
        file2=[],
    output:
        res1="res.txt",
    shell:
        "cat {input.file} {input.file2} > {output.res1}"
