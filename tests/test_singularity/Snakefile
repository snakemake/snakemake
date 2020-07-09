

rule a:
    input:
        "test.txt"
    output:
        "test.out"
    singularity:
        "docker://bash"
    shell:
        'ls {input}; echo "test" > {output}'


rule b:
    output:
        "test.txt"
    singularity:
        "docker://python:3.8.3"
    script:
        "script.py"


rule c:
    output:
        "invalid.txt"
    singularity:
        "docker://invalid/image url even with whitespace"
    shell:
        "touch {output}"
