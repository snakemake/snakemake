rule a:
    container:
        "docker://rule_image"
    output:
        "result.txt",
    shell:
        "echo module > {output}"
