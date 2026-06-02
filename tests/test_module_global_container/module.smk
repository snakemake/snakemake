container: "docker://module_image"


rule a:
    output:
        "result.txt",
    shell:
        "echo module > {output}"
