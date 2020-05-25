rule a:
    input:
        "test.txt"
    log:
        "logs/a.log"
    conda:
        "envs/a.yaml"
	output:
    	"test.out"
    shell:
        "cat {input} > {output}"
