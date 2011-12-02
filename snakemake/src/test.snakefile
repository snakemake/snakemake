import subprocess
reads = ["5100", "A1"]

rule mapreads:
	input: ["sample.{}.fastq".format(xy) for xy in reads] + ["hg19.fasta"]
	output: "mappped.bam"
	run:
		open(output[0], "w").write("test")

rule prepare_fasta:
	output: "sample.5100.fastq"
	run:
		open(output[0], "w").write("test")

rule prepare_fasta:
	output: "sample.A1.fastq"
	run:
		open(output[0], "w").write("test")

rule create_hg19:
	output: "hg19.fasta"
	run:
		open(output[0], "w").write("test")
