import os,subprocess,time
reads = ["5100", "A1"]

workdir: "test"

rule mapreads:
	input: ["sample.{}.fastq".format(xy) for xy in reads],
	       "hg19.fasta",
	       "5100.fasta"
	output: "mappped.bam"
	run:
		open(output[0], "w").write("test")
		shell("ech `echo $PATH`")

rule prepare_fasta:
	output: "sample.5100.fastq", "5100.fasta"
	run:
		open(output[0], "w").write("test")
		for i in range(3):
			print("X")
			time.sleep(1)
		open(output[1], "w").write("test")


rule prepare_fasta2:
	output: "sample.A1.fastq"
	run:
		open(output[0], "w").write("test")
		for i in range(3):
			print("Y")
			time.sleep(1)

rule create_hg19:
	output: "hg19.fasta"
	run:
		open(output[0], "w").write("test")
