
reads = ["5100", "A1"]

rule mapreads:
	input: ["sample.{}.fastq".format(xy) for xy in reads] + ["hg19.fasta"]
	output: "mappped.bam", "logfile"
	run:
		print(output)

rule prepare_fasta:
	output: "sample.{read}.fastq"
	run:
		open(output[0], "w").write("test")
		print("written fastq for {}".format(wildcards['read']))

rule create_hg19:
	output: "hg19.fasta"
	run:
		open(output[0], "w").write("test")