greeting = "hello"
message = greeting + "\n"
with open(snakemake.output[0], "w") as f:
    f.write(message)
