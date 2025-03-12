println("Julia script executing!")

@assert snakemake.config["test"] == true
@assert snakemake.config["foo' bar"] == "let's go"
@assert snakemake.params["integer"] == 123
@assert snakemake.params["astring"] == "foo\n'\\\" "
@assert snakemake.output[1] == "julia.out"
@assert snakemake.input[1] == "test.in"
@assert snakemake.input["named_input"] == "test.in"

f = open(snakemake.output[1], "w")
println(f, "Julia test succeded!")
close(f)
