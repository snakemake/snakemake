shell.executable("bash")

rule All: input: "C"

rule A:
    output: "outputA/A1", "outputA/A2", directory("outputA")
    shell: "mkdir -p outputA; touch outputA/A1; touch outputA/A2"

rule B:
    input: "outputA"
    output: "B"
    shell: "touch B"

rule C:
    input: test="outputA/A1", testb="B"
    output: "C"
    shell: "cp {input.test} C"
