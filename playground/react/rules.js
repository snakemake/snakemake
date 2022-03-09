var rules = {
    "a": {
        "input": [
            "x.txt",
            "y.txt"
        ],
        "output": [
            "z.txt",
            "y2.txt"
        ],
        "conda": [
            "bwa =0.74"
        ],
        "container": "docker://miniconda",
        "code": [
            "codeblock1(1.4)",
            "codeblock2(1.5)"
        ]
    }
}