checkpoint A:
    output:
        touch("a")


use rule A as B with:
    output:
        touch("b")

checkpoint a_A:
    output:
        touch("c")
