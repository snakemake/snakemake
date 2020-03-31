rule a:
    input:
        "test.csv"
    output:
        "test.out.csv"
    log:
        "logs/a.log"
    run:
        # a few lines of code are still ok
        import pandas as pd

        t = pd.read_csv(snakemake.input[0])
        t.sort(inplace=True)

        t.to_csv(snakemake.output[0])