rule a:
    input:
        "test.csv"
    output:
        "test.out.csv"
    log:
        "logs/a.log"
    run:
        # this becomes too much and should be migrated into a script directive
        import pandas as pd

        t = pd.read_csv(snakemake.input[0])
        t.sort(inplace=True)

        t.loc["foo"] += 16
        t.loc["bar"] /= 237879.4

        t.to_csv(snakemake.output[0])