import marimo

__generated_with = "0.19.4"
app = marimo.App(width="medium")


@app.cell
def _(snakemake):
    with open(snakemake.output[0], "w") as out_file:
        print("Snakemake + marimo!", file=out_file)
    return


if __name__ == "__main__":
    app.run()
