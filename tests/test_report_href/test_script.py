import textwrap

with open(snakemake.output[0], "w") as f:
    print(
        textwrap.dedent(f"""
        <html>
        <head>
        <title>Report</title>
        </head>
        <body>
            <a href={snakemake.report_href("test.html").url_args(foo=4).anchor("bar")}>Link to test.html</a>
            <a href={snakemake.report_href("subdir").child_path("subdir/test3.html")}>Link to subdir/test3.html</a>
        </body>
        </html>
          """
        ),
        file=f,
    )