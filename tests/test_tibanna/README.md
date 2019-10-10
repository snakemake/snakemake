```
# the following has been performed by aws admin already
tibanna deploy_unicorn --usergroup=johannes --buckets=snakemake-tibanna-test,snakemake-tibanna-test2
```
```
# run the following to do a test run
pip install -U tibanna
python cleanup.py  # first delete output files that are already on s3
export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=tibanna_unicorn_johannes
snakemake --tibanna --use-conda --configfile=config.json --default-remote-prefix=snakemake-tibanna-test/1
```

