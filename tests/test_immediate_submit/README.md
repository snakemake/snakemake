```bash
snakemake --profile slurm all
touch output/output.touch
```
- only after create `output/output.touch`, snakemake can run as expected.
