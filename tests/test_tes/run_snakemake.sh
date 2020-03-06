sleep 30
cd "$(dirname "$0")"
python -m snakemake --snakefile Snakefile all
sleep 10