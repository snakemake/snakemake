#!/bin/sh
# properties = {properties}
{workflow.snakemakepath} --snakefile {workflow.snakefile} \
--force -j{cores} \
--directory {workdir} --nocolor --notemp --quiet --nolock {job.output} \
&& touch "{jobfinished}" || touch "{jobfailed}"
exit 0
