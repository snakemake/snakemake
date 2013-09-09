#!/bin/sh
#rule: {job}
#input: {job.input}
#output: {job.output}
#jobid: {jobid}
{workflow.snakemakepath} --snakefile {workflow.snakefile} \
--force -j{cores} \
--directory {workdir} --nocolor --notemp --quiet --nolock {job.output} \
&& touch "{jobfinished}" || touch "{jobfailed}"
exit 0
