#!/bin/sh
#rule: {job}
#input: {job.input}
#output: {job.output}
{self.workflow.snakemakepath} --snakefile {self.workflow.snakefile} \
--force -j{self.cores} \
--directory {workdir} --nocolor --notemp --quiet --nolock {job.output} \
> /dev/null && touch "{jobfinished}" || touch "{jobfailed}"
exit 0
