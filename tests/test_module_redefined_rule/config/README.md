
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, `alias`, `platform`, `datatype`, `calling` and `group` have to be defined. 
* Samples within the same `group` can be referenced in a joint [Calling scenario](#calling-scenario) via their `alias`es.
* `alias`es represent the name of the sample within its group. They are meant to be some abstract description of the sample type to be used in the [Calling scenario](#calling-scenario), and should thus be used consistently across groups. A classic example would be a combination of the `tumor` and `normal` aliases.
* The `platform` column needs to contain the used sequencing plaform (one of 'CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'ONT', 'PACBIO’).
* The purity column is required when being used with the default scenario. If it is unknown, it can be set to `1.0`.
* The same `sample_name` entry can be used multiple times within a `samples.tsv` sample sheet, with only the value in the `group` column differing between repeated rows. This way, you can use the same sample for variant calling in different groups, for example if you use a panel of normal samples when you don't have matched normal samples for tumor variant calling.
* The `datatype` column specifies what kind of data each sample corresponds to. This can either be `rna` or `dna`.
* The `calling` column sets the kind of analysis to be performed. This can be either `fusions`, `variants` or both (comma separated). Fusion calling is still under developement and should be considered as experimental. 

If mutational burdens shall be estimated for a sample, the to be used ``events`` from the calling scenario (see below) have to be specified in an additional column ``mutational_burden_events``. Multiple events have to be separated by commas within that column.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# Unit sheet

For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`.
* Each unit has a `unit_name`. This can be a running number, or an actual run, lane or replicate id.
* Each unit has a `sample_name`, which associates it with the biological sample it comes from. This information is used to merged all the units of a sample before read mapping and duplicate marking.
* For each unit, you need to specify either of these columns:
  * `fq1` only for single end reads. This can point to any FASTQ file on your system
  * `fq1` and `fq2` for paired end reads. These can point to any FASTQ files on your system
  * `sra` only: specify an SRA (sequence read archive) accession (starting with e.g. ERR or SRR). The pipeline will automatically download the corresponding paired end reads from SRA.
  * If both local files (`fq1`, `fq2`) and SRA accession (`sra`) are available, the local files will be used.
* Define adapters in the `adapters` column, by putting [cutadapt arguments](https://cutadapt.readthedocs.org) in quotation marks (e.g. `"-a ACGCGATCG -A GCTAGCGTACT"`).

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# Calling scenario

Varlociraptor supports integrated uncertainty aware calling and filtering of variants for arbitrary scenarios. These are defined as so-called scenarios, via a [variant calling grammar](https://varlociraptor.github.io/docs/calling#generic-variant-calling).
* For each group, a scenario is rendered via [YTE](https://yte-template-engine.github.io).
* Therefore, edit the template scenario (`scenario.yaml`) according to your needs. The sample sheet is available for YTE rendering as a pandas data frame in the variable `samples`. This allows to customize the scenario according to the contents of the sample sheet. You can therefore add additional columns to the sample sheet (e.g. purity) and access them in the scenario template, in order to pass the information to Varlociraptor.
* Example scenarios for various use cases can be found in the [scenario catalog](https://varlociraptor.github.io/varlociraptor-scenarios).

# Primer trimming

For panel data the pipeline allows trimming of amplicon primers on both ends of a fragment but also on a single end only. 
In case of single end primers these are supposed to be located at the left end of a read.
When primer trimming is enabled, primers have to be defined either directly in the `config.yaml` or in a seperate tsv-file.
Defining primers directly in the config file is prefered when all samples come from the same primer set.
In case of different panels, primers have to be set panel-wise in a seperate tsv-file.
For each panel the following columns need to be set: `panel`, `fa1` and `fa2` (optional).
Additionally, for each sample the corresponding panel must be defined in `samples.tsv` (column `panel`).
If a panel is not provided for a sample, trimming will not be performed on that sample. 
For single primer trimming only, the first entry in the config (respective in the tsv file) needs to be defined.

# Annotating UMIS

For annotating UMIs two additional columns in `sample.tsv` must be set:
* `umi_read`: this can be either of the following options:
  * `fq1` if the UMIs are part of read 1
  * `fq2` if the UMIs are part of read 2
  * `both` if there are UMIs in both paired end reads
  * the path to an additional fastq file containing just the UMI of each fragment in fq1 and fq2 (with the same read names)
* `umi_read_structure`: A read structure defining the UMI position in each UMI record (see https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures). If `both` reads contain a UMI, specify a read structure for both with whitespace in between (for example, `8M142T 8M142T`). In case a separate fastq file only containg UMI sequences is set the read structure needs to be `+M`.
Read names of UMI records must match the corresponding records of the sample fastqs.

