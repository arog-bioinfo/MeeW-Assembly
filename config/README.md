## Workflow overview

This workflow routes exactly one configured metagenomic sample by assembly type and runs assembly, assembly quality assessment, and read-to-assembly mapping.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Route the configured sample by `assembly_type` (`short`, `long`, `hybrid`)
2. Assemble using:
   - `metaSPAdes` for `short`
   - local `Flye` wrapper for `long`
   - `metaSPAdes` hybrid mode for `hybrid`
3. Run `QUAST` for all assemblies
4. Map reads back to assemblies:
   - `BWA-MEM2` for `short`
   - `minimap2` for `long`
   - `BWA-MEM2` for `hybrid`; its long reads assist assembly but are not mapped, and no hybrid minimap2 route is provided

## Running the workflow

### Configuration organization

`config/config.yaml` is organized into commented sections that match the major workflow stages:

- Workflow Input: sample sheet location
- Execution Resources: tool threads and memory settings
- Assembly: metaSPAdes and Flye options
- Assembly Quality Assessment: QUAST options
- Read Mapping: BWA-MEM2 and minimap2 options

### Input data

The workflow consumes user-provided reads listed in the sample sheet referenced by `sample_sheet` in `config/config.yaml`. The sheet must contain its header and exactly one sample row; header-only and multi-row sheets are rejected.
The sample-sheet path and read paths expand `~` and environment variables.
Relative paths remain relative to the Snakemake working directory (including
when `--directory` is used); they are not rebased to the YAML or TSV location.
Use exactly one of the following mode-specific row shapes per run. Do not combine these examples into one sheet.

Short-read run:

| sample       | assembly_type | short_read1                | short_read2                | long_reads | long_read_type |
| ------------ | ------------- | -------------------------- | -------------------------- | ---------- | -------------- |
| short_sample | short         | data/short_sample_R1.fastq | data/short_sample_R2.fastq |            |                |

Long Nanopore run:

| sample      | assembly_type | short_read1 | short_read2 | long_reads             | long_read_type |
| ----------- | ------------- | ----------- | ----------- | ---------------------- | -------------- |
| long_sample | long          |             |             | data/long_sample.fastq | nano-raw       |

Hybrid PacBio run:

| sample        | assembly_type | short_read1                 | short_read2                 | long_reads                    | long_read_type |
| ------------- | ------------- | --------------------------- | --------------------------- | ----------------------------- | -------------- |
| hybrid_sample | hybrid        | data/hybrid_sample_R1.fastq | data/hybrid_sample_R2.fastq | data/hybrid_sample_long.fastq | pacbio-raw     |

Rules for required fields:

- `short`: requires `short_read1` and `short_read2`
- `long`: requires `long_reads` and `long_read_type`
- `hybrid`: requires `short_read1`, `short_read2`, `long_reads`, and `long_read_type`

Valid `long_read_type` values are Flye-compatible: `nano-raw`, `nano-corr`, `nano-hq`, `pacbio-raw`, `pacbio-corr`, and `pacbio-hifi`.
For hybrid metaSPAdes assembly, the long-read input key is inferred from `long_read_type`: `nano-*` values use `nanopore`, and `pacbio-*` values use `pacbio`.
For minimap2 mapping of long assemblies, the preset is derived per sample: all `nano-*` values use `map-ont`, `pacbio-raw` and `pacbio-corr` use `map-pb`, and `pacbio-hifi` uses `map-hifi`. `minimap2_extra` accepts additional trusted minimap2 arguments, but must not contain any `-x` preset option.

`assembly_type` and `long_read_type` are normalized by trimming whitespace and lowercasing before routing and validation.

The configured sample row is loaded, validated, normalized, and cached only when sample metadata is first requested. Listing or selectively importing rules does not read the sample sheet. An importing workflow that supplies sample-dependent behavior directly can use Snakemake's module validation controls; standalone configuration continues to require `sample_sheet`.

### BWA-MEM2 BAM contract

For `short` and `hybrid` samples, BWA-MEM2 always produces the coordinate-sorted BAM `results/{sample}/mapping-bwamem2/alignment.bam`, followed by its index at `results/{sample}/mapping-bwamem2/alignment.bam.bai`. This coordinate-sorted, indexed BAM is the mapping input expected by MeeW-Binning. Long samples instead use the corresponding BAM and index under `results/{sample}/mapping-minimap2/`.

`bwamem2_sort` accepts `samtools` or `picard`; omitting it uses `samtools`. `bwamem2_sort_order` accepts only `coordinate`, which is also its default. Unsorted (`none`) and query-name-sorted (`queryname`) configurations are rejected during configuration validation.

Binning owns any derived tool-specific inputs such as `depth` and `kmer`; standalone Binning input remains a separately user-authored sheet.
