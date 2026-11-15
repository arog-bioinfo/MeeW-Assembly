# Snakemake pipeline: MeeW-Assembly

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A59.12.0-039475?logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/conda/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<!-- TEMPORARY CATALOG OPT-OUT: later replace every standalone "pipeline" with the standard term formed by "work" + "flow", then remove this comment. -->

A Snakemake pipeline for metagenome assembly, assembly QC, and read-to-assembly mapping.

> [!WARNING]
> **Research preview:** This pipeline is under active development and is not a deployment-ready release. Validate its results and resource requirements for your use case.

## Usage

Each invocation accepts exactly one sample-sheet row in `short`, `long`, or `hybrid` mode. Short reads are assembled with metaSPAdes and mapped with BWA-MEM2; long reads are assembled with Flye and mapped with minimap2. Hybrid long reads assist metaSPAdes assembly, but only the short reads are mapped with BWA-MEM2. QUAST assesses every assembly. The current routes require no reference database.

Install Snakemake `>=9.12.0` and Conda, then edit [`config/config.yaml`](config/config.yaml) and its one-row sample sheet. Supply reads for the selected mode and review resource settings before execution; the shipped metaSPAdes memory request is substantial. Results below `results/{sample}/` contain the selected assembler output, QUAST report, and coordinate-sorted indexed mapping BAM used by downstream binning.

See the [configuration guide](config/README.md) for the sample-sheet contract, mode-specific routing, resources, and output paths. Run the commands below from the repository root.

## Deployment options

Conda is the supported software deployment method. After configuring `config/config.yaml`, inspect the planned jobs and then execute them with:

```bash
snakemake --cores 2 --dry-run
snakemake \
  --cores 2 \
  --software-deployment-method conda
```

Adjust cores and resource overrides for the data and execution environment. This pipeline does not supply a container deployment.

## Execution profiles

No pipeline-specific profile is supplied. A user-maintained [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) may provide executor and site-specific resource settings.

## License

The pipeline is distributed under the [MIT License](LICENSE).

## References

No module-specific DOI is declared. Cite the repository version used, the selected route's tools—[SPAdes](https://github.com/ablab/spades), [Flye](https://github.com/mikolmogorov/Flye), [QUAST](https://github.com/ablab/quast), [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2), and [minimap2](https://github.com/lh3/minimap2)—and Snakemake:

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
