"""Snakemake wrapper for metaspades."""

__author__ = "Silas Kieser @silask"
__copyright__ = "Copyright 2021, Silas Kieser"
__email__ = "silas.kieser@gmail.com"
__license__ = "MIT"

import os
import shutil
from snakemake.shell import shell

# infer output directory

if hasattr(snakemake.output, "dir"):
    output_dir = snakemake.output.dir

else:
    # get output_dir file from output
    if hasattr(snakemake.output, "contigs"):
        output_file = snakemake.output.contigs
    elif hasattr(snakemake.output, "scaffolds"):
        output_file = snakemake.output.scaffolds
    else:
        output_file = snakemake.output[0]

    output_dir = os.path.split(output_file)[0]


# parse params
extra = snakemake.params.get("extra", "")
kmers = snakemake.params.get("k", "'auto'")
genome_context = snakemake.params.get("genome_context", "metagenome")
if genome_context not in {"metagenome", "isolate"}:
    raise ValueError(
        "`genome_context` must be exactly `metagenome` or `isolate`; "
        f"got {genome_context!r}."
    )
meta_arg = "--meta" if genome_context == "metagenome" else ""
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if hasattr(snakemake.resources, "mem_mb"):
    mem_gb = snakemake.resources.mem_mb // 1000
    memory_option = "--memory"
else:
    mem_gb = ""
    memory_option = ""

if not os.path.exists(os.path.join(output_dir, "params.txt")):
    # parse short reads
    if hasattr(snakemake.input, "reads"):
        reads = snakemake.input.reads
    else:
        reads = snakemake.input

    assert (
        len(reads) > 1
    ), "Metaspades needs a paired end library. This means you should supply at least 2 fastq files in the rule input."

    assert (
        type(reads[0]) == str
    ), f"Metaspades allows only 1 library. Therefore reads need to be strings got {reads}"

    paired_1_option = "--pe1-1"
    paired_1 = reads[0]
    paired_2_option = "--pe1-2"
    paired_2 = reads[1]

    merged = reads[2] if len(reads) >= 3 and reads[2] else ""
    merged_option = "--pe1-m" if merged else ""

    singleton = reads[3] if len(reads) >= 4 and reads[3] else ""
    singleton_option = "--pe1-s" if singleton else ""

    # parse long reads
    pacbio = getattr(snakemake.input, "pacbio", "") or ""
    pacbio_option = "--pacbio" if pacbio else ""
    nanopore = getattr(snakemake.input, "nanopore", "") or ""
    nanopore_option = "--nanopore" if nanopore else ""

    shell(
        "export OMP_NUM_THREADS={snakemake.threads}; spades.py {meta_arg} "
        " --threads {snakemake.threads} "
        " {memory_option} {mem_gb:q} "
        " -o {output_dir:q} "
        " -k {kmers:q} "
        " {paired_1_option} {paired_1:q} "
        " {paired_2_option} {paired_2:q} "
        " {merged_option} {merged:q} "
        " {singleton_option} {singleton:q} "
        " {pacbio_option} {pacbio:q} "
        " {nanopore_option} {nanopore:q} "
        " {extra} "
        " > {snakemake.log[0]:q} 2>&1 "
    )


else:
    # params.txt file exitst already I restart from previous run

    shell(
        "echo '\n\nRestart Spades \n Remove pipline_state file copy files to force copy files if necessary.' >> {snakemake.log[0]:q}"
    )

    shell(
        "rm -f {output_dir:q}/pipeline_state/stage_*_copy_files "
        "2>> {snakemake.log[0]:q}"
    )

    shell(
        "export OMP_NUM_THREADS={snakemake.threads}; spades.py "
        " --restart-from last "
        " --threads {snakemake.threads} "
        " {memory_option} {mem_gb:q} "
        " -o {output_dir:q} "
        " >> {snakemake.log[0]:q} 2>&1 "
    )


# Rename/ move output files

Output_key_mapping = {
    "contigs": "contigs.fasta",
    "scaffolds": "scaffolds.fasta",
    "graph": "assembly_graph_with_scaffolds.gfa",
}

has_named_output = False
for key in Output_key_mapping:
    if hasattr(snakemake.output, key):
        has_named_output = True
        file_produced = os.path.join(output_dir, Output_key_mapping[key])
        file_renamed = getattr(snakemake.output, key)

        if file_produced != file_renamed:
            shutil.move(file_produced, file_renamed)


if not has_named_output:
    file_produced = os.path.join(output_dir, "contigs.fasta")
    file_renamed = snakemake.output[0]

    if file_produced != file_renamed:
        shutil.move(file_produced, file_renamed)
