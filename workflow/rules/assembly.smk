# ----------------------------------------------------- #
# Assembly workflow                                     #
# ----------------------------------------------------- #

import os


rule assemble_spades:
    input:
        unpack(spades_inputs),
    output:
        assembly="<results>/assembly-spades/assembly.fasta",
    log:
        "<results>/assembly-spades/run.log",
    conda:
        "../envs/spades.yaml"
    threads: config.get("metaspades_threads", 4)
    resources:
        mem_mb=lambda wc, threads: config.get("metaspades_memory_gb", 110) * 1000,
    params:
        k=config.get("metaspades_kmers", "auto"),
        genome_context=GENOME_CONTEXT,
        extra=config.get("metaspades_extra", ""),
    message:
        """--- Running SPAdes assembly for sample {wildcards.sample}."""
    script:
        "../scripts/spades.py"


rule assemble_flye_long:
    input:
        reads=lambda wc: get_long_reads(wc.sample),
    output:
        assembly="<results>/assembly-flye/assembly.fasta",
    log:
        "<results>/assembly-flye/run.log",
    conda:
        "../envs/flye.yaml"
    threads: config.get("flye_threads", 4)
    params:
        out_dir=lambda wc, output: os.path.dirname(output.assembly),
        seq_type=lambda wc: get_long_read_type(wc.sample),
        genome_context=GENOME_CONTEXT,
        extra=config.get("flye_extra", ""),
    message:
        """--- Running Flye assembly for long-read sample {wildcards.sample}."""
    script:
        "../scripts/flye.py"
