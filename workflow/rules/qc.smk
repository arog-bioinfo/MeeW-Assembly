# ----------------------------------------------------- #
# QC and mapping workflow                               #
# ----------------------------------------------------- #


rule quast:
    input:
        fasta=selected_assembly,
    output:
        report_txt="<results>/qc-quast/report.txt",
    log:
        "<results>/qc-quast/run.log",
    threads: config.get("quast_threads", 2)
    params:
        tool=quast_tool(GENOME_CONTEXT),
        extra=lambda wc: (
            lambda qe: " ".join(
                x
                for x in [
                    (
                        ""
                        if "--min-contig " in qe or "--min-contig=" in qe
                        else f"--min-contig {config.get('quast_min_contig',500)}"
                    ),
                    (
                        ""
                        if "--min-identity " in qe or "--min-identity=" in qe
                        else f"--min-identity {config.get('quast_min_identity',95.0)}"
                    ),
                    qe,
                ]
                if x
            )
        )(config.get("quast_extra", "").strip()),
    message:
        """--- Running QUAST on assembly for sample {wildcards.sample}."""
    wrapper:
        "v9.9.0/bio/quast"


rule bwamem2_index:
    input:
        selected_assembly,
    output:
        num="<results>/mapping-bwamem2/assembly.0123",
        amb="<results>/mapping-bwamem2/assembly.amb",
        ann="<results>/mapping-bwamem2/assembly.ann",
        bwt_2bit_64="<results>/mapping-bwamem2/assembly.bwt.2bit.64",
        pac="<results>/mapping-bwamem2/assembly.pac",
    log:
        "<results>/mapping-bwamem2/index.log",
    params:
        prefix=lambda wc, output: output.amb.rsplit(".amb", 1)[0],
        extra=config.get("bwamem2_index_extra", ""),
    message:
        """--- Building BWA-MEM2 index for sample {wildcards.sample}."""
    wrapper:
        "v9.8.0/bio/bwa-mem2/index"


rule bwamem2_map:
    input:
        idx=rules.bwamem2_index.output,
        reads=lambda wc: get_short_reads(wc.sample),
    output:
        bam="<results>/mapping-bwamem2/alignment.bam",
    log:
        "<results>/mapping-bwamem2/map.log",
    threads: config.get("bwamem2_threads", 4)
    params:
        sort=config.get("bwamem2_sort", "samtools"),
        sort_order=config.get("bwamem2_sort_order", "coordinate"),
        sort_extra=config.get("bwamem2_sort_extra", ""),
        extra=config.get("bwamem2_mem_extra", ""),
    message:
        """--- Mapping short reads with BWA-MEM2 for sample {wildcards.sample}."""
    wrapper:
        "v9.8.0/bio/bwa-mem2/mem"


rule bwamem2_map_index:
    input:
        bam="<results>/mapping-bwamem2/alignment.bam",
    output:
        bai="<results>/mapping-bwamem2/alignment.bam.bai",
    log:
        "<results>/mapping-bwamem2/index_bam.log",
    threads: 2
    message:
        """--- Indexing BWA-MEM2 BAM for sample {wildcards.sample}."""
    wrapper:
        "v9.8.0/bio/samtools/index"


rule minimap2_map:
    input:
        target=selected_assembly,
        query=lambda wc: get_minimap2_long_reads(wc.sample),
    output:
        bam="<results>/mapping-minimap2/alignment.bam",
    log:
        "<results>/mapping-minimap2/map.log",
    threads: config.get("minimap2_threads", 4)
    params:
        sorting="coordinate",
        extra=lambda wc: compose_minimap2_extra(
            get_minimap2_long_read_type(wc.sample), config.get("minimap2_extra", "")
        ),
    message:
        """--- Mapping long reads with minimap2 for sample {wildcards.sample}."""
    wrapper:
        "v9.8.0/bio/minimap2/aligner"


rule minimap2_map_index:
    input:
        bam="<results>/mapping-minimap2/alignment.bam",
    output:
        bai="<results>/mapping-minimap2/alignment.bam.bai",
    log:
        "<results>/mapping-minimap2/index_bam.log",
    threads: 2
    message:
        """--- Indexing minimap2 BAM for sample {wildcards.sample}."""
    wrapper:
        "v9.8.0/bio/samtools/index"
