# ----------------------------------------------------- #
# Imports and constants                                 #
# ----------------------------------------------------- #

import csv
import os
import sys
from pathlib import Path
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

_SCRIPT_SOURCE_PATHS = [
    Path(workflow.source_path(f"../scripts/{module}.py"))
    for module in (
        "genome_context",
        "minimap2_presets",
    )
]
SCRIPTS_DIR = _SCRIPT_SOURCE_PATHS[0].parent
sys.path.insert(0, str(SCRIPTS_DIR))

from genome_context import quast_tool, reject_meta_extra, resolve_config_genome_context
from minimap2_presets import compose_minimap2_extra, validate_minimap2_extra

SAMPLE_ID_PATTERN = r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}"


def expand_configured_path(value):
    """Expand environment variables and ``~`` without rebasing relative paths."""
    if value is None or str(value).strip() == "":
        return ""
    return os.path.expanduser(os.path.expandvars(str(value)))


SAMPLES_SCHEMA = workflow.current_basedir.join(
    "../schemas/samples.schema.yaml"
).get_path_or_uri(secret_free=False)


# Validate raw path values before expanding them. This keeps type errors at the
# configuration boundary instead of passing them to os.path expansion.
validate(config, schema="../schemas/config.schema.yaml")

try:
    GENOME_CONTEXT = resolve_config_genome_context(config)
    reject_meta_extra(config.get("metaspades_extra", ""), "metaspades_extra")
    reject_meta_extra(config.get("flye_extra", ""), "flye_extra")
    validate_minimap2_extra(config.get("minimap2_extra", ""))
except ValueError as error:
    raise WorkflowError(str(error)) from error

# ----------------------------------------------------- #
# Sample sheet loading                                  #
# ----------------------------------------------------- #


def _configured_sample_sheet_path():
    value = config.get("sample_sheet")
    if value is None or str(value).strip() == "":
        raise WorkflowError(
            "Assembly sample metadata requires a non-empty `sample_sheet` "
            "configuration value."
        )
    path = expand_configured_path(value)
    config["sample_sheet"] = path
    return path


def _load_sample_record(path):
    sample_sheet = Path(path)
    if not sample_sheet.exists():
        raise WorkflowError(f"Sample sheet does not exist: {path}")

    try:
        with sample_sheet.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t", strict=True)
            fieldnames = [field.strip() for field in (reader.fieldnames or [])]
            if not fieldnames or any(not field for field in fieldnames):
                raise WorkflowError(
                    f"Sample sheet must have a non-empty TSV header: {path}"
                )

            duplicate_headers = sorted(
                {field for field in fieldnames if fieldnames.count(field) > 1}
            )
            if duplicate_headers:
                raise WorkflowError(
                    "Sample sheet contains duplicate header(s): "
                    + ", ".join(duplicate_headers)
                )
            reader.fieldnames = fieldnames

            rows = []
            for row in reader:
                if None in row:
                    raise WorkflowError(
                        f"Sample sheet row {reader.line_num} has more fields than the header"
                    )
                if any(value is None for value in row.values()):
                    raise WorkflowError(
                        f"Sample sheet row {reader.line_num} has fewer fields than the header"
                    )

                # Sample identifiers are deliberately preserved byte-for-byte.
                # Established categorical values remain case-insensitive.
                raw_row = {
                    key: (
                        value
                        if key == "sample"
                        else (
                            value.strip().lower()
                            if key in {"assembly_type", "long_read_type"}
                            else value.strip()
                        )
                    )
                    for key, value in row.items()
                }
                validate(raw_row, schema=SAMPLES_SCHEMA)

                # Expand paths only after the complete row, including blank
                # optional cells, has passed schema validation.
                normalized = dict(raw_row)
                for column in ["short_read1", "short_read2", "long_reads"]:
                    if normalized[column]:
                        normalized[column] = expand_configured_path(normalized[column])
                rows.append(normalized)
    except csv.Error as error:
        raise WorkflowError(f"Malformed TSV sample sheet {path}: {error}") from error

    if len(rows) != 1:
        raise WorkflowError(
            f"Sample sheet must contain exactly one sample row; found {len(rows)}: "
            f"{path}"
        )

    return rows[0]


# ----------------------------------------------------- #
# Lazy normalized singleton state                       #
# ----------------------------------------------------- #


_SAMPLE_RECORD = None
SAMPLE_PATTERN = SAMPLE_ID_PATTERN
RESULTS_ROOT = f"results/{{sample,{SAMPLE_PATTERN}}}"
LOGS_ROOT = "<results>/logs"


def _sample_record():
    global _SAMPLE_RECORD
    if _SAMPLE_RECORD is None:
        _SAMPLE_RECORD = _load_sample_record(_configured_sample_sheet_path())
    return _SAMPLE_RECORD


# ----------------------------------------------------- #
# Sample sheet accessors                                #
# ----------------------------------------------------- #


def _record_for_sample(sample):
    record = _sample_record()
    configured_sample = record["sample"]
    if sample != configured_sample:
        raise WorkflowError(
            f"Unknown Assembly sample '{sample}'; the configured sample is "
            f"'{configured_sample}'. Exactly one matching sample is available per run."
        )
    return record


def _sample_sheet_input(wildcards):
    _record_for_sample(wildcards.sample)
    return _configured_sample_sheet_path()


SAMPLE_SHEET = _sample_sheet_input


def _value_from_samplesheet(sample, column):
    return _record_for_sample(sample)[column] or None


def get_short_reads(sample):
    read1 = _value_from_samplesheet(sample, "short_read1")
    read2 = _value_from_samplesheet(sample, "short_read2")
    if not read1 or not read2:
        raise WorkflowError(
            f"Sample '{sample}' requires both short_read1 and short_read2."
        )
    return [read1, read2]


def get_long_reads(sample):
    long_reads = _value_from_samplesheet(sample, "long_reads")
    if not long_reads:
        raise WorkflowError(f"Sample '{sample}' requires long_reads.")
    return long_reads


def get_long_read_type(sample):
    seq_type = _value_from_samplesheet(sample, "long_read_type")
    if not seq_type:
        raise WorkflowError(f"Sample '{sample}' requires long_read_type.")
    return seq_type


def _require_long_minimap2_sample(sample):
    assembly_type = _value_from_samplesheet(sample, "assembly_type")
    if assembly_type != "long":
        raise WorkflowError(
            "Minimap2 mapping is only available for samples with "
            f"`assembly_type=long`; sample '{sample}' has "
            f"`assembly_type={assembly_type}`. Use the BWA-MEM2 mapping outputs "
            "for short and hybrid assemblies."
        )


def get_minimap2_long_reads(sample):
    _require_long_minimap2_sample(sample)
    return get_long_reads(sample)


def get_minimap2_long_read_type(sample):
    _require_long_minimap2_sample(sample)
    return get_long_read_type(sample)


def get_hybrid_long_read_key(sample):
    seq_type = get_long_read_type(sample)
    if seq_type.startswith("nano-"):
        return "nanopore"
    if seq_type.startswith("pacbio-"):
        return "pacbio"

    raise WorkflowError(
        f"Cannot route unsupported long_read_type '{seq_type}' for hybrid sample "
        f"'{sample}'; expected a nano-* or pacbio-* value."
    )


def spades_inputs(wildcards):
    assembly_type = _value_from_samplesheet(wildcards.sample, "assembly_type")
    if assembly_type == "short":
        return {"reads": get_short_reads(wildcards.sample)}
    if assembly_type == "hybrid":
        return {
            "reads": get_short_reads(wildcards.sample),
            get_hybrid_long_read_key(wildcards.sample): get_long_reads(
                wildcards.sample
            ),
        }
    raise WorkflowError(
        "SPAdes assembly is only available for samples with "
        f"`assembly_type=short` or `assembly_type=hybrid`; sample "
        f"'{wildcards.sample}' has `assembly_type={assembly_type}`. Use the "
        "Flye assembly output for long assemblies."
    )


def _selected_route_outputs(sample):
    """Resolve the configured assembly, BAM, and BAI after rules are defined."""
    assembly_type = _record_for_sample(sample)["assembly_type"]
    if assembly_type in {"short", "hybrid"}:
        route = {
            "assembly": rules.assemble_spades.output.assembly,
            "bam": rules.bwamem2_map.output.bam,
            "bai": rules.bwamem2_map_index.output.bai,
        }
    elif assembly_type == "long":
        route = {
            "assembly": rules.assemble_flye_long.output.assembly,
            "bam": rules.minimap2_map.output.bam,
            "bai": rules.minimap2_map_index.output.bai,
        }
    else:
        raise WorkflowError(
            f"Cannot route sample '{sample}' with assembly_type {assembly_type!r}."
        )

    wildcards = {"sample": sample}
    return {
        name: str(output.apply_wildcards(wildcards)) for name, output in route.items()
    }


def selected_assembly(wildcards):
    return _selected_route_outputs(wildcards.sample)["assembly"]


# ----------------------------------------------------- #
# Final workflow targets                                #
# ----------------------------------------------------- #


def final_targets(wildcards):
    sample = _sample_record()["sample"]
    route = _selected_route_outputs(sample)
    targets = [
        route["assembly"],
        str(rules.quast.output.report_txt.apply_wildcards({"sample": sample})),
        route["bam"],
        route["bai"],
    ]

    return targets


ALL_TARGETS = final_targets
