import csv
import io
import json
import subprocess
from collections import Counter
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
TEST_DIRECTORY = ROOT / ".test"
COLUMNS = [
    "sample",
    "assembly_type",
    "short_read1",
    "short_read2",
    "long_reads",
    "long_read_type",
]
CASES = {
    "short": {
        "config": TEST_DIRECTORY / "config" / "short.yaml",
        "sample": "short_sample",
        "assembly_rule": "assemble_spades",
        "assembly_dir": "assembly-spades",
        "mapping_rule": "bwamem2_map",
        "mapping_index_rule": "bwamem2_map_index",
        "mapping_dir": "mapping-bwamem2",
    },
    "long-nanopore": {
        "config": TEST_DIRECTORY / "config" / "long-nanopore.yaml",
        "sample": "long_sample",
        "assembly_rule": "assemble_flye_long",
        "assembly_dir": "assembly-flye",
        "mapping_rule": "minimap2_map",
        "mapping_index_rule": "minimap2_map_index",
        "mapping_dir": "mapping-minimap2",
    },
    "hybrid-pacbio": {
        "config": TEST_DIRECTORY / "config" / "hybrid-pacbio.yaml",
        "sample": "hybrid_sample",
        "assembly_rule": "assemble_spades",
        "assembly_dir": "assembly-spades",
        "mapping_rule": "bwamem2_map",
        "mapping_index_rule": "bwamem2_map_index",
        "mapping_dir": "mapping-bwamem2",
    },
}


def _snakemake(case, *arguments, configfile=None, directory=TEST_DIRECTORY):
    selected_config = configfile or CASES[case]["config"]
    return subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(directory),
            "--configfile",
            str(selected_config),
            "--nolock",
            *arguments,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def _successful_result(case, *arguments, **kwargs):
    result = _snakemake(case, *arguments, **kwargs)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


@pytest.mark.parametrize("case", CASES)
def test_one_sample_case_has_selected_route(case):
    expected = CASES[case]
    result = _successful_result(case, "--d3dag")
    graph = json.loads(result.stdout)
    rules = {node["value"]["rule"] for node in graph["nodes"]}

    assert expected["assembly_rule"] in rules
    assert expected["mapping_rule"] in rules
    assert expected["mapping_index_rule"] in rules

    opposite_rules = (
        {"assemble_spades", "bwamem2_map", "bwamem2_map_index"}
        if case == "long-nanopore"
        else {"assemble_flye_long", "minimap2_map", "minimap2_map_index"}
    )
    assert rules.isdisjoint(opposite_rules)


@pytest.mark.parametrize("case", CASES)
def test_one_sample_case_has_exact_mapping_targets(case):
    expected = CASES[case]
    sample = expected["sample"]
    result = _successful_result(case, "--summary")
    rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter="\t"))
    outputs = {row["output_file"] for row in rows}

    assert f"results/{sample}/{expected['assembly_dir']}/assembly.fasta" in outputs
    assert f"results/{sample}/qc-quast/report.txt" in outputs
    assert f"results/{sample}/{expected['mapping_dir']}/alignment.bam" in outputs
    assert f"results/{sample}/{expected['mapping_dir']}/alignment.bam.bai" in outputs
    opposite_mapping = (
        "mapping-minimap2"
        if expected["mapping_dir"] == "mapping-bwamem2"
        else "mapping-bwamem2"
    )
    assert not any(
        output.startswith(f"results/{sample}/{opposite_mapping}/") for output in outputs
    )
    opposite_assembler = (
        "assembly-flye"
        if expected["assembly_dir"] == "assembly-spades"
        else "assembly-spades"
    )
    assert not any(
        output.startswith(f"results/{sample}/{opposite_assembler}/")
        for output in outputs
    )
    assert not any(
        output.startswith(f"results/{sample}/{obsolete}/")
        for output in outputs
        for obsolete in ("assembly", "mapping", "quast")
    )


@pytest.mark.parametrize(
    ("case", "sample", "assembly_type"),
    [
        ("short", "short_sample", "short"),
        ("hybrid-pacbio", "hybrid_sample", "hybrid"),
    ],
)
def test_non_long_minimap2_target_has_actionable_error(case, sample, assembly_type):
    result = _snakemake(
        case, "--dry-run", f"results/{sample}/mapping-minimap2/alignment.bam"
    )
    output = result.stdout + result.stderr
    diagnostic = (
        "Minimap2 mapping is only available for samples with "
        f"`assembly_type=long`; sample '{sample}' has "
        f"`assembly_type={assembly_type}`. Use the BWA-MEM2 mapping outputs "
        "for short and hybrid assemblies."
    )

    assert result.returncode != 0
    assert diagnostic in output


def _valid_short(sample="sample"):
    return {
        "sample": sample,
        "assembly_type": "short",
        "short_read1": "reads_R1.fastq",
        "short_read2": "reads_R2.fastq",
        "long_reads": "",
        "long_read_type": "",
    }


def _write_sheet(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _run_temporary_sheet(
    tmp_path, rows=None, sheet_text=None, arguments=("--dry-run",)
):
    sample_sheet = tmp_path / "samples.tsv"
    if sheet_text is None:
        _write_sheet(sample_sheet, rows or [])
    else:
        sample_sheet.write_text(sheet_text, encoding="utf-8")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample_sheet": str(sample_sheet),
                "genome_context": "metagenome",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    result = _snakemake("short", *arguments, configfile=config_path, directory=tmp_path)
    return result, result.stdout + result.stderr


def _short_with_fixture_reads(sample, assembly_type="short"):
    row = _valid_short(sample)
    row.update(
        {
            "assembly_type": assembly_type,
            "short_read1": str(TEST_DIRECTORY / "data" / "short_sample_R1.fastq"),
            "short_read2": str(TEST_DIRECTORY / "data" / "short_sample_R2.fastq"),
        }
    )
    return row


def test_single_sample_record_is_normalized_before_route_selection(tmp_path):
    result, output = _run_temporary_sheet(
        tmp_path,
        rows=[_short_with_fixture_reads("normalized_sample", " SHORT ")],
        arguments=("--d3dag",),
    )

    assert result.returncode == 0, output
    graph = json.loads(result.stdout)
    rules = Counter(node["value"]["rule"] for node in graph["nodes"])
    assert rules["assemble_spades"] == 1
    assert rules["bwamem2_map"] == 1
    assert rules["minimap2_map"] == 0


def test_unknown_schema_compatible_sample_has_actionable_error(tmp_path):
    row = _short_with_fixture_reads("sample.one")
    configured_target = "results/sample.one/assembly-spades/assembly.fasta"
    result, output = _run_temporary_sheet(
        tmp_path, rows=[row], arguments=("--dry-run", configured_target)
    )
    assert result.returncode == 0, output

    nonconfigured_target = "results/sampleXone/assembly-spades/assembly.fasta"
    result, output = _run_temporary_sheet(
        tmp_path, rows=[row], arguments=("--dry-run", nonconfigured_target)
    )
    assert result.returncode != 0
    assert (
        "Unknown Assembly sample 'sampleXone'; the configured sample is "
        "'sample.one'. Exactly one matching sample is available per run."
    ) in output


@pytest.mark.parametrize(
    "rows",
    [
        [],
        [_valid_short("sample_one"), _valid_short("sample_two")],
    ],
)
def test_sample_sheet_requires_exactly_one_row(tmp_path, rows):
    result, output = _run_temporary_sheet(tmp_path, rows=rows)

    assert result.returncode != 0
    assert (
        f"Sample sheet must contain exactly one sample row; found {len(rows)}" in output
    )


def test_malformed_row_shape_precedes_cardinality_error(tmp_path):
    header = "\t".join(COLUMNS)
    valid = "sample_one\tshort\treads_R1.fastq\treads_R2.fastq\t\t"
    malformed = "sample_two\tshort\treads_R1.fastq\treads_R2.fastq\t\t\textra"
    result, output = _run_temporary_sheet(
        tmp_path, sheet_text=f"{header}\n{valid}\n{malformed}\n"
    )

    assert result.returncode != 0
    assert "row 3 has more fields than the header" in output
    assert "exactly one sample row" not in output


def test_schema_error_precedes_cardinality_error(tmp_path):
    invalid = {**_valid_short("sample_two"), "assembly_type": "unsupported"}
    result, output = _run_temporary_sheet(
        tmp_path, rows=[_valid_short("sample_one"), invalid]
    )

    assert result.returncode != 0
    assert "assembly_type" in output
    assert "unsupported" in output
    assert "exactly one sample row" not in output
