import json
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator, ValidationError

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
CONFIG_SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"


def _run_snakemake(snakefile, workdir, *arguments):
    return subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(snakefile),
            "--directory",
            str(workdir),
            "--nolock",
            *arguments,
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_standalone_schema_requires_sample_sheet():
    schema = yaml.safe_load(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    validator = Draft202012Validator(schema)

    errors = list(validator.iter_errors({"genome_context": "metagenome"}))

    assert len(errors) == 1
    assert errors[0].validator == "required"
    assert "sample_sheet" in errors[0].message


def test_standalone_schema_rejects_removed_stage_sheets_setting():
    schema = yaml.safe_load(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config["stage_sheets"] = {"enabled": False}

    with pytest.raises(ValidationError):
        Draft202012Validator(schema).validate(config)


def test_selective_module_import_does_not_require_sample_metadata(tmp_path):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    harness = tmp_path / "Snakefile"
    harness.write_text(
        f"""module assembly:
    pathvars:
        results="results/{{sample}}/assembly",
    snakefile:
        {json.dumps(str(SNAKEFILE))}
    config:
        {{
            "genome_context": "metagenome",
        }}
    skip_validation:
        True

use rule assemble_spades from assembly as imported_assemble_spades with:
    input:
        reads=["reads_R1.fastq", "reads_R2.fastq"],

rule imported_target:
    default_target: True
    input:
        "results/imported/assembly/assembly-spades/assembly.fasta",
""",
        encoding="utf-8",
    )
    (workdir / "reads_R1.fastq").touch()
    (workdir / "reads_R2.fastq").touch()

    result = _run_snakemake(harness, workdir, "--cores", "1", "--dry-run")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "imported_assemble_spades" in result.stdout

    listed = _run_snakemake(harness, workdir, "--list-rules")
    assert listed.returncode == 0, listed.stdout + listed.stderr
    rules = set(listed.stdout.splitlines())
    assert rules == {"imported_assemble_spades", "imported_target"}
    assert "all" not in rules


def test_standalone_rule_listing_does_not_load_nonexistent_sheet(tmp_path):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample_sheet": "missing/assembly_samples.tsv",
                "genome_context": "metagenome",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    result = _run_snakemake(
        SNAKEFILE,
        workdir,
        "--configfile",
        str(config_path),
        "--replace-workflow-config",
        "--list-rules",
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "all" in result.stdout.splitlines()


def test_standalone_missing_sample_sheet_fails_schema_validation(tmp_path):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    config_path = tmp_path / "config.yaml"
    config_path.write_text("genome_context: metagenome\n", encoding="utf-8")

    result = _run_snakemake(
        SNAKEFILE,
        workdir,
        "--configfile",
        str(config_path),
        "--replace-workflow-config",
        "--cores",
        "1",
        "--dry-run",
    )
    output = result.stdout + result.stderr

    assert result.returncode != 0
    assert "sample_sheet" in output
    assert "required" in output


def test_standalone_nonexistent_sheet_fails_when_targets_are_requested(tmp_path):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    config_path = tmp_path / "config.yaml"
    missing_sheet = tmp_path / "missing.tsv"
    config_path.write_text(
        yaml.safe_dump(
            {
                "sample_sheet": str(missing_sheet),
                "genome_context": "metagenome",
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    result = _run_snakemake(
        SNAKEFILE,
        workdir,
        "--configfile",
        str(config_path),
        "--replace-workflow-config",
        "--cores",
        "1",
        "--dry-run",
    )
    output = result.stdout + result.stderr

    assert result.returncode != 0
    assert f"Sample sheet does not exist: {missing_sheet}" in output


def test_standalone_invalid_sheet_fails_when_targets_are_requested(tmp_path):
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    sample_sheet = tmp_path / "invalid.tsv"
    sample_sheet.write_text(
        "sample\tassembly_type\tshort_read1\tshort_read2\tlong_reads\tlong_read_type\n"
        "invalid sample\tshort\treads_R1.fastq\treads_R2.fastq\t\t\n",
        encoding="utf-8",
    )
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

    result = _run_snakemake(
        SNAKEFILE,
        workdir,
        "--configfile",
        str(config_path),
        "--replace-workflow-config",
        "--cores",
        "1",
        "--dry-run",
    )
    output = result.stdout + result.stderr

    assert result.returncode != 0
    assert "invalid sample" in output
    assert "does not match" in output
