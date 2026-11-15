import json
import subprocess
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflow" / "Snakefile"
SCHEMA = ROOT / "workflow" / "schemas" / "config.schema.yaml"
TEST_DIRECTORY = ROOT / ".test"
SHORT_SAMPLE_SHEET = TEST_DIRECTORY / "config" / "samples.short.tsv"
SORT_CONFIG_VALIDATOR = Draft202012Validator(
    yaml.safe_load(SCHEMA.read_text(encoding="utf-8"))
)


def _config(**overrides):
    config = {
        "sample_sheet": str(SHORT_SAMPLE_SHEET),
        "genome_context": "metagenome",
    }
    config.update(overrides)
    return config


def _run_snakemake(tmp_path, config):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--directory",
            str(TEST_DIRECTORY),
            "--configfile",
            str(config_path),
            "--replace-workflow-config",
            "--nolock",
            "--d3dag",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("sort_engine", ["samtools", "picard"])
def test_coordinate_sorting_engines_are_accepted_by_schema(sort_engine):
    SORT_CONFIG_VALIDATOR.validate(
        _config(bwamem2_sort=sort_engine, bwamem2_sort_order="coordinate")
    )


def test_omitted_sort_settings_construct_bwamem2_route(tmp_path):
    result = _run_snakemake(tmp_path, _config())
    assert result.returncode == 0, result.stdout + result.stderr
    graph = json.loads(result.stdout)
    rules = {node["value"]["rule"] for node in graph["nodes"]}
    assert {"bwamem2_map", "bwamem2_map_index"} <= rules


@pytest.mark.parametrize(
    ("field", "invalid_value", "accepted_values"),
    [
        ("bwamem2_sort", "none", ["samtools", "picard"]),
        ("bwamem2_sort_order", "queryname", ["coordinate"]),
    ],
)
def test_unsafe_sort_settings_are_rejected_by_schema(
    field, invalid_value, accepted_values
):
    errors = list(SORT_CONFIG_VALIDATOR.iter_errors(_config(**{field: invalid_value})))

    assert len(errors) == 1
    error = errors[0]
    assert list(error.absolute_path) == [field]
    assert error.validator == "enum"
    assert error.instance == invalid_value
    assert error.validator_value == accepted_values
