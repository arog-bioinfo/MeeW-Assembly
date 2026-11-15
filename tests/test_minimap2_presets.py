import shlex
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from minimap2_presets import compose_minimap2_extra, minimap2_preset


@pytest.mark.parametrize(
    ("long_read_type", "expected_preset"),
    [
        ("nano-raw", "map-ont"),
        ("nano-corr", "map-ont"),
        ("nano-hq", "map-ont"),
        ("pacbio-raw", "map-pb"),
        ("pacbio-corr", "map-pb"),
        ("pacbio-hifi", "map-hifi"),
    ],
)
def test_exact_preset_mapping_and_composition(long_read_type, expected_preset):
    assert minimap2_preset(long_read_type) == expected_preset
    arguments = shlex.split(compose_minimap2_extra(long_read_type, "-N 2"))

    assert arguments[:2] == ["-x", expected_preset]
    assert arguments.count("-x") == 1
    assert arguments[2:] == ["-N", "2"]


def test_valid_extra_arguments_are_preserved_after_preset():
    extra = '--secondary=no -R "@RG\\tID:sample-x"'

    composed = compose_minimap2_extra("pacbio-corr", extra)

    assert composed == f"-x map-pb {extra}"
    assert shlex.split(composed) == [
        "-x",
        "map-pb",
        "--secondary=no",
        "-R",
        "@RG\\tID:sample-x",
    ]


@pytest.mark.parametrize("extra", ["", " "])
def test_empty_or_whitespace_extra_adds_only_the_preset(extra):
    assert compose_minimap2_extra("nano-raw", extra) == "-x map-ont"


def test_non_option_value_containing_dash_x_is_not_rejected():
    extra = "--junc-bed sample-x-regions.bed"

    assert (
        compose_minimap2_extra("pacbio-hifi", extra)
        == "-x map-hifi --junc-bed sample-x-regions.bed"
    )


def test_unsupported_long_read_type_is_rejected():
    with pytest.raises(ValueError, match="Unsupported `long_read_type`.*pacbio"):
        compose_minimap2_extra("pacbio", "")


def test_malformed_shell_quoting_is_rejected():
    with pytest.raises(ValueError, match="valid POSIX shell argument string"):
        compose_minimap2_extra("nano-raw", '--secondary="no')


@pytest.mark.parametrize(
    ("extra", "error_form"),
    [
        ("-x map-pb", "separated preset option"),
        ("-x", "bare `-x` with no value"),
        ("-xmap-pb", "attached preset option"),
        ("-axmap-pb", "bundled preset option"),
    ],
)
def test_user_supplied_preset_forms_are_rejected(extra, error_form):
    with pytest.raises(ValueError, match=error_form):
        compose_minimap2_extra("nano-raw", extra)
