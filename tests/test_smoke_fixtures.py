from pathlib import Path

from smoke_validation import DATA_ROOT, validate_fixtures


def test_synthetic_fastq_fixtures_are_valid():
    validate_fixtures()


def test_fixture_documentation_exists():
    assert (Path(DATA_ROOT) / "README.md").is_file()
