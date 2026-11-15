from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_enabled_ci_workflow_matches_template_plus_contract():
    workflow = yaml.load(
        (ROOT / ".github" / "workflows" / "main.yaml").read_text(encoding="utf-8"),
        Loader=yaml.BaseLoader,
    )

    assert {"push", "pull_request", "workflow_dispatch"} <= workflow["on"].keys()
    assert workflow["on"]["push"]["branches"] == ["main"]
    assert workflow["on"]["pull_request"]["branches"] == ["main"]
    assert workflow["permissions"] == {"contents": "read"}

    jobs = workflow["jobs"]
    assert set(jobs) == {"formatting", "linting", "python-tests", "testing"}

    testing = jobs["testing"]
    assert set(testing["needs"]) == {"formatting", "linting", "python-tests"}
    suites = {
        (suite["id"], suite["config"])
        for suite in testing["strategy"]["matrix"]["suite"]
    }
    assert suites == {
        ("short", "short.yaml"),
        ("long-nanopore", "long-nanopore.yaml"),
        ("hybrid-pacbio", "hybrid-pacbio.yaml"),
    }

    fixture_step = next(
        step for step in testing["steps"] if step.get("name") == "Run fixture workflow"
    )
    assert fixture_step["uses"] == "snakemake/snakemake-github-action@v2"
    fixture_options = fixture_step["with"]
    assert fixture_options["snakemake-version"] == "9.19.0"
    assert fixture_options["args"] == (
        "--configfile .test/config/${{ matrix.suite.config }} "
        "--cores 2 "
        "--software-deployment-method conda "
        "--show-failed-logs "
        "--conda-cleanup-pkgs cache"
    )
    assert "assembly.fasta" not in fixture_options["args"]

    validation_step = next(
        step
        for step in testing["steps"]
        if step.get("name") == "Assert complete route outputs"
    )
    assert validation_step["run"] == (
        "python tests/smoke_validation.py complete --case ${{ matrix.suite.id }}"
    )

    report_step = next(
        step
        for step in testing["steps"]
        if step.get("name") == "Generate Snakemake report"
    )
    assert report_step["uses"] == "snakemake/snakemake-github-action@v2"
    report_options = report_step["with"]
    assert report_options["snakemake-version"] == "9.19.0"
    assert report_options["args"] == (
        "--report snakemake-report-${{ matrix.suite.id }}.zip "
        "--configfile .test/config/${{ matrix.suite.config }} "
        "--software-deployment-method conda"
    )

    artifact_step = next(
        step
        for step in testing["steps"]
        if step.get("name") == "Upload Snakemake report"
    )
    assert artifact_step["uses"] == "actions/upload-artifact@v4"
    assert artifact_step["with"] == {
        "name": "snakemake-report-${{ matrix.suite.id }}",
        "path": ".test/snakemake-report-${{ matrix.suite.id }}.zip",
        "if-no-files-found": "error",
        "retention-days": "7",
        "include-hidden-files": "true",
    }
