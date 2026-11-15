"""Semantic validation for committed smoke fixtures and their generated outputs."""

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
TEST_ROOT = ROOT / ".test"
DATA_ROOT = TEST_ROOT / "data"
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")

CASES = {
    "short": {
        "sample": "short_sample",
        "assembly_type": "short",
        "sample_sheet": "samples.short.tsv",
        "short_reads": ("short_sample_R1.fastq", "short_sample_R2.fastq"),
        "long_reads": None,
        "long_read_type": "",
        "mapping": "bwamem2",
        "assembler": "assembly-spades",
    },
    "long-nanopore": {
        "sample": "long_sample",
        "assembly_type": "long",
        "sample_sheet": "samples.long-nanopore.tsv",
        "short_reads": (),
        "long_reads": "long_sample.fastq",
        "long_read_type": "nano-raw",
        "mapping": "minimap2",
        "assembler": "assembly-flye",
    },
    "hybrid-pacbio": {
        "sample": "hybrid_sample",
        "assembly_type": "hybrid",
        "sample_sheet": "samples.hybrid-pacbio.tsv",
        "short_reads": ("hybrid_sample_R1.fastq", "hybrid_sample_R2.fastq"),
        "long_reads": "hybrid_sample_long.fastq",
        "long_read_type": "pacbio-raw",
        "mapping": "bwamem2",
        "assembler": "assembly-spades",
    },
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def read_fastq(path):
    """Return identifiers from a structurally valid, uncompressed FASTQ."""
    path = Path(path)
    _require(path.is_file(), f"FASTQ does not exist: {path}")
    lines = path.read_text(encoding="ascii").splitlines()
    _require(lines, f"FASTQ is empty: {path}")
    _require(len(lines) % 4 == 0, f"FASTQ has an incomplete record: {path}")

    identifiers = []
    for offset in range(0, len(lines), 4):
        header, sequence, plus, quality = lines[offset : offset + 4]
        record = offset // 4 + 1
        _require(header.startswith("@"), f"FASTQ record {record} lacks @: {path}")
        _require(plus.startswith("+"), f"FASTQ record {record} lacks +: {path}")
        _require(sequence, f"FASTQ record {record} has no sequence: {path}")
        _require(
            set(sequence.upper()) <= DNA_ALPHABET,
            f"FASTQ record {record} contains non-DNA symbols: {path}",
        )
        _require(
            len(sequence) == len(quality),
            f"FASTQ record {record} sequence/quality lengths differ: {path}",
        )
        _require(
            all(33 <= ord(character) <= 126 for character in quality),
            f"FASTQ record {record} has invalid quality characters: {path}",
        )
        identifier = header[1:].split()[0]
        _require(identifier, f"FASTQ record {record} has no identifier: {path}")
        identifiers.append(identifier)

    _require(
        len(identifiers) == len(set(identifiers)),
        f"FASTQ contains duplicate identifiers: {path}",
    )
    return identifiers


def validate_pair(read1, read2):
    """Validate ordered /1 and /2 FASTQ mates without comparing sequences."""
    identifiers1 = read_fastq(read1)
    identifiers2 = read_fastq(read2)
    _require(
        len(identifiers1) == len(identifiers2),
        f"Paired FASTQs have different record counts: {read1}, {read2}",
    )
    for record, (identifier1, identifier2) in enumerate(
        zip(identifiers1, identifiers2, strict=True), start=1
    ):
        _require(
            identifier1.endswith("/1") and identifier2.endswith("/2"),
            f"FASTQ pair record {record} lacks /1 and /2 mate suffixes",
        )
        _require(
            identifier1[:-2] == identifier2[:-2],
            f"FASTQ pair record {record} has mismatched identifiers",
        )


def _load_single_sample(case_name):
    case = CASES[case_name]
    sheet = TEST_ROOT / "config" / case["sample_sheet"]
    with sheet.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    _require(len(rows) == 1, f"Fixture sheet must have exactly one row: {sheet}")
    return rows[0]


def validate_fixtures():
    """Validate all six committed FASTQs and their mode assignments."""
    expected_files = {
        filename
        for case in CASES.values()
        for filename in (*case["short_reads"], case["long_reads"])
        if filename
    }
    actual_files = {path.name for path in DATA_ROOT.glob("*.fastq")}
    _require(
        actual_files == expected_files,
        f"Expected exactly the six declared FASTQs; found {sorted(actual_files)}",
    )

    for case_name, case in CASES.items():
        row = _load_single_sample(case_name)
        _require(row["sample"] == case["sample"], f"Wrong sample for {case_name}")
        _require(
            row["assembly_type"].strip().lower() == case["assembly_type"],
            f"Wrong assembly type for {case_name}",
        )
        _require(
            row["long_read_type"].strip().lower() == case["long_read_type"],
            f"Wrong long-read type for {case_name}",
        )

        short_reads = case["short_reads"]
        expected_short_paths = tuple(f"data/{name}" for name in short_reads)
        observed_short_paths = (row["short_read1"], row["short_read2"])
        if short_reads:
            _require(
                observed_short_paths == expected_short_paths,
                f"Wrong paired-read assignment for {case_name}",
            )
            validate_pair(*(DATA_ROOT / name for name in short_reads))
        else:
            _require(
                observed_short_paths == ("", ""),
                f"Unexpected short reads for {case_name}",
            )

        long_reads = case["long_reads"]
        expected_long_path = f"data/{long_reads}" if long_reads else ""
        _require(
            row["long_reads"] == expected_long_path,
            f"Wrong long-read assignment for {case_name}",
        )
        if long_reads:
            read_fastq(DATA_ROOT / long_reads)

        config_path = TEST_ROOT / "config" / f"{case_name}.yaml"
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
        quast_extra = shlex.split(config.get("quast_extra", ""))
        expected_reference = (
            f"results/{case['sample']}/{case['assembler']}/assembly.fasta"
        )
        _require(
            quast_extra == ["-r", expected_reference, "--fast"],
            f"{config_path} must use its synthetic assembly as the offline QUAST reference",
        )
        if case["assembly_type"] in {"short", "hybrid"}:
            extra = shlex.split(config.get("metaspades_extra", ""))
            _require(
                extra == ["--phred-offset", "33"],
                f"{config_path} must set only --phred-offset 33 as SPAdes extra",
            )
            _require(
                config.get("metaspades_memory_gb") == 4,
                f"{config_path} must use the 4 GB SPAdes fixture limit",
            )


def read_fasta(path):
    """Return contig lengths from a nonempty DNA FASTA."""
    path = Path(path)
    _require(path.is_file(), f"FASTA does not exist: {path}")
    _require(path.stat().st_size > 0, f"FASTA is empty: {path}")

    contigs = []
    current_length = None
    with path.open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            _require(line, f"FASTA has a blank line at {line_number}: {path}")
            if line.startswith(">"):
                _require(len(line) > 1, f"FASTA has an empty header: {path}")
                if current_length is not None:
                    _require(current_length > 0, f"FASTA has an empty contig: {path}")
                    contigs.append(current_length)
                current_length = 0
            else:
                _require(
                    current_length is not None,
                    f"FASTA sequence precedes its header at line {line_number}: {path}",
                )
                _require(
                    set(line.upper()) <= DNA_ALPHABET,
                    f"FASTA has non-DNA symbols at line {line_number}: {path}",
                )
                current_length += len(line)

    _require(current_length is not None, f"FASTA has no contig header: {path}")
    _require(current_length > 0, f"FASTA has an empty final contig: {path}")
    contigs.append(current_length)
    return contigs


def _result_paths(case_name):
    case = CASES[case_name]
    sample_root = TEST_ROOT / "results" / case["sample"]
    return case, sample_root, sample_root / case["assembler"] / "assembly.fasta"


def validate_assembly(case_name):
    """Validate a case's assembler output without exact-sequence assertions."""
    _, _, assembly = _result_paths(case_name)
    return read_fasta(assembly)


def _run_samtools(*arguments):
    samtools = shutil.which("samtools")
    if samtools is None:
        candidates = sorted((TEST_ROOT / ".snakemake" / "conda").glob("*/bin/samtools"))
        _require(
            candidates,
            "samtools is unavailable on PATH and in generated Snakemake environments",
        )
        samtools = str(candidates[0])
    result = subprocess.run(
        [samtools, *map(str, arguments)],
        check=False,
        capture_output=True,
        text=True,
    )
    _require(
        result.returncode == 0,
        f"samtools {' '.join(map(str, arguments))} failed:\n{result.stderr}",
    )
    return result.stdout


def validate_complete(case_name):
    """Validate a complete local route, mapping index, and logs."""
    case, sample_root, assembly = _result_paths(case_name)
    contigs = validate_assembly(case_name)

    report = sample_root / "qc-quast" / "report.txt"
    _require(report.is_file() and report.stat().st_size > 0, f"Missing QUAST: {report}")

    mapping = case["mapping"]
    mapping_root = sample_root / f"mapping-{mapping}"
    bam = mapping_root / "alignment.bam"
    bai = mapping_root / "alignment.bam.bai"
    _require(bam.is_file() and bam.stat().st_size > 0, f"Missing BAM: {bam}")
    _require(bai.is_file() and bai.stat().st_size > 0, f"Missing BAI: {bai}")
    opposite_mapping = "minimap2" if mapping == "bwamem2" else "bwamem2"
    _require(
        not (sample_root / f"mapping-{opposite_mapping}").exists(),
        f"Unexpected {opposite_mapping} route for {case_name}",
    )
    opposite_assembler = (
        "assembly-flye" if case["assembler"] == "assembly-spades" else "assembly-spades"
    )
    _require(
        not (sample_root / opposite_assembler).exists(),
        f"Unexpected {opposite_assembler} route for {case_name}",
    )

    expected_root_entries = {
        case["assembler"],
        f"mapping-{mapping}",
        "qc-quast",
    }
    actual_root_entries = {path.name for path in sample_root.iterdir()}
    _require(
        actual_root_entries == expected_root_entries,
        f"Unexpected result tree for {case_name}: "
        f"{sorted(actual_root_entries)} != {sorted(expected_root_entries)}",
    )
    for obsolete in ("assembly", "mapping", "quast"):
        _require(
            not (sample_root / obsolete).exists(),
            f"Obsolete output path remains for {case_name}: {sample_root / obsolete}",
        )

    _run_samtools("quickcheck", "-v", bam)
    header = _run_samtools("view", "-H", bam)
    hd_lines = [line for line in header.splitlines() if line.startswith("@HD\t")]
    _require(len(hd_lines) == 1, f"BAM must have one @HD header: {bam}")
    _require(
        "SO:coordinate" in hd_lines[0].split("\t"),
        f"BAM is not coordinate sorted: {bam}",
    )
    _require(_run_samtools("idxstats", bam).strip(), f"BAI is not usable: {bai}")

    expected_mapping_entries = {
        "alignment.bam",
        "alignment.bam.bai",
        "map.log",
        "index_bam.log",
    }
    if mapping == "bwamem2":
        expected_mapping_entries.update(
            {
                "assembly.0123",
                "assembly.amb",
                "assembly.ann",
                "assembly.bwt.2bit.64",
                "assembly.pac",
                "index.log",
            }
        )
    actual_mapping_entries = {path.name for path in mapping_root.iterdir()}
    _require(
        actual_mapping_entries == expected_mapping_entries,
        f"Unexpected mapping tree for {case_name}: "
        f"{sorted(actual_mapping_entries)} != {sorted(expected_mapping_entries)}",
    )

    logs = [
        sample_root / case["assembler"] / "run.log",
        sample_root / "qc-quast" / "run.log",
        mapping_root / "map.log",
        mapping_root / "index_bam.log",
    ]
    if mapping == "bwamem2":
        logs.append(mapping_root / "index.log")
    missing_logs = [str(path) for path in logs if not path.is_file()]
    _require(not missing_logs, f"Missing workflow logs: {missing_logs}")

    native_log = (
        sample_root
        / case["assembler"]
        / ("spades.log" if case["assembler"] == "assembly-spades" else "flye.log")
    )
    _require(
        native_log.is_file() and native_log.stat().st_size > 0,
        f"Missing native assembler log: {native_log}",
    )
    return contigs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scope", choices=("fixtures", "assembly", "complete"))
    parser.add_argument("--case", choices=CASES)
    arguments = parser.parse_args()

    if arguments.scope == "fixtures":
        _require(arguments.case is None, "--case is not used with fixture validation")
        validate_fixtures()
        print("Validated all six synthetic FASTQ fixtures.")
        return

    _require(arguments.case is not None, f"--case is required for {arguments.scope}")
    if arguments.scope == "assembly":
        contigs = validate_assembly(arguments.case)
        print(f"Validated {arguments.case} assembly: {len(contigs)} contig(s).")
    else:
        contigs = validate_complete(arguments.case)
        print(f"Validated {arguments.case} complete results: {len(contigs)} contig(s).")


if __name__ == "__main__":
    main()
