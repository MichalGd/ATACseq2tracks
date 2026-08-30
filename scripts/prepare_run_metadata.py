#!/usr/bin/env python3
"""Write persistent sequencing-unit, technical-merge, plan and resource tables."""

from __future__ import annotations

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


STAGES = [
    ("1", "fastqc_raw"), ("2", "trim"), ("3", "fastqc_trimmed"),
    ("4", "alignment"), ("5", "dedup"), ("6", "filtering"),
    ("6s", "spikein"), ("7", "coverage"), ("8", "merged_tracks"),
    ("9", "peaks"), ("10", "qc"), ("10b", "filtering_sensitivity"),
    ("11", "diffbind_prep"), ("12", "diffbind"),
    ("12a", "deseq2atac"), ("13", "browser"), ("14", "report"),
]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def integer(values: dict[str, str], key: str, default: int) -> int:
    try:
        return int(values.get(key, str(default)))
    except ValueError:
        return default


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--resolved-config-tsv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.samplesheet.open(encoding="utf-8-sig", newline="") as handle:
        units = list(csv.DictReader(handle))
    if not units:
        raise SystemExit("ERROR: samplesheet contains no data rows")

    unit_fields = list(units[0]) + ["biological_library_key"]
    unit_rows = []
    groups: OrderedDict[tuple[str, str], list[dict[str, str]]] = OrderedDict()
    for row in units:
        key = (row["sample_id"], row["replicate"])
        groups.setdefault(key, []).append(row)
        unit_rows.append({**row, "biological_library_key": f"{key[0]}_bioR{key[1]}"})
    write_tsv(args.output_dir / "validated_sequencing_units.tsv", unit_rows, unit_fields)

    libraries = []
    audits = []
    for (sample_id, replicate), rows in groups.items():
        first = rows[0]
        key = f"{sample_id}_bioR{replicate}"
        libraries.append({
            "biological_library_key": key, "sample_id": sample_id,
            "replicate": replicate, "condition": first.get("condition", ""),
            "genome": first.get("genome", ""), "layout": first.get("layout", ""),
            "technical_unit_count": len(rows),
        })
        audits.append({
            "biological_library_key": key,
            "technical_replicates": ",".join(row.get("tech_replicate", "") for row in rows),
            "technical_unit_count": len(rows),
            "ordered_fastq_1": ";".join(row.get("fastq_1", "") for row in rows),
            "ordered_fastq_2": ";".join(row.get("fastq_2", "") for row in rows),
            "merge_stage": "before_trim",
            "downstream_library_count": 1,
        })
    write_tsv(args.output_dir / "biological_libraries.tsv", libraries, list(libraries[0]))
    write_tsv(args.output_dir / "technical_merge_audit.tsv", audits, list(audits[0]))
    write_tsv(
        args.output_dir / "planned_stages.tsv",
        [{"step": step, "stage": stage, "default_order": index} for index, (step, stage) in enumerate(STAGES, 1)],
        ["step", "stage", "default_order"],
    )

    with args.resolved_config_tsv.open(encoding="utf-8", newline="") as handle:
        config = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    resources = [
        ("alignment", "THREADS_PARALLEL_JOBS", "THREADS_ALIGN", 8, 16),
        ("sample_qc", "QC_SAMPLE_PARALLEL_JOBS", "THREADS_DEEPTOOLS", 4, 16),
        ("ataqv", "ATAQV_PARALLEL_JOBS", "THREADS_ATAQV", 4, 8),
        ("coverage", "TRACK_PARALLEL_JOBS", "THREADS_BIGWIG", 2, 16),
        ("spikein", "SPIKEIN_PARALLEL_JOBS", "THREADS_ALIGN", 2, 16),
    ]
    resource_rows = []
    budget = integer(config, "TOTAL_CPU_BUDGET", 0)
    for stage, jobs_key, threads_key, jobs_default, threads_default in resources:
        jobs = integer(config, jobs_key, jobs_default)
        threads = integer(config, threads_key, threads_default)
        demand = jobs * threads
        resource_rows.append({
            "stage": stage, "jobs_setting": jobs_key, "jobs": jobs,
            "threads_setting": threads_key, "threads_per_job": threads,
            "maximum_requested_threads": demand, "configured_cpu_budget": budget,
            "exceeds_budget": "yes" if budget and demand > budget else "no",
        })
    write_tsv(args.output_dir / "resource_budget.tsv", resource_rows, list(resource_rows[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
