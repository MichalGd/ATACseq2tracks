#!/usr/bin/env python3
"""Combine DiffBind and DESeq2ATAC pair-level summaries into TSV and HTML."""

from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


REQUIRED = [
    "module", "peak_type", "comparison_id", "numerator", "reference",
    "numerator_replicates", "reference_replicates", "consensus_regions",
    "tested_sites", "significant_sites", "higher_in_numerator",
    "higher_in_reference", "alpha", "min_abs_log2fc", "status",
    "results_all", "results_significant", "summary_file", "message",
]


def discover(output_dir: Path) -> list[Path]:
    roots = [output_dir / "diffbind_results", output_dir / "deseq2atac"]
    return sorted(
        path
        for root in roots
        if root.exists()
        for path in root.rglob("differential_accessibility_comparisons.tsv")
    )


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = sorted(set(REQUIRED) - set(reader.fieldnames or []))
            if missing:
                raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
            for row in reader:
                normalized = {column: (row.get(column) or "NA") for column in REQUIRED}
                normalized["source_table"] = str(path)
                rows.append(normalized)
    return rows


def write_tsv(rows: list[dict[str, str]], path: Path) -> None:
    columns = REQUIRED + ["source_table"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def write_html(rows: list[dict[str, str]], path: Path) -> None:
    completed = [row for row in rows if row["status"] == "SUCCESS"]
    total_significant = sum(numeric(row["significant_sites"]) for row in completed)
    displayed = [
        "module", "peak_type", "comparison_id", "numerator", "reference",
        "numerator_replicates", "reference_replicates", "tested_sites",
        "significant_sites", "higher_in_numerator", "higher_in_reference", "status",
    ]
    lines = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8">',
        "<title>ATACseq2tracks differential accessibility summary</title>",
        "<style>body{font-family:Arial,sans-serif;margin:2rem;color:#222}",
        "table{border-collapse:collapse;width:100%;font-size:.9rem}",
        "th,td{border:1px solid #ccc;padding:.4rem;text-align:left}",
        "th{background:#eef3f8;position:sticky;top:0}",
        ".SUCCESS{background:#edf8ed}.FAILED{background:#fdeaea}",
        ".SKIPPED{background:#fff8df}code{font-size:.85em}</style></head><body>",
        "<h1>ATACseq2tracks differential accessibility summary</h1>",
        f"<p>{len(rows)} analysis rows; {len(completed)} completed; "
        f"{total_significant} significant site calls across rows. "
        "Counts are per method, peak type and contrast and must not be added as unique loci.</p>",
        "<table><thead><tr>",
        *(f"<th>{html.escape(column)}</th>" for column in displayed),
        "</tr></thead><tbody>",
    ]
    for row in rows:
        status_class = html.escape(row["status"])
        lines.append(f'<tr class="{status_class}">')
        for column in displayed:
            lines.append(f"<td>{html.escape(row[column])}</td>")
        lines.append("</tr>")
    if not rows:
        lines.append(f'<tr><td colspan="{len(displayed)}">No differential-analysis summaries found.</td></tr>')
    lines.extend([
        "</tbody></table>",
        "<p>Positive/increased counts are relative to the named numerator; decreased counts "
        "are higher in the named reference. FDR correction is performed independently for "
        "each method, peak type and contrast.</p>",
        "</body></html>",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("report_dir", type=Path)
    args = parser.parse_args()
    args.report_dir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(discover(args.output_dir))
    write_tsv(rows, args.report_dir / "differential_accessibility_summary.tsv")
    write_html(rows, args.report_dir / "differential_accessibility_summary.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
