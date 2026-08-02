#!/usr/bin/env python3
"""Create a strand-aware BED6 TSS reference from a GTF annotation."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if len(parts) == 2:
            result[parts[0]] = parts[1].strip().strip('"')
    return result


def chrom_key(chrom: str):
    core = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    if core.isdigit():
        return (0, int(core), chrom)
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if core.upper() in order:
        return (0, order[core.upper()], chrom)
    return (1, 0, chrom)


def read_features(path: Path, wanted: str) -> dict[tuple[str, int, str], str]:
    records: dict[tuple[str, int, str], str] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != wanted:
                continue
            chrom, start_text, end_text, strand = fields[0], fields[3], fields[4], fields[6]
            if strand not in {"+", "-"}:
                continue
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(f"invalid GTF coordinates at line {line_number}") from exc
            tss0 = start - 1 if strand == "+" else end - 1
            attrs = attributes(fields[8])
            name = (
                attrs.get("gene_name")
                or attrs.get("transcript_id")
                or attrs.get("gene_id")
                or f"TSS_{chrom}_{tss0 + 1}_{strand}"
            )
            records.setdefault((chrom, tss0, strand), name)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gtf", type=Path)
    parser.add_argument("output_bed", type=Path)
    args = parser.parse_args()
    if not args.gtf.is_file():
        parser.error(f"GTF not found: {args.gtf}")

    records = read_features(args.gtf, "transcript")
    source_feature = "transcript"
    if not records:
        records = read_features(args.gtf, "gene")
        source_feature = "gene"
    if not records:
        raise SystemExit("ERROR: no strand-aware transcript or gene records found in GTF")

    args.output_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.output_bed.open("w", encoding="utf-8", newline="\n") as out:
        for (chrom, tss0, strand), name in sorted(
            records.items(), key=lambda item: (chrom_key(item[0][0]), item[0][1], item[0][2])
        ):
            safe_name = re.sub(r"[\t\r\n ]+", "_", name)
            out.write(f"{chrom}\t{tss0}\t{tss0 + 1}\t{safe_name}\t0\t{strand}\n")

    print(f"Wrote {len(records)} unique TSS records from GTF {source_feature} features: {args.output_bed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
