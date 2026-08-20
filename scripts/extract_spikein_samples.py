#!/usr/bin/env python3
"""Emit one tab-separated row per dm6-declared biological sample key."""

from __future__ import annotations

import argparse
import csv
import math
import sys


REQUIRED = {
    "sample_id",
    "replicate",
    "layout",
    "genome",
    "is_control",
    "blacklist",
    "spikein_genome",
    "spikein_to_host_ratio",
    "spikein_stage",
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("samplesheet")
    args = parser.parse_args()

    with open(args.samplesheet, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(REQUIRED.difference(reader.fieldnames or []))
        if missing:
            raise SystemExit("Spike-in mode requires samplesheet columns: " + ", ".join(missing))
        rows = list(reader)

    seen: dict[tuple[str, str], tuple[str, ...]] = {}
    emitted = 0
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ["key", "sample_id", "replicate", "layout", "host_genome", "host_blacklist",
         "spikein_genome", "spikein_to_host_ratio", "spikein_stage", "is_control"]
    )
    for row_number, row in enumerate(rows, start=2):
        sample_id = row["sample_id"].strip()
        replicate = row["replicate"].strip()
        key = (sample_id, replicate)
        values = (
            row["layout"].strip().upper(),
            row["genome"].strip(),
            row["blacklist"].strip(),
            row["spikein_genome"].strip().lower(),
            row["spikein_to_host_ratio"].strip(),
            row["spikein_stage"].strip().lower(),
            row["is_control"].strip().lower(),
        )
        if key in seen:
            if values != seen[key]:
                raise SystemExit(f"Inconsistent spike-in metadata for {sample_id}_bioR{replicate} at row {row_number}")
            continue
        seen[key] = values
        layout, host_genome, host_blacklist, spike_genome, ratio_text, stage, is_control = values
        if not spike_genome:
            continue
        if spike_genome != "dm6":
            raise SystemExit(f"Unsupported spike-in genome for {sample_id}_bioR{replicate}: {spike_genome}")
        try:
            ratio = float(ratio_text)
        except ValueError as exc:
            raise SystemExit(f"Invalid spike-in ratio for {sample_id}_bioR{replicate}: {ratio_text}") from exc
        if not math.isfinite(ratio) or ratio <= 0:
            raise SystemExit(f"Spike-in ratio must be positive for {sample_id}_bioR{replicate}")
        if stage != "pre_tagmentation_nuclei":
            raise SystemExit(
                f"Spike-in calibration requires pre_tagmentation_nuclei for {sample_id}_bioR{replicate}"
            )
        writer.writerow(
            [f"{sample_id}_bioR{replicate}", sample_id, replicate, layout, host_genome,
             host_blacklist, spike_genome, format(ratio, ".12g"), stage, is_control]
        )
        emitted += 1

    if emitted == 0:
        raise SystemExit("Spike-in mode is enabled but no sample declares spikein_genome=dm6")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
