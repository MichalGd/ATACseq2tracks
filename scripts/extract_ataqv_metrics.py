#!/usr/bin/env python3
"""Extract selected scalar TSS/fragment metrics from an ataqv JSON file."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path


TOKENS = (
    "tss_enrichment",
    "short_mononucleosomal_ratio",
    "fragment_length_distance",
    "fragment_length_mean",
    "median_fragment_length",
)


def open_json(path: Path):
    with path.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8") if compressed else path.open(encoding="utf-8")


def normalized(path: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", path.lower()).strip("_")


def scalars(value, path=""):
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else str(key)
            yield from scalars(child, child_path)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            if isinstance(child, dict):
                yield from scalars(child, f"{path}[{index}]")
    elif value is None or isinstance(value, (str, int, float, bool)):
        yield path, value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics_json", type=Path)
    parser.add_argument("sample_id")
    parser.add_argument("output_tsv", type=Path)
    args = parser.parse_args()

    with open_json(args.metrics_json) as handle:
        document = json.load(handle)
    selected = []
    for metric_path, value in scalars(document):
        if any(token in normalized(metric_path) for token in TOKENS):
            selected.append((metric_path, value))

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", encoding="utf-8", newline="\n") as out:
        out.write("sample_id\tmetric_path\tvalue\n")
        for metric_path, value in selected:
            out.write(f"{args.sample_id}\t{metric_path}\t{value}\n")
    if not selected:
        print(f"WARNING: no selected scalar metrics found in {args.metrics_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
