#!/usr/bin/env python3
"""Calculate compact paired-end ATAC fragment metrics and plot periodicity."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REGIONS = {
    "nucleosome_free": (1, 100),
    "mononucleosome": (180, 247),
    "dinucleosome": (315, 473),
    "trinucleosome": (558, 615),
}


def collect_lengths(bam: Path, maximum: int) -> tuple[list[int], int, int]:
    histogram = [0] * (maximum + 1)
    total = 0
    overflow = 0
    # One record per proper pair (read 1); reject unmapped, secondary,
    # QC-failed, duplicate and supplementary records.
    command = ["samtools", "view", "-f", "66", "-F", "3852", str(bam)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        try:
            length = abs(int(fields[8]))
        except ValueError:
            continue
        if length == 0:
            continue
        total += 1
        if length <= maximum:
            histogram[length] += 1
        else:
            overflow += 1
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"samtools view failed with exit code {return_code}: {bam}")
    return histogram, total, overflow


def region_count(histogram: list[int], bounds: tuple[int, int]) -> int:
    lower, upper = bounds
    upper = min(upper, len(histogram) - 1)
    return sum(histogram[lower : upper + 1]) if lower <= upper else 0


def peak_length(histogram: list[int], lower: int, upper: int) -> int | None:
    upper = min(upper, len(histogram) - 1)
    if lower > upper or not any(histogram[lower : upper + 1]):
        return None
    return max(range(lower, upper + 1), key=lambda value: histogram[value])


def ratio(numerator: int, denominator: int) -> str:
    return "NA" if denominator == 0 else f"{numerator / denominator:.6f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--max-fragment-length", type=int, default=1000)
    args = parser.parse_args()
    if not args.bam.is_file():
        parser.error(f"BAM not found: {args.bam}")
    if args.max_fragment_length < 615:
        parser.error("--max-fragment-length must be at least 615 bp")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    histogram, total, overflow = collect_lengths(args.bam, args.max_fragment_length)
    if total == 0:
        raise SystemExit(f"ERROR: no valid paired-end fragments found in {args.bam}")

    counts = {name: region_count(histogram, bounds) for name, bounds in REGIONS.items()}
    mono_peak = peak_length(histogram, 150, 260)
    di_peak = peak_length(histogram, 300, 500)
    spacing = di_peak - mono_peak if mono_peak is not None and di_peak is not None else None
    periodic_count = counts["mononucleosome"] + counts["dinucleosome"] + counts["trinucleosome"]

    prefix = args.out_dir / args.sample
    histogram_path = Path(f"{prefix}.fragment_length_histogram.tsv")
    metrics_path = Path(f"{prefix}.nucleosome_periodicity_metrics.tsv")
    with histogram_path.open("w", encoding="utf-8", newline="\n") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["fragment_length_bp", "fragment_count"])
        writer.writerows((length, histogram[length]) for length in range(1, len(histogram)))

    metrics = [
        ("proper_pair_fragments", str(total)),
        ("fragments_above_plot_max", str(overflow)),
        ("nucleosome_free_1_100_count", str(counts["nucleosome_free"])),
        ("nucleosome_free_1_100_fraction", ratio(counts["nucleosome_free"], total)),
        ("mononucleosome_180_247_count", str(counts["mononucleosome"])),
        ("mononucleosome_180_247_fraction", ratio(counts["mononucleosome"], total)),
        ("dinucleosome_315_473_count", str(counts["dinucleosome"])),
        ("dinucleosome_315_473_fraction", ratio(counts["dinucleosome"], total)),
        ("trinucleosome_558_615_count", str(counts["trinucleosome"])),
        ("trinucleosome_558_615_fraction", ratio(counts["trinucleosome"], total)),
        ("nfr_to_mononucleosome_ratio", ratio(counts["nucleosome_free"], counts["mononucleosome"])),
        ("periodic_fragment_fraction", ratio(periodic_count, total)),
        ("mononucleosome_local_peak_bp", str(mono_peak) if mono_peak is not None else "NA"),
        ("dinucleosome_local_peak_bp", str(di_peak) if di_peak is not None else "NA"),
        ("mono_to_di_peak_spacing_bp", str(spacing) if spacing is not None else "NA"),
        ("mono_to_di_spacing_deviation_from_200_bp", str(abs(spacing - 200)) if spacing is not None else "NA"),
    ]
    with metrics_path.open("w", encoding="utf-8", newline="\n") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "metric", "value"])
        writer.writerows((args.sample, name, value) for name, value in metrics)

    x = list(range(1, len(histogram)))
    y = histogram[1:]
    scale = max(total, 1) / 1_000_000
    y_rpm = [value / scale for value in y]
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    colors = {
        "nucleosome_free": "#4C78A8",
        "mononucleosome": "#F58518",
        "dinucleosome": "#54A24B",
        "trinucleosome": "#E45756",
    }
    for name, (lower, upper) in REGIONS.items():
        ax.axvspan(lower, upper, color=colors[name], alpha=0.12, linewidth=0)
    ax.plot(x, y_rpm, color="#263238", linewidth=1.15)
    for position in (200, 400, 600):
        ax.axvline(position, color="#777777", linestyle="--", linewidth=0.7, alpha=0.55)
    ax.set(xlabel="Paired-end fragment length (bp)", ylabel="Fragments per million", title=args.sample)
    ax.set_xlim(0, args.max_fragment_length)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{prefix}.fragment_length_periodicity.png", dpi=300)
    fig.savefig(f"{prefix}.fragment_length_periodicity.pdf")
    plt.close(fig)
    print(f"Wrote fragment periodicity metrics and plots for {args.sample}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
