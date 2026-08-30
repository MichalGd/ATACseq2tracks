#!/usr/bin/env python3
"""
ATACseq2tracks v4.3.0 - sample sheet validator
Usage: python3 scripts/validate_samplesheet.py [--check-files] samplesheet.csv

Schema changes:
  v3.1.0: chipqc_annotation column removed (17-column schema).
          Old 18-column samplesheets (with chipqc_annotation) are still accepted
          and validated correctly — the column is ignored by the pipeline.
  v3.0.4: Duplicate sample_id is ALLOWED when tech_replicate values differ.
          The enforced unique key is (sample_id, replicate, tech_replicate).
          Within a tech-replicate group (same sample_id + replicate) all metadata
          columns except fastq_1, fastq_2, and tech_replicate must be identical.
          True duplicate rows (all fields identical) still fail.
          Duplicate composite keys still fail.
          control_id cross-reference uses the unique_sids set (unchanged semantics).
"""
import csv, sys, os, argparse, math

# Current 17-column schema (v3.1.0+, unchanged in v4.1.0)
REQUIRED_COLS_V31 = [
    "sample_id","fastq_1","fastq_2","layout","genome",
    "assay","factor","condition","treatment","cell_type",
    "replicate","tech_replicate","is_control","control_id",
    "macs2_mode","blacklist","output_prefix"
]

# Legacy 18-column schema (v3.0.x) — chipqc_annotation present but ignored by pipeline
REQUIRED_COLS_V30 = [
    "sample_id","fastq_1","fastq_2","layout","genome",
    "assay","factor","condition","treatment","cell_type",
    "replicate","tech_replicate","is_control","control_id",
    "macs2_mode","blacklist","chipqc_annotation","output_prefix"
]

# Optional v4.2.0 Drosophila spike-in columns. Existing 17/18-column sheets
# remain valid when the spike-in track family is disabled.
SPIKEIN_COLS = ["spikein_genome", "spikein_to_host_ratio", "spikein_stage"]

# Columns that must be consistent within a (sample_id, replicate) tech-rep group
CONSISTENT_COLS_V31 = [
    "layout","genome","assay","factor","condition","treatment",
    "cell_type","is_control","control_id","macs2_mode","blacklist",
    "output_prefix"
]
CONSISTENT_COLS_V30 = CONSISTENT_COLS_V31 + ["chipqc_annotation"]

VALID_LAYOUTS = {"PE","SE"}
VALID_GENOMES = {"hg38","mm39"}
VALID_MACS2   = {"narrow","broad","both","none"}
VALID_ISCTR   = {"TRUE","FALSE","true","false","1","0","yes","no"}
VALID_SPIKEIN_GENOMES = {"", "dm6"}
VALID_SPIKEIN_STAGES = {"", "pre_tagmentation_nuclei"}

def err(row, msg):  print(f"  [ERROR] row {row}: {msg}", file=sys.stderr); return 1
def warn(row, msg): print(f"  [WARN]  row {row}: {msg}")


def detect_schema(fieldnames):
    """Return (REQUIRED_COLS, CONSISTENT_COLS, schema_version_str)."""
    if "chipqc_annotation" in fieldnames:
        return REQUIRED_COLS_V30, CONSISTENT_COLS_V30, "v3.0.x (18-column, chipqc_annotation present — accepted, ignored by pipeline)"
    suffix = "; dm6 spike-in columns present" if all(c in fieldnames for c in SPIKEIN_COLS) else ""
    return REQUIRED_COLS_V31, CONSISTENT_COLS_V31, f"v4.3.0 (17-column core; compatible with v3.1.0{suffix})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("samplesheet")
    ap.add_argument("--check-files", action="store_true",
                    help="Verify FASTQ/BED files exist on disk")
    args = ap.parse_args()
    errors = 0
    print(f"Validating: {args.samplesheet}")

    with open(args.samplesheet, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []

        required_cols, consistent_cols, schema_ver = detect_schema(fieldnames)
        present_spikein_cols = [c for c in SPIKEIN_COLS if c in fieldnames]
        if present_spikein_cols and len(present_spikein_cols) != len(SPIKEIN_COLS):
            missing_spikein_cols = [c for c in SPIKEIN_COLS if c not in fieldnames]
            print(f"[FATAL] Partial spike-in schema: add missing columns {missing_spikein_cols}", file=sys.stderr)
            sys.exit(1)
        if len(present_spikein_cols) == len(SPIKEIN_COLS):
            consistent_cols = consistent_cols + SPIKEIN_COLS
        print(f"Schema detected      : {schema_ver}")

        missing_cols = [c for c in required_cols if c not in fieldnames]
        if missing_cols:
            print(f"[FATAL] Missing columns: {missing_cols}", file=sys.stderr)
            sys.exit(1)

        rows = list(reader)

    composite_keys       = {}
    seen_raw_rows        = {}
    unique_sids          = set()
    control_ids_declared = set()
    ip_to_control        = {}
    group_meta           = {}
    run_layouts          = set()

    # File-check columns: always check fastq_1, fastq_2, blacklist
    # chipqc_annotation is intentionally excluded even in 18-col sheets —
    # the file is no longer required by the pipeline
    FILE_CHECK_COLS = ["fastq_1", "fastq_2", "blacklist"]

    for i, row in enumerate(rows, start=2):
        sid      = row["sample_id"].strip()
        rep      = row["replicate"].strip()
        tech_rep = row["tech_replicate"].strip()

        # True duplicate row check
        raw_key = frozenset((k, v.strip()) for k, v in row.items())
        if raw_key in seen_raw_rows:
            errors += err(i, f"Exact duplicate row (identical to row {seen_raw_rows[raw_key]})")
        else:
            seen_raw_rows[raw_key] = i

        # Composite key uniqueness
        comp = (sid, rep, tech_rep)
        if comp in composite_keys:
            errors += err(i, f"Duplicate composite key (sample_id='{sid}', "
                             f"replicate='{rep}', tech_replicate='{tech_rep}') "
                             f"— first seen at row {composite_keys[comp]}")
        else:
            composite_keys[comp] = i

        unique_sids.add(sid)

        layout = row["layout"].strip().upper()
        if layout not in VALID_LAYOUTS:
            errors += err(i, f"Invalid layout '{layout}' (PE or SE)")
        else:
            run_layouts.add(layout)
        if not row["fastq_1"].strip():
            errors += err(i, "fastq_1 is empty")
        if layout == "PE" and not row["fastq_2"].strip():
            errors += err(i, "PE sample missing fastq_2")
        if layout == "SE" and row["fastq_2"].strip():
            warn(i, "SE sample has fastq_2 — will be ignored")
        if row["genome"].strip() not in VALID_GENOMES:
            errors += err(i, f"Invalid genome '{row['genome']}'")
        if row["macs2_mode"].strip().lower() not in VALID_MACS2:
            errors += err(i, f"Invalid macs2_mode '{row['macs2_mode']}'")

        if all(c in row for c in SPIKEIN_COLS):
            spike_genome = row["spikein_genome"].strip().lower()
            spike_ratio = row["spikein_to_host_ratio"].strip()
            spike_stage = row["spikein_stage"].strip().lower()
            if spike_genome not in VALID_SPIKEIN_GENOMES:
                errors += err(i, f"Invalid spikein_genome '{spike_genome}' (blank or dm6)")
            if spike_stage not in VALID_SPIKEIN_STAGES:
                errors += err(i, "spikein_stage must be blank or pre_tagmentation_nuclei")
            if spike_genome == "dm6":
                try:
                    ratio_value = float(spike_ratio)
                    if not math.isfinite(ratio_value) or ratio_value <= 0:
                        raise ValueError
                except ValueError:
                    errors += err(i, "dm6 spike-in requires a positive spikein_to_host_ratio")
                if spike_stage != "pre_tagmentation_nuclei":
                    errors += err(i, "dm6 spike-in requires spikein_stage=pre_tagmentation_nuclei")
            elif spike_ratio or spike_stage:
                errors += err(i, "spike-in ratio/stage supplied while spikein_genome is blank")

        is_ctr = row["is_control"].strip()
        if is_ctr not in VALID_ISCTR:
            errors += err(i, f"Invalid is_control '{is_ctr}'")
        is_ctr_bool = is_ctr.lower() in {"true","1","yes"}
        if is_ctr_bool:
            control_ids_declared.add(sid)
            if row["control_id"].strip():
                warn(i, "Control sample should have empty control_id")
            if row["macs2_mode"].strip().lower() != "none":
                warn(i, "Control sample: macs2_mode should be 'none'")
        else:
            ctrl = row["control_id"].strip()
            if not ctrl:
                if not row["assay"].strip().lower().startswith("atac"):
                    warn(i, "Non-ATAC sample has no control_id; MACS3 will run without control")
            else:
                ip_to_control[sid] = ctrl

        try:
            int(rep)
        except ValueError:
            errors += err(i, f"replicate must be integer, got '{rep}'")

        try:
            int(tech_rep)
        except ValueError:
            errors += err(i, f"tech_replicate must be integer, got '{tech_rep}'")

        if args.check_files:
            for col in FILE_CHECK_COLS:
                p = row.get(col, "").strip()
                if p and not os.path.exists(p):
                    errors += err(i, f"File not found: {col}={p}")

        # Metadata consistency within tech-replicate group
        grp = (sid, rep)
        meta_snap = {c: row[c].strip() for c in consistent_cols}
        if grp not in group_meta:
            group_meta[grp] = {"first_row": i, **meta_snap}
        else:
            ref = group_meta[grp]
            mismatches = [c for c in consistent_cols if meta_snap[c] != ref[c]]
            if mismatches:
                errors += err(i,
                    f"Inconsistent metadata within tech-replicate group "
                    f"(sample_id='{sid}', replicate='{rep}'): "
                    f"columns {mismatches} differ from row {ref['first_row']}")

    for sid, ctrl in ip_to_control.items():
        if ctrl not in unique_sids:
            errors += err("?", f"Sample '{sid}' references control_id '{ctrl}' "
                               f"not found in samplesheet")

    if len(run_layouts) > 1:
        errors += err("?", "Mixed PE/SE runs are not supported: use separate "
                           "samplesheets, output directories and normalization cohorts")

    n_ip  = sum(1 for r in rows if r["is_control"].strip().lower() not in ("true","1","yes"))
    n_grp = len(group_meta)
    print(f"\nRows parsed          : {len(rows)}")
    print(f"Unique sample groups : {n_grp}  (sample_id + replicate)")
    print(f"IP sample rows       : {n_ip}")
    print(f"Control sample IDs   : {len(control_ids_declared)}")

    if errors:
        print(f"\n[FAIL] {errors} error(s) — fix before running workflow.",
              file=sys.stderr)
        sys.exit(1)
    print("\n[OK] Sample sheet is valid.")


if __name__ == "__main__":
    main()
