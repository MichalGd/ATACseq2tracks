#!/usr/bin/env python3
"""
fastq2tracks v3.0.2 — Sample sheet validator
Usage: python3 scripts/validate_samplesheet.py <samplesheet.csv> [--check-files]
"""
import csv, sys, os, argparse

REQUIRED_COLS = [
    "sample_id","fastq_1","fastq_2","layout","genome",
    "assay","factor","condition","treatment","cell_type",
    "replicate","tech_replicate","is_control","control_id",
    "macs2_mode","blacklist","chipqc_annotation","output_prefix"
]
VALID_LAYOUTS = {"PE","SE"}
VALID_GENOMES = {"hg38","mm39"}
VALID_MACS2   = {"narrow","broad","both","none"}
VALID_ISCTR   = {"TRUE","FALSE","true","false","1","0","yes","no"}

def err(row, msg):  print(f"  [ERROR] row {row}: {msg}", file=sys.stderr); return 1
def warn(row, msg): print(f"  [WARN]  row {row}: {msg}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("samplesheet")
    ap.add_argument("--check-files", action="store_true",
                    help="Verify FASTQ/BED/RDS files exist on disk")
    args = ap.parse_args()
    errors = 0
    print(f"Validating: {args.samplesheet}")

    with open(args.samplesheet, newline="") as fh:
        reader = csv.DictReader(fh)
        missing_cols = [c for c in REQUIRED_COLS if c not in (reader.fieldnames or [])]
        if missing_cols:
            print(f"[FATAL] Missing columns: {missing_cols}", file=sys.stderr); sys.exit(1)
        rows = list(reader)

    sample_ids, control_ids_declared, ip_to_control = {}, set(), {}
    for i, row in enumerate(rows, start=2):
        sid = row["sample_id"].strip()
        if sid in sample_ids: errors += err(i, f"Duplicate sample_id '{sid}'")
        sample_ids[sid] = i
        layout = row["layout"].strip().upper()
        if layout not in VALID_LAYOUTS: errors += err(i, f"Invalid layout '{layout}' (PE or SE)")
        if not row["fastq_1"].strip(): errors += err(i, "fastq_1 is empty")
        if layout == "PE" and not row["fastq_2"].strip(): errors += err(i, "PE sample missing fastq_2")
        if layout == "SE" and row["fastq_2"].strip(): warn(i, "SE sample has fastq_2 — will be ignored")
        if row["genome"].strip() not in VALID_GENOMES: errors += err(i, f"Invalid genome '{row['genome']}'")
        if row["macs2_mode"].strip().lower() not in VALID_MACS2: errors += err(i, f"Invalid macs2_mode '{row['macs2_mode']}'")
        is_ctr = row["is_control"].strip()
        if is_ctr not in VALID_ISCTR: errors += err(i, f"Invalid is_control '{is_ctr}'")
        is_ctr_bool = is_ctr.lower() in {"true","1","yes"}
        if is_ctr_bool:
            control_ids_declared.add(sid)
            if row["control_id"].strip(): warn(i, "Control sample should have empty control_id")
            if row["macs2_mode"].strip().lower() != "none": warn(i, "Control sample: macs2_mode should be 'none'")
        else:
            ctrl = row["control_id"].strip()
            if not ctrl: warn(i, "IP sample has no control_id — MACS2 will run without control")
            else: ip_to_control[sid] = ctrl
        try: int(row["replicate"])
        except ValueError: errors += err(i, f"replicate must be integer, got '{row['replicate']}'")
        if args.check_files:
            for col in ["fastq_1","fastq_2","blacklist","chipqc_annotation"]:
                p = row[col].strip()
                if p and not os.path.exists(p): errors += err(i, f"File not found: {col}={p}")

    for sid, ctrl in ip_to_control.items():
        if ctrl not in sample_ids:
            errors += err("?", f"Sample '{sid}' references control_id '{ctrl}' not in samplesheet")

    print(f"\nSamples parsed  : {len(rows)}")
    print(f"IP samples      : {len([r for r in rows if r['is_control'].lower() not in ('true','1','yes')])}")
    print(f"Control samples : {len(control_ids_declared)}")
    if errors:
        print(f"\n[FAIL] {errors} error(s) — fix before running workflow.", file=sys.stderr); sys.exit(1)
    print("\n[OK] Sample sheet is valid.")

if __name__ == "__main__":
    main()
