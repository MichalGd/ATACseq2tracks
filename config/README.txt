fastq2tracks v3.0.3 — Config / Samplesheet directory
=====================================================

Files:
  config.conf                Main configuration — the ONLY file you edit per run.
  samplesheet_template.csv   Column headers only; copy and fill in.
  samplesheet_example.csv    Filled-in examples (hg38 PE + mm39 SE).
  README.txt                 This file.

How to set up a new project
─────────────────────────────────────────────────────
1. Create a project config folder anywhere:
     mkdir -p /home/USER/myproject/config

2. Copy and edit config.conf:
     cp /path/to/fastq2tracks/config/config.conf /home/USER/myproject/config/
     nano /home/USER/myproject/config/config.conf
     # Change: SAMPLESHEET, OUTPUT_DIR, and THREADS_PARALLEL_JOBS as needed

3. Prepare your samplesheet (path must match SAMPLESHEET in config.conf):
     cp /path/to/fastq2tracks/config/samplesheet_template.csv \
        /home/USER/myproject/config/samplesheet.csv
     # Fill in one row per FASTQ file / tech-rep pair

4. Run — only ONE argument needed:
     bash /path/to/fastq2tracks/fastq2tracks.sh \
          --config /home/USER/myproject/config/config.conf

Multiple users on the same server
──────────────────────────────────
Each user keeps their own config.conf + samplesheet.csv anywhere they like.
The scripts/ directory is installed once and shared read-only.
No user-specific paths exist in any script.
