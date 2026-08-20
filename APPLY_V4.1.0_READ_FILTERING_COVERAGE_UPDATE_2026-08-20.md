# Apply ATACseq2tracks v4.1.0

Do not overwrite the active shared installation in place. Install the complete
repository in a versioned directory, validate it, then update the `current`
symlink atomically.

## Shared deployment outline

```bash
sudo mkdir -p /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0
sudo cp -a /absolute/path/to/ATACseq2tracks_v4.1.0/. \
  /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0/

sudo find /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0 \
  -type d -exec chmod 755 {} +
sudo find /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0/scripts \
  /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0/tests \
  /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0/utilities \
  -type f -name '*.sh' -exec chmod 755 {} +
sudo chmod 755 \
  /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0/atacseq2tracks.sh \
  /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0/fastq2tracks.sh
sudo find /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0 \
  -type f ! -name '*.sh' -exec chmod 644 {} +
sudo chown -R root:root /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0
```

Run static checks and representative PE/SE validation before switching:

```bash
cd /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0
bash tests/check_bash_syntax.sh

bash scripts/smoke_test.sh \
  /absolute/user-owned/path/samplesheet.csv \
  /absolute/user-owned/path/config.conf
```

After validation:

```bash
sudo ln -sfn /opt/bioinformatics/workflows/ATACseq2tracks/4.1.0 \
  /opt/bioinformatics/workflows/ATACseq2tracks/current.new
sudo mv -Tf /opt/bioinformatics/workflows/ATACseq2tracks/current.new \
  /opt/bioinformatics/workflows/ATACseq2tracks/current
```

Verify as a non-owner account that `VERSION` reads `4.1.0`, all scripts are
readable/executable, preflight passes, and output paths remain user-owned.

Existing v4.0.0 project configs continue to run because code fallbacks enable
the new families. To preserve v4.0.0 output volume when intentionally reusing an
old config, explicitly set the new permissive/intermediate switches and CPM
bedGraph format switch to `false`.
