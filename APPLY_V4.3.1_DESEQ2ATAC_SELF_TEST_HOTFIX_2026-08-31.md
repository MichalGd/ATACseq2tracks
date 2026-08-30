# Apply ATACseq2tracks v4.3.1

This release corrects only the v4.3.0 DESeq2ATAC synthetic self-test. Production
analysis behavior is unchanged.

## Validate the source

```bash
cd /absolute/path/to/ATACseq2tracks_v4.3.1
bash tests/check_bash_syntax.sh
```

The validation must finish without `FAIL` lines and must include:

```text
OK   DESeq2ATAC synthetic significant and zero-significant analyses
OK   differential accessibility regression checks
```

## Install the immutable shared release

```bash
sudo bash /absolute/path/to/ATACseq2tracks_v4.3.1/utilities/install_shared_release.sh
```

The installer validates the staged copy before atomically changing `current`.
After success, verify:

```bash
readlink -f /opt/bioinformatics/workflows/ATACseq2tracks/current
cat /opt/bioinformatics/workflows/ATACseq2tracks/current/VERSION
/usr/local/bin/atacseq2tracks --version
```

All three checks should identify version `4.3.1`.
