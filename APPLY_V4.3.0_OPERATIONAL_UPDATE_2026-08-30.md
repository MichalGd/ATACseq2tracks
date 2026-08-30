# Apply ATACseq2tracks v4.3.0

## Review and validate locally

```bash
cd /path/to/ATACseq2tracks_v4.3.0
bash tests/check_bash_syntax.sh
bin/atacseq2tracks --version
```

Review `VERSION`, `CHANGELOG.md`, `docs/v4.3.0_OPERATIONAL_UPDATE_2026-08-30.md`
and `V4.3.0_OPERATIONAL_MANIFEST_2026-08-30.tsv` before deployment.

## Install for all users

The installer copies rather than mutates the reviewed source into a staging
directory, runs regression checks, promotes the tested copy to a root-owned
immutable release, and only then updates the stable workflow and command links
atomically:

```bash
sudo bash /absolute/path/to/ATACseq2tracks_v4.3.0/utilities/install_shared_release.sh
```

Default runtime paths are:

```text
/opt/bioinformatics/workflows/ATACseq2tracks/releases/4.3.0
/opt/bioinformatics/workflows/ATACseq2tracks/current
/usr/local/bin/atacseq2tracks
/opt/miniconda/envs/ATACseq2tracks
/opt/miniconda/envs/ataqv-tools
```

## Validate as an ordinary user

```bash
command -v atacseq2tracks
atacseq2tracks --version
atacseq2tracks --config /absolute/path/user/config.conf --plan
atacseq2tracks --config /absolute/path/user/config.conf --preflight-only
```

No user environment activation is needed. Do not overwrite or delete v4.2.0;
retaining the prior versioned release provides a rollback target. To roll back,
atomically repoint `current` and `/usr/local/bin/atacseq2tracks` to the reviewed
prior release.
