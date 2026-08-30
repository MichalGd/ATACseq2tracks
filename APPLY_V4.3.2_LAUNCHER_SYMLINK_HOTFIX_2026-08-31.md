# Apply ATACseq2tracks v4.3.2

Version 4.3.2 corrects activation-free invocation through
`/usr/local/bin/atacseq2tracks` and retains the v4.3.1 self-test fix.

```bash
cd /absolute/path/to/ATACseq2tracks_v4.3.2
sudo bash utilities/install_shared_release.sh
```

The installer runs the complete staged validation suite before atomically
changing the shared links. After success, verify both direct and shared launch:

```bash
/opt/bioinformatics/workflows/ATACseq2tracks/current/bin/atacseq2tracks --version
/usr/local/bin/atacseq2tracks --version
sudo -u kmk /usr/local/bin/atacseq2tracks --version
```

All commands should print `4.3.2`.
