# How to upload this repository to GitHub

This folder is ready to become a GitHub repository. Choose one method below.

## Before upload

1. Unzip the package.
2. Rename the folder if desired.
3. Add a license if the repository will be public.
4. Search for private paths, usernames, server URLs, or project names you do not want public:

```bash
grep -R -E "micgdu|dysk2|your-server|http://" -n .
```

5. Decide whether to keep the PowerPoint deck in `presentations/`. It may contain screenshots or local file paths from previous analyses.

## Method A. Upload using the GitHub website

1. Sign in to GitHub.
2. Create a new repository.
3. Do not initialize it with a README, `.gitignore`, or license if you plan to upload this prepared folder as-is.
4. Open the new empty repository.
5. Use **Add file** and then **Upload files**.
6. Drag the contents of this folder into the upload area.
7. Commit the upload.

This method is simple for small repositories. The command-line method is more reliable for many files.

## Method B. Upload with Git command line

From inside the unzipped folder:

```bash
git init
git add .
git commit -m "Initial fastq2tracks workflow documentation"
git branch -M main
git remote add origin https://github.com/<OWNER>/<REPOSITORY>.git
git push -u origin main
```

Replace `<OWNER>` and `<REPOSITORY>` with your GitHub user or organization and repository name.

## Method C. Upload with GitHub CLI

Install and authenticate GitHub CLI, then run:

```bash
gh auth login
git init
git add .
git commit -m "Initial fastq2tracks workflow documentation"
gh repo create <OWNER>/<REPOSITORY> --source=. --public --push
```

Use `--private` instead of `--public` if the repository should be private.

## After upload

1. Open the repository on GitHub.
2. Confirm that `README.md` renders correctly.
3. Confirm that the schematic image appears.
4. Open `docs/KNOWN_ISSUES.md` and decide which deployment notes should be fixed in code before release.
5. Add repository topics such as `bioinformatics`, `illumina`, `bowtie2`, `fastqc`, `multiqc`, `bigwig`, `genomics`.
6. Add a release tag if this will be versioned:

```bash
git tag v2.1
git push origin v2.1
```

## Large-file caution

Do not commit generated BAM, FASTQ, bedGraph, BigWig, MultiQC output folders, or full run outputs unless you intentionally use Git LFS or another data repository. The included `.gitignore` excludes typical generated analysis files.
