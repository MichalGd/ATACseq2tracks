#!/usr/bin/env python3
"""Safely normalize Windows text artifacts in workflow input files.

Affected files are backed up byte-for-byte beside the original and then
atomically replaced with UTF-8 text using LF line endings. Clean UTF-8/LF files
are never rewritten and do not receive a backup.
"""

from __future__ import annotations

import argparse
import codecs
import os
import shutil
import stat
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


class SanitizationError(RuntimeError):
    """Raised when an input cannot be normalized safely."""


def decode_and_normalize(payload: bytes) -> tuple[bytes, list[str]]:
    artifacts: list[str] = []

    # UTF-32 BOMs begin with the same bytes as UTF-16 BOMs, so test them first.
    if payload.startswith((codecs.BOM_UTF32_LE, codecs.BOM_UTF32_BE)):
        text = payload.decode("utf-32")
        artifacts.append("UTF-32 BOM/encoding")
    elif payload.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        text = payload.decode("utf-16")
        artifacts.append("UTF-16 BOM/encoding")
    elif payload.startswith(codecs.BOM_UTF8):
        text = payload.decode("utf-8-sig")
        artifacts.append("UTF-8 BOM")
    else:
        if b"\x00" in payload:
            raise SanitizationError(
                "NUL bytes detected without a recognized Unicode BOM; "
                "refusing to guess the encoding"
            )
        try:
            text = payload.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise SanitizationError(
                "file is not valid UTF-8 and has no supported UTF BOM"
            ) from exc

    if "\r\n" in text:
        artifacts.append("CRLF line endings")
    without_crlf = text.replace("\r\n", "\n")
    if "\r" in without_crlf:
        artifacts.append("bare CR line endings")
    normalized = without_crlf.replace("\r", "\n")

    if normalized.endswith("\x1a"):
        normalized = normalized[:-1]
        artifacts.append("DOS EOF marker (Ctrl-Z)")

    return normalized.encode("utf-8"), artifacts


def unique_backup_path(path: Path) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return path.with_name(
        f"{path.name}.windows-artifact-backup.{stamp}.{os.getpid()}"
    )


def atomic_replace(path: Path, payload: bytes, mode: int) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.sanitize.", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def sanitize(path_argument: str) -> bool:
    supplied_path = Path(path_argument).expanduser()
    try:
        path = supplied_path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise SanitizationError(f"file not found: {supplied_path}") from exc
    if not path.is_file():
        raise SanitizationError(f"not a regular file: {path}")

    original = path.read_bytes()
    normalized, artifacts = decode_and_normalize(original)
    if not artifacts:
        print(f"[INPUT] OK: {path} (UTF-8/LF; unchanged)")
        return False

    original_stat = path.stat()
    backup = unique_backup_path(path)
    # The backup must succeed before the original is changed.
    shutil.copy2(path, backup)
    try:
        atomic_replace(path, normalized, stat.S_IMODE(original_stat.st_mode))
    except Exception as exc:
        raise SanitizationError(
            f"normalization failed; original backup is preserved at {backup}"
        ) from exc

    print(f"[INPUT] Corrected: {path}")
    print(f"[INPUT] Artifacts: {', '.join(artifacts)}")
    print(f"[INPUT] Original backup: {backup}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Back up and normalize Windows text artifacts in-place"
    )
    parser.add_argument("files", nargs="+", help="Configuration/CSV files")
    args = parser.parse_args()

    failed = False
    for filename in args.files:
        try:
            sanitize(filename)
        except (OSError, SanitizationError) as exc:
            print(f"[INPUT] ERROR: {filename}: {exc}", file=sys.stderr)
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
