"""Provenance tracking — SHA-256 checksums and source lineage."""

import hashlib
import json
from pathlib import Path


def compute_checksum(file_path: Path) -> str:
    """Compute SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def build_provenance(raw_dir: Path) -> dict:
    """Build provenance record for all JSON files in raw_dir."""
    files = sorted(raw_dir.glob("*.json"))
    checksums = {}
    for fp in files:
        checksums[fp.name] = compute_checksum(fp)
    return {
        "source_directory": str(raw_dir),
        "file_count": len(files),
        "checksums": checksums,
    }


def verify_provenance(provenance: dict, raw_dir: Path) -> list:
    """Verify checksums match current files. Returns list of failures."""
    failures = []
    checksums = provenance.get("checksums", {})
    for fname, expected_hash in checksums.items():
        fp = raw_dir / fname
        if not fp.exists():
            failures.append(f"Missing file: {fname}")
            continue
        actual = compute_checksum(fp)
        if actual != expected_hash:
            failures.append(f"Checksum mismatch: {fname} (expected {expected_hash[:12]}..., got {actual[:12]}...)")
    return failures
