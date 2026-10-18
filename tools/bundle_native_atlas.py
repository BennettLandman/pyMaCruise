"""Preserve supplied ibsr11objV3 source assets without inventing geometry.

Run: python tools/bundle_native_atlas.py /path/to/Atlas/ibsr11objV3
This maintenance tool is never called during installation or algorithm runs.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path


def bundle(source: Path, destination: Path) -> dict:
    source = source.resolve()
    if source.name != "ibsr11objV3":
        raise ValueError("Expected the supplied ibsr11objV3 source folder")
    sources = sorted(p for p in source.iterdir() if p.suffix in {".raw", ".xml", ".txt"})
    descriptor = source / "LesionToads-11-obj-atlasV3.txt"
    revisions = {}
    entries = source / "CVS/Entries"
    if entries.is_file():
        for row in entries.read_text().splitlines():
            parts = row.split("/")
            if len(parts) >= 4 and parts[0] == "":
                revisions[parts[1]] = {"revision": parts[2], "timestamp": parts[3]}
    references = [line.partition(":")[2].strip() for line in descriptor.read_text().splitlines()
                  if line.startswith("File:")]
    for reference in references:
        path = Path(reference)
        if path.parts[0] != source.name or len(path.parts) != 2 or not (source / path.name).is_file():
            raise ValueError(f"Missing or unsafe source payload: {reference}")
    manifest = {
        "schema_version": 1, "atlas": "ibsr11objV3",
        "provenance": "Supplied VandyCRUISE/src/Atlas/ibsr11objV3; per-file CVS revisions recorded below",
        "descriptor": descriptor.name, "referenced_payloads": references,
        "storage": "Lossless original bytes; native indexing; no affine assignment or resampling",
        "runtime_status": "Source assets available; AtlasEM descriptor adaptation and checked RAS conversion remain open",
        "files": {},
    }
    for path in sources:
        payload = path.read_bytes()
        compressed = True  # Preserve XML line endings across Git autocrlf settings too.
        stored_path = Path("native") / (path.name + ".gz" if compressed else path.name)
        output = destination / stored_path
        output.parent.mkdir(parents=True, exist_ok=True)
        if compressed:
            stream = io.BytesIO()
            with gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0) as writer:
                writer.write(payload)
            stored = stream.getvalue()
        else:
            stored = payload
        output.write_bytes(stored)
        manifest["files"][path.name] = {
            "stored_path": stored_path.as_posix(), "encoding": "gzip" if compressed else "identity",
            "source_bytes": len(payload), "source_sha256": hashlib.sha256(payload).hexdigest(),
            "stored_bytes": len(stored), "stored_sha256": hashlib.sha256(stored).hexdigest(),
            "source_cvs": revisions.get(path.name),
        }
    (destination / "SOURCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    destination = Path(__file__).resolve().parents[1] / "data/atlases/ibsr11objV3"
    result = bundle(args.source, destination)
    print(f"Preserved {len(result['files'])} source files, including {len(result['referenced_payloads'])} referenced payloads")
