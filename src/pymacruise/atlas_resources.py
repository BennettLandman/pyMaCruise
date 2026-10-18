"""Resolve the atlas bundle in installed packages and source checkouts.

Wheels stage the single repository atlas tree as ``pymacruise/data/atlases``
during build. In a source checkout, this resolver falls back to the canonical
repository-level ``data/atlases`` directory. It never searches neighboring
VandyCRUISE, MIPAV, or JIST checkouts.
"""

from __future__ import annotations

from importlib.resources import files
import gzip
import hashlib
import json
from pathlib import Path


class AtlasResourceError(FileNotFoundError):
    """Raised when a bundled atlas resource cannot be found or resolved."""


def get_bundled_atlas_root() -> Path:
    """Return the installed or source-checkout atlas directory as a filesystem path."""

    packaged = files("pymacruise").joinpath("data", "atlases")
    if packaged.is_dir():
        # Standard wheel installations expose importlib resources as real paths.
        # The public loaders accept Path objects and perform relative NIfTI lookup.
        if isinstance(packaged, Path):
            return packaged
        raise AtlasResourceError(
            "The atlas resources are inside a non-filesystem package importer; "
            "install PyMaCRUISE from its wheel to access them as paths."
        )

    source_checkout = Path(__file__).resolve().parents[2] / "data" / "atlases"
    if source_checkout.is_dir():
        return source_checkout
    raise AtlasResourceError(
        "Bundled atlases are unavailable. Install the PyMaCRUISE wheel or provide "
        "an atlas descriptor from a source checkout containing data/atlases."
    )


def get_bundled_atlas_descriptor(relative_path: str | Path) -> Path:
    """Resolve a bundled descriptor under the atlas root, rejecting path escapes."""

    relative = Path(relative_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise AtlasResourceError(f"Atlas descriptor path must be relative to the bundle: {relative}")
    if relative.suffix.lower() != ".txt":
        raise AtlasResourceError(f"Atlas descriptor must be a .txt file: {relative}")
    return get_bundled_atlas_asset(relative)


def get_bundled_atlas_asset(relative_path: str | Path) -> Path:
    """Resolve any installed atlas asset, including preserved native resources."""
    relative = Path(relative_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise AtlasResourceError(f"Atlas asset path must be relative to the bundle: {relative}")
    root = get_bundled_atlas_root().resolve()
    asset = (root / relative).resolve()
    if not asset.is_relative_to(root):
        raise AtlasResourceError(f"Atlas asset path escapes the bundle: {relative}")
    if not asset.is_file():
        raise AtlasResourceError(f"Bundled atlas asset does not exist: {relative}")
    return asset


def read_bundled_atlas_source(atlas: str, filename: str) -> bytes:
    """Read and checksum an original source file from a preserved atlas bundle.

    Native voxels retain their original indexing and metadata. Reading them
    does not convert them to RAS or enable them as an AtlasEM runtime atlas.
    Neither the checkout nor a neighboring legacy installation is consulted.
    """
    if (atlas in {"", ".", ".."} or filename in {"", ".", ".."}
            or Path(atlas).name != atlas or Path(filename).name != filename):
        raise AtlasResourceError("Atlas and source filename must each be a single path component")
    manifest = json.loads(get_bundled_atlas_asset(f"{atlas}/SOURCE_MANIFEST.json").read_text())
    try:
        record = manifest["files"][filename]
    except KeyError as exc:
        raise AtlasResourceError(f"Unknown source asset: {atlas}/{filename}") from exc
    path = get_bundled_atlas_asset(Path(atlas) / record["stored_path"])
    stored = path.read_bytes()
    if hashlib.sha256(stored).hexdigest() != record["stored_sha256"]:
        raise AtlasResourceError(f"Stored asset checksum mismatch: {path}")
    if record["encoding"] not in {"gzip", "identity"}:
        raise AtlasResourceError(f"Unsupported source encoding: {record['encoding']}")
    payload = gzip.decompress(stored) if record["encoding"] == "gzip" else stored
    if len(payload) != record["source_bytes"] or hashlib.sha256(payload).hexdigest() != record["source_sha256"]:
        raise AtlasResourceError(f"Source asset checksum mismatch: {path}")
    return payload
