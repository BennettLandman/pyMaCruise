"""Resolve the atlas bundle in installed packages and source checkouts.

Wheels stage the single repository atlas tree as ``pymacruise/data/atlases``
during build. In a source checkout, this resolver falls back to the canonical
repository-level ``data/atlases`` directory. It never searches neighboring
VandyCRUISE, MIPAV, or JIST checkouts.
"""

from __future__ import annotations

from importlib.resources import files
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
    root = get_bundled_atlas_root().resolve()
    descriptor = (root / relative).resolve()
    if not descriptor.is_relative_to(root):
        raise AtlasResourceError(f"Atlas descriptor path escapes the bundle: {relative}")
    if not descriptor.is_file():
        raise AtlasResourceError(f"Bundled atlas descriptor does not exist: {relative}")
    return descriptor
