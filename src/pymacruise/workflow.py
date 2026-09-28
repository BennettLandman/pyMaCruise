"""Typed workflow settings and provenance manifest helpers."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.contracts import ValidatedInputs, WorkflowProfile

LOGGER = logging.getLogger("pymacruise")


@dataclass(frozen=True)
class WorkflowConfig:
    """Validated input paths, profile, and source-derived defaults."""

    profile: WorkflowProfile | str
    t1_path: Path | str
    maseg_path: Path | str
    brain_mask_path: Path | str
    output_dir: Path | str
    flair_path: Path | str | None = None
    bias_correction: str = "N4"

    def __post_init__(self) -> None:
        try:
            profile = WorkflowProfile(self.profile)
        except ValueError as exc:
            raise ValueError(f"Unsupported workflow profile: {self.profile!r}") from exc
        object.__setattr__(self, "profile", profile)
        for name in ("t1_path", "maseg_path", "brain_mask_path", "output_dir"):
            object.__setattr__(self, name, Path(getattr(self, name)))
        if self.flair_path is not None:
            object.__setattr__(self, "flair_path", Path(self.flair_path))
        for name in ("t1_path", "maseg_path", "brain_mask_path", "flair_path"):
            path = getattr(self, name)
            if path is not None and not path.name.lower().endswith(".nii.gz"):
                raise ValueError(f"{name} must use compressed NIfTI (.nii.gz)")
        if profile is WorkflowProfile.PLUS and self.flair_path is None:
            raise ValueError("The plus profile requires flair_path")
        if profile is WorkflowProfile.BASE and self.flair_path is not None:
            raise ValueError("The base profile does not accept flair_path")
        if self.bias_correction != "N4":
            raise ValueError("Only the source-selected N4 bias correction is supported")

    def as_dict(self) -> dict[str, Any]:
        result = asdict(self)
        return {
            key: (str(value) if isinstance(value, Path) else value.value if isinstance(value, WorkflowProfile) else value)
            for key, value in result.items()
        }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def create_initial_manifest(
    config: WorkflowConfig,
    validated: ValidatedInputs,
    *,
    source_revision: str | None = None,
) -> dict[str, Any]:
    """Create a no-op run manifest after input validation, before processing."""
    if validated.profile is not config.profile:
        raise ValueError("Validated input profile does not match workflow configuration")
    paths = {
        "t1": config.t1_path,
        "maseg": config.maseg_path,
        "brain_mask": config.brain_mask_path,
    }
    if config.flair_path is not None:
        paths["flair"] = config.flair_path
    input_records: dict[str, Any] = {}
    images = {
        "t1": validated.t1,
        "maseg": validated.maseg,
        "brain_mask": validated.brain_mask,
        "flair": validated.flair,
    }
    for role, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Cannot write run manifest; missing {role} input: {path}")
        image = images[role]
        assert image is not None
        input_records[role] = {
            "path": str(path),
            "sha256": _sha256(path),
            "shape": list(image.shape),
            "dtype": str(image.get_data_dtype()),
            "affine_ras": np.asarray(image.affine, dtype=np.float64).tolist(),
            "input_axis_codes": list(validated.input_orientations[role]),
            "output_axis_codes": list(nib.aff2axcodes(image.affine)),
        }

    config_data = config.as_dict()
    signature_data = {
        "config": config_data,
        "inputs": {role: row["sha256"] for role, row in input_records.items()},
    }
    signature = hashlib.sha256(
        json.dumps(signature_data, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    dependencies = {
        package: _package_version(package)
        for package in ("pymacruise", "numpy", "nibabel", "scipy", "antspyx")
    }
    manifest = {
        "schema_version": 1,
        "run_signature": signature,
        "workflow_profile": config.profile.value,
        "config": config_data,
        "source_revision": source_revision,
        "dependencies": dependencies,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "inputs": input_records,
        "transforms": [],
        "stages": {
            "input_validation": "passed",
            "pipeline": "not_started",
        },
        "outputs": {},
        "failures": [],
    }
    LOGGER.info("Initialized %s workflow run %s", config.profile.value, signature[:12])
    return manifest


def write_run_manifest(path: str | Path, manifest: dict[str, Any]) -> Path:
    """Write a run manifest as JSON under the caller's output directory."""
    destination = Path(path)
    if destination.suffix.lower() != ".json":
        raise ValueError("Run manifests must use the .json suffix")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    LOGGER.info("Wrote run manifest to %s", destination)
    return destination
