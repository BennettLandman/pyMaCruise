"""Typed workflow settings and provenance manifest helpers."""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import importlib.util
import json
import logging
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np

from pymacruise.contracts import ValidatedInputs, WorkflowProfile
from pymacruise.contracts import validate_workflow_files
from pymacruise.plus_pipeline import PlusPipelineOptions
from pymacruise.atlas_resources import get_bundled_atlas_descriptor


_WORKFLOW_CONFIG_KEYS = {
    "schema_version",
    "profile",
    "t1_path",
    "maseg_path",
    "brain_mask_path",
    "output_dir",
    "flair_path",
    "bias_correction",
    "plus_options",
}
_ATLAS_DESCRIPTORS = {
    WorkflowProfile.BASE: "cruise-atlas-10obj-toads2008.txt",
    WorkflowProfile.PLUS: (
        "lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
    ),
}


class WorkflowPreflightError(RuntimeError):
    """Raised when a dry-run workflow preflight cannot be completed."""

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
    plus_options: PlusPipelineOptions | None = None

    def __post_init__(self) -> None:
        try:
            profile = WorkflowProfile(self.profile)
        except ValueError as exc:
            raise ValueError(f"Unsupported workflow profile: {self.profile!r}") from exc
        object.__setattr__(self, "profile", profile)
        if profile is WorkflowProfile.PLUS:
            options = self.plus_options or PlusPipelineOptions()
            if not isinstance(options, PlusPipelineOptions):
                raise ValueError("plus_options must be PlusPipelineOptions")
            object.__setattr__(self, "plus_options", options)
        elif self.plus_options is not None:
            raise ValueError("plus_options are only valid for the plus profile")
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
        if self.plus_options is None:
            # Keep the established base manifest/config schema unchanged.
            result.pop("plus_options")
        else:
            result["plus_options"] = self.plus_options.as_dict()
        return {
            key: (str(value) if isinstance(value, Path) else value.value if isinstance(value, WorkflowProfile) else value)
            for key, value in result.items()
        }


def load_workflow_config_file(path: str | Path) -> WorkflowConfig:
    """Read the version-1 JSON workflow config and resolve paths beside it.

    The config file contains input/output paths and profile settings only; it
    does not select or override bundled atlas resources.
    """
    config_path = Path(path).expanduser().resolve()
    try:
        raw = json.loads(config_path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read workflow config {config_path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError("Workflow config must be a JSON object")
    unknown = sorted(set(raw) - _WORKFLOW_CONFIG_KEYS)
    if unknown:
        raise ValueError(f"Unknown workflow config field(s): {', '.join(unknown)}")
    required = {
        "schema_version", "profile", "t1_path", "maseg_path", "brain_mask_path", "output_dir"
    }
    missing = sorted(required - set(raw))
    if missing:
        raise ValueError(f"Missing workflow config field(s): {', '.join(missing)}")
    if type(raw["schema_version"]) is not int or raw["schema_version"] != 1:
        raise ValueError(f"Unsupported workflow config schema_version: {raw['schema_version']!r}")

    def resolve_config_path(value: Any, name: str) -> Path | None:
        if value is None and name == "flair_path":
            return None
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"Workflow config {name} must be a nonempty path string")
        candidate = Path(value).expanduser()
        return candidate.resolve() if candidate.is_absolute() else (config_path.parent / candidate).resolve()

    plus_options = raw.get("plus_options")
    if plus_options is not None:
        if not isinstance(plus_options, dict):
            raise ValueError("Workflow config plus_options must be a JSON object")
        option_names = set(PlusPipelineOptions.__dataclass_fields__)
        unknown_options = sorted(set(plus_options) - option_names)
        if unknown_options:
            raise ValueError(
                f"Unknown plus_options field(s): {', '.join(unknown_options)}"
            )
        plus_options = dict(plus_options)
        if "boundary_overlook_labels" in plus_options:
            labels = plus_options["boundary_overlook_labels"]
            if not isinstance(labels, list):
                raise ValueError("plus_options.boundary_overlook_labels must be a JSON array")
            plus_options["boundary_overlook_labels"] = tuple(labels)
        try:
            plus_options = PlusPipelineOptions(**plus_options)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid workflow config plus_options: {exc}") from exc

    return WorkflowConfig(
        profile=raw["profile"],
        t1_path=resolve_config_path(raw["t1_path"], "t1_path"),
        maseg_path=resolve_config_path(raw["maseg_path"], "maseg_path"),
        brain_mask_path=resolve_config_path(raw["brain_mask_path"], "brain_mask_path"),
        output_dir=resolve_config_path(raw["output_dir"], "output_dir"),
        flair_path=resolve_config_path(raw.get("flair_path"), "flair_path"),
        bias_correction=raw.get("bias_correction", "N4"),
        plus_options=plus_options,
    )


def _ants_available() -> bool:
    """Return whether the optional ANTsPy module imports successfully."""
    try:
        if importlib.util.find_spec("ants") is None:
            return False
        importlib.import_module("ants")
        return True
    except (ImportError, ValueError):
        return False


def _workflow_atlas_preflight(profile: WorkflowProfile) -> dict[str, Any]:
    """Validate the bundled profile descriptor and referenced NIfTI headers."""
    descriptor = get_bundled_atlas_descriptor(_ATLAS_DESCRIPTORS[profile])
    atlas_root = descriptor.parent.resolve()
    assets: list[dict[str, Any]] = []
    for line in descriptor.read_text(encoding="utf-8").splitlines():
        if not line.startswith("File:"):
            continue
        relative = Path(line.partition(":")[2].strip())
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe atlas asset path in {descriptor}: {relative}")
        asset_path = (atlas_root / relative).resolve()
        if not asset_path.is_relative_to(atlas_root):
            raise ValueError(f"Atlas asset path escapes descriptor directory: {relative}")
        if asset_path.suffixes[-2:] != [".nii", ".gz"] or not asset_path.is_file():
            raise FileNotFoundError(f"Missing or uncompressed atlas asset: {asset_path}")
        try:
            image = nib.load(str(asset_path))
        except (OSError, nib.filebasedimages.ImageFileError) as exc:
            raise ValueError(f"Cannot load atlas asset header {asset_path}: {exc}") from exc
        if len(image.shape) != 3 or nib.aff2axcodes(image.affine) != ("R", "A", "S"):
            raise ValueError(
                f"Atlas asset must be a 3D RAS NIfTI image: {asset_path}"
            )
        assets.append(
            {
                "path": str(asset_path),
                "shape": [int(size) for size in image.shape],
                "axis_codes": list(nib.aff2axcodes(image.affine)),
            }
        )
    if not assets:
        raise ValueError(f"Atlas descriptor contains no NIfTI assets: {descriptor}")
    return {
        "descriptor": str(descriptor),
        "descriptor_sha256": _sha256(descriptor),
        "assets_checked": len(assets),
        "assets": assets,
        "payloads_materialized": False,
    }


def _check_workflow_dependencies(config: WorkflowConfig) -> dict[str, str]:
    """Check importability of core requirements and selected optional N4."""
    dependencies = {
        package: _package_version(package)
        for package in ("pymacruise", "numpy", "nibabel", "scipy", "scikit-image")
    }
    missing_core = [
        package
        for package in ("numpy", "nibabel", "scipy", "skimage")
        if importlib.util.find_spec(package) is None
    ]
    if missing_core:
        raise WorkflowPreflightError(
            "Missing core PyMaCRUISE dependency/dependencies: "
            + ", ".join(missing_core)
            + ". Reinstall PyMaCRUISE in this environment."
        )
    if config.bias_correction == "N4" and not _ants_available():
        raise WorkflowPreflightError(
            "The selected workflow requires N4, but ANTsPy is not available. "
            "Install it with `python -m pip install 'pymacruise[ants]'` and retry."
        )
    if _ants_available():
        try:
            dependencies["antspyx"] = importlib.metadata.version("antspyx")
        except importlib.metadata.PackageNotFoundError:
            dependencies["antspyx"] = "importable; distribution version unavailable"
    else:
        dependencies["antspyx"] = "not-installed"
    return dependencies


def preflight_workflow(
    config: WorkflowConfig,
    *,
    overwrite: bool = False,
    source_revision: str | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Validate workflow prerequisites and write a manifest without execution.

    This checks image inputs, the bundled profile atlas, dependencies, and
    destination behavior. It deliberately invokes no algorithm stage.
    """
    dependencies = _check_workflow_dependencies(config)
    validated = validate_workflow_files(
        config.profile,
        config.t1_path,
        config.maseg_path,
        config.flair_path,
        brain_mask_path=config.brain_mask_path,
    )
    atlas = _workflow_atlas_preflight(config.profile)

    output_dir = Path(config.output_dir).expanduser()
    if output_dir.exists() and not output_dir.is_dir():
        raise NotADirectoryError(f"Workflow output path is not a directory: {output_dir}")
    manifest_path = output_dir / "run.json"
    if manifest_path.exists() and not overwrite:
        raise FileExistsError(
            f"Preflight manifest already exists: {manifest_path}; use --overwrite to replace it"
        )
    if output_dir.is_dir():
        other_entries = [entry for entry in output_dir.iterdir() if entry != manifest_path]
        if other_entries and not manifest_path.exists():
            raise FileExistsError(
                f"Workflow output directory is not empty: {output_dir}; choose an empty directory"
            )

    manifest = create_initial_manifest(
        config, validated, source_revision=source_revision
    )
    manifest["dependencies"].update(dependencies)
    manifest["dry_run"] = True
    manifest["configuration_schema_version"] = 1
    manifest["preflight"] = {
        "status": "passed",
        "mode": "validation_only",
        "workflow_executed": False,
        "dependencies": dependencies,
        "atlas": atlas,
        "orientation_conversions": {
            role: {
                "from": list(codes),
                "to": ["R", "A", "S"],
                "reoriented": tuple(codes) != ("R", "A", "S"),
            }
            for role, codes in validated.input_orientations.items()
        },
        "notes": [
            "No N4 correction or segmentation stage was executed.",
            "Atlas image headers were checked; atlas voxel payloads were not loaded.",
        ],
    }
    manifest["stages"]["pipeline"] = "not_started"
    output_dir.mkdir(parents=True, exist_ok=True)
    write_run_manifest(manifest_path, manifest)
    return manifest_path, manifest


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
