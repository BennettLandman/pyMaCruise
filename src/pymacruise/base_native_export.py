"""Restore BASE deliverables to the original input lattice and physical frame."""
from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.processing import resample_from_to

from pymacruise.image_io import normalize_image_to_ras
from pymacruise.surfaces import read_legacy_vtk, transform_surface_ras, write_legacy_vtk


NATIVE_OUTPUT_NAMES = (
    "inner.nii.gz", "central.nii.gz", "outer.nii.gz", "SEG_refine.nii.gz",
    "inner_surface.vtk", "central_surface.vtk", "outer_surface.vtk",
    "inner_surface_thickness.vtk", "thickness_levelset.nii.gz",
    "thickness_summary.json", "provenance.json",
)


def _on_original_reference(reference, values, dtype):
    """Retain the original storage axes, affine and both spatial forms."""
    header = reference.header.copy()
    header.set_data_dtype(dtype)
    header.set_slope_inter(None, None)
    kind = nib.Nifti2Image if isinstance(reference, nib.Nifti2Image) else nib.Nifti1Image
    output = kind(np.asarray(values, dtype=dtype), reference.affine, header=header)
    output.set_qform(reference.get_qform(), code=int(reference.header["qform_code"]))
    output.set_sform(reference.get_sform(), code=int(reference.header["sform_code"]))
    return output


def restore_image_to_original_grid(
    image, original_reference, working_reference, *, orientation_only=False,
    labels=False, outside_value=0.0,
):
    """Use physical resampling, or invert the axis-only compatibility mapping.

    Label sampling is nearest neighbor. Continuous maps use linear sampling.
    The orientation-only route reverses permutations/flips without interpolation
    and reinstates the original direction cosines; interpreting its discarded
    oblique affine as physical geometry would rotate the voxel data incorrectly.
    """
    for name, candidate in (("source", image), ("original", original_reference),
                            ("working", working_reference)):
        if not isinstance(candidate, (nib.Nifti1Image, nib.Nifti2Image)) or len(candidate.shape) != 3:
            raise ValueError(f"{name} must be a 3D NIfTI image")
    if image.shape != working_reference.shape or not np.allclose(
        image.affine, working_reference.affine, rtol=0, atol=1e-4
    ):
        raise ValueError("Native export source must share the recorded working T1 grid")
    values = np.asanyarray(image.dataobj)
    if not np.isfinite(values).all():
        raise ValueError("Native export source contains nonfinite values")
    if labels and not np.equal(values, np.rint(values)).all():
        raise ValueError("Native export segmentation must contain integer labels")
    dtype = (values.dtype if np.issubdtype(values.dtype, np.integer) else np.dtype("int32")) if labels else np.dtype("float32")
    if orientation_only:
        orientation = nib.orientations.ornt_transform(
            nib.orientations.io_orientation(working_reference.affine),
            nib.orientations.io_orientation(original_reference.affine),
        )
        values = nib.orientations.apply_orientation(values, orientation)
        if values.shape != original_reference.shape:
            raise ValueError("Orientation-only working grid does not map to the original dimensions")
    elif image.shape != original_reference.shape or not np.allclose(
        image.affine, original_reference.affine, rtol=0, atol=1e-6
    ):
        # Float continuous source avoids integer quantization in scipy when
        # exporting a continuous image that happens to use integer storage.
        source = image if labels else nib.Nifti1Image(values.astype(np.float32), image.affine)
        sampled = resample_from_to(
            source, (original_reference.shape, original_reference.affine),
            order=0 if labels else 1, mode="constant", cval=outside_value,
        )
        values = np.asanyarray(sampled.dataobj)
    return _on_original_reference(original_reference, values, dtype)


def export_base_native_outputs(
    *, reference_path, working_reference_path, segmentation_path,
    inner_path, central_path, outer_path, surface_dir, output_dir,
    orientation_only=False,
):
    """Export final products; retain intermediates and existing mesh topology."""
    reference = nib.load(str(reference_path))
    working = nib.load(str(working_reference_path))
    # Validate both forms while retaining the reference's original storage axes.
    canonical_reference = normalize_image_to_ras(reference, source_name="original output reference")
    normalize_image_to_ras(working, source_name="working output reference")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    paths = {name: root / name for name in NATIVE_OUTPUT_NAMES}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite native BASE export products")
    sources = {"SEG_refine.nii.gz": Path(segmentation_path),
               "inner.nii.gz": Path(inner_path), "central.nii.gz": Path(central_path),
               "outer.nii.gz": Path(outer_path),
               "thickness_levelset.nii.gz": Path(surface_dir) / "thickness_levelset.nii.gz"}
    for name, source in sources.items():
        image = restore_image_to_original_grid(
            nib.load(str(source)), reference, working, orientation_only=orientation_only,
            labels=name == "SEG_refine.nii.gz",
            outside_value=20.0 if name in ("inner.nii.gz", "central.nii.gz", "outer.nii.gz") else 0.0,
        )
        nib.save(image, str(paths[name]))
        saved = nib.load(str(paths[name]))
        if saved.shape != reference.shape or not np.allclose(saved.affine, reference.affine, rtol=0, atol=1e-6):
            raise RuntimeError(f"Native product does not retain the original grid: {name}")
    # Regridded surfaces already use physical RAS mm. Only axis-only working
    # geometry needs its discarded direction cosines restored. This preserves
    # the extracted surface and per-vertex thickness instead of remeshing it.
    transform = (
        canonical_reference.affine @ np.linalg.inv(working.affine)
        if orientation_only else np.eye(4)
    )
    for name in ("inner_surface.vtk", "central_surface.vtk", "outer_surface.vtk", "inner_surface_thickness.vtk"):
        mesh = read_legacy_vtk(Path(surface_dir) / name)
        write_legacy_vtk(transform_surface_ras(mesh, transform), paths[name])
    summary = json.loads((Path(surface_dir) / "thickness_summary.json").read_text())
    summary["native_export"] = {
        "reference": str(Path(reference_path).resolve()),
        "volume_statistics": "retained from working-grid measurement; native volume is sampled from that map",
        "surface_statistics": "retained from the same mesh measurement in physical millimetres",
    }
    paths["thickness_summary.json"].write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    provenance = {
        "schema_version": 1, "command": "base-native-export",
        "original_reference": str(Path(reference_path).resolve()),
        "working_reference": str(Path(working_reference_path).resolve()),
        "original_shape": list(reference.shape), "original_affine": reference.affine.tolist(),
        "volume_axis_codes": list(nib.aff2axcodes(reference.affine)),
        "coordinate_frame": "RAS", "surface_units": "mm",
        "geometry_mapping": "inverse-axis-reorientation" if orientation_only else "physical-affine",
        "label_interpolation": "none" if orientation_only else "nearest-neighbor",
        "continuous_interpolation": "none" if orientation_only else "linear",
        "surface_transform_ras_mm": transform.tolist(),
        "surface_topology": "preserved; native fields sampled without remeshing",
        "sources": {name: str(path.resolve()) for name, path in sources.items()},
        "outputs": {name: str(path.resolve()) for name, path in paths.items()},
    }
    paths["provenance.json"].write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return paths
