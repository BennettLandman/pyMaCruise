from __future__ import annotations

import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.skeleton import thin_ace, thin_ace_image
from pymacruise.surfaces import extract_levelset_surface
from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics


FIXTURES = Path(__file__).parent / "data" / "synthetic" / "skeleton"


def test_thin_ace_matches_vandy_java_reference_for_synthetic_shapes() -> None:
    spec = json.loads((FIXTURES / "spec.json").read_text(encoding="utf-8"))

    for name, case in spec["cases"].items():
        source_image = nib.load(FIXTURES / case["input"])
        reference_image = nib.load(FIXTURES / case["java_reference"])
        candidate = thin_ace(np.asarray(source_image.dataobj))
        reference = np.asarray(reference_image.dataobj)

        assert nib.aff2axcodes(source_image.affine) == ("R", "A", "S"), name
        assert nib.aff2axcodes(reference_image.affine) == ("R", "A", "S"), name
        assert set(np.unique(candidate)) <= {0, 254}, name
        assert dice_coefficient(reference > 0, candidate > 0) >= 0.97, name
        np.testing.assert_array_equal(candidate, reference, err_msg=name)

        reference_surface = extract_levelset_surface(reference_image, level=127)
        candidate_image = nib.Nifti1Image(candidate.astype(np.uint8), source_image.affine)
        candidate_surface = extract_levelset_surface(candidate_image, level=127)
        metrics = mesh_surface_distance_metrics(
            reference_surface.vertices,
            reference_surface.faces,
            candidate_surface.vertices,
            candidate_surface.faces,
        )
        assert metrics["hausdorff_mm"] <= 1.0, name
        assert metrics["mean_absolute_distance_mm"] <= 0.01, name


def test_thin_ace_reference_nifti_adapter_reports_ras() -> None:
    source = nib.load(FIXTURES / "nested_gm_wm_input.nii.gz")
    result = thin_ace_image(source)
    reference = nib.load(FIXTURES / "nested_gm_wm_thinace_java.nii.gz")

    assert nib.aff2axcodes(result.skeleton.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(result.skeleton.dataobj, reference.dataobj)
