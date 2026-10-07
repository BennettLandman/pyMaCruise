from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_dura import dura_initial_level_sets, run_dura_removal


FIXTURE = Path(__file__).parent / "data/synthetic/dura_attached_no_dura_41"


def _image(case: str, folder: str, name: str) -> nib.spatialimages.SpatialImage:
    path = FIXTURE / case / folder / f"{name}.nii.gz"
    image = nib.load(path)
    assert image.shape == (41, 41, 41)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return image


def _reference(case: str, name: str) -> np.ndarray:
    return np.asarray(_image(case, "java_reference", name).dataobj)


def test_attached_and_clean_dura_41_match_selected_java_wrapper_sequence() -> None:
    provenance = json.loads((FIXTURE / "provenance.json").read_text())
    assert provenance["scope"] == "synthetic selected-wrapper-sequence comparison"
    assert provenance["dimensions_xyz"] == [41, 41, 41]
    assert provenance["harness_sha256"] == "e1f3030ea31884b1b3ac3aae9eaf42216876bd55be9f1ec88358f00d391dc7aa"
    for line in (FIXTURE / "SHA256SUMS").read_text().splitlines():
        digest, relative = line.split("  ", 1)
        assert hashlib.sha256((FIXTURE / relative).read_bytes()).hexdigest() == digest

    for case in ("attached_dura", "no_dura"):
        t1, wm, gm, mask = (
            _image(case, "inputs", name) for name in ("t1", "wm", "gm", "paint")
        )
        result = run_dura_removal(t1, wm, gm, mask, record_passes=True)

        np.testing.assert_array_equal(result.wm_topology.dataobj, _reference(case, "wm_tc"))
        np.testing.assert_array_equal(result.gm_topology.dataobj, _reference(case, "gm_tc"))
        initial = dura_initial_level_sets(
            np.asarray(result.wm_topology.dataobj, dtype=np.float32),
            np.asarray(result.gm_topology.dataobj, dtype=np.float32),
        )
        np.testing.assert_allclose(initial.inner_signed_distance, _reference(case, "inner_distance"), rtol=0.0, atol=2e-3)
        np.testing.assert_allclose(initial.outer_signed_distance, _reference(case, "outer_distance"), rtol=0.0, atol=2e-3)

        for index, values in enumerate(result.evolution.pass_level_sets, start=1):
            expected_pass = _reference(case, f"tgdm_pass_{index:02d}")
            assert np.array_equal(np.signbit(values), np.signbit(expected_pass))
        np.testing.assert_allclose(
            result.evolution.level_set,
            _reference(case, "level_set"),
            rtol=0.0,
            atol=1e-3,
        )
        for image, name in (
            (result.stripped_t1, "stripped_t1"),
            (result.corrected_white_matter, "corrected_wm"),
            (result.corrected_gray_matter, "corrected_gm"),
            (result.brain_mask_level_set, "brain_mask_level_set"),
            (result.brain_mask, "brain_mask"),
        ):
            expected = _reference(case, name)
            if name == "brain_mask_level_set":
                np.testing.assert_allclose(image.dataobj, expected, rtol=0.0, atol=1e-3)
            else:
                np.testing.assert_array_equal(image.dataobj, expected)
            assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
            np.testing.assert_array_equal(image.affine, t1.affine)
