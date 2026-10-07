from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_dura import dura_initial_level_sets, run_dura_removal


FIXTURE = Path(__file__).parent / "data" / "synthetic" / "dura_attached_no_dura"
HARNESS = Path(__file__).parent / "java_reference" / "DuraAttachedNoDuraReferenceHarness.java"
EXPECTED_HARNESS_SHA256 = "e1f3030ea31884b1b3ac3aae9eaf42216876bd55be9f1ec88358f00d391dc7aa"
EXPECTED_SOURCE_HASHES = {
    "MedicAlgorithmRemoveDura09": "1ed180daf4a4334c995ad310828ec91cfd5593413eec17ab15404a421d7805a4",
    "AlgorithmTopologyCorrection": "9dd7cd0db8c83e47a21059fe5acc52c7737cb0001222d6019619bf62526159a8",
    "GenericTGDM": "70f5b21541e8a6c6139f1bd227c4ece12e8818db717f2a78ef50d69fe89f08d1",
    "DynamicPressureTGDM": "53b872a5828539bb51d70db0c048d1fd8162f6b7a6b63176266ef7e2e796b6fc",
    "DynamicPressureCortexProfile": "f85912d621198de297e65d6b9205505b3ac2de72e67e6840f3968c58a850f505",
    "DistanceField": "92294be2db74351446c3e5898fc35b58497f89ccb38772d4d565885f3f0c35b7",
}


def _load_reference(case_dir: Path, filename: str) -> np.ndarray:
    return np.asarray(nib.load(str(case_dir / "java_reference" / filename)).dataobj)


def _assert_ras_output(image: nib.spatialimages.SpatialImage, affine: np.ndarray) -> None:
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(image.affine, affine)


def test_attached_and_clean_dura_cases_match_selected_java_wrapper_sequence() -> None:
    provenance = json.loads((FIXTURE / "provenance.json").read_text())
    assert provenance["scope"] == "synthetic selected-wrapper-sequence comparison"
    assert provenance["source_sha256"] == EXPECTED_SOURCE_HASHES
    assert provenance["harness_sha256"] == EXPECTED_HARNESS_SHA256
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == EXPECTED_HARNESS_SHA256

    expected_fixture_hashes = json.loads((FIXTURE / "fixture_hashes.json").read_text())
    actual_fixture_hashes = {
        str(path.relative_to(FIXTURE)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(FIXTURE.rglob("*.nii.gz"))
    }
    assert actual_fixture_hashes == expected_fixture_hashes

    for case_name in ("attached_dura", "no_dura"):
        case_dir = FIXTURE / case_name
        t1, wm, gm, brain_mask = (
            nib.load(str(case_dir / f"{name}.nii.gz"))
            for name in ("t1", "wm", "gm", "brain_mask")
        )
        result = run_dura_removal(t1, wm, gm, brain_mask, record_passes=True)
        expected_wm_topology = _load_reference(case_dir, "wm_topology.nii.gz")
        expected_gm_topology = _load_reference(case_dir, "gm_topology.nii.gz")
        expected_levels = dura_initial_level_sets(
            expected_wm_topology.astype(np.float32),
            expected_gm_topology.astype(np.float32),
        )

        np.testing.assert_array_equal(result.wm_topology.dataobj, expected_wm_topology)
        np.testing.assert_array_equal(result.gm_topology.dataobj, expected_gm_topology)
        _assert_ras_output(result.wm_topology, t1.affine)
        _assert_ras_output(result.gm_topology, t1.affine)
        np.testing.assert_array_equal(
            expected_levels.inner_signed_distance,
            _load_reference(case_dir, "inner_distance.nii.gz"),
        )
        np.testing.assert_array_equal(
            expected_levels.outer_signed_distance,
            _load_reference(case_dir, "outer_distance.nii.gz"),
        )

        trace = json.loads((case_dir / "java_reference" / "trace.json").read_text())
        assert result.evolution.completed_outer_passes == len(trace["trace_steps"]) == 5
        assert result.evolution.completed_steps == sum(trace["trace_steps"]) == 150
        assert len(result.evolution.pass_level_sets) == len(trace["trace_levelset_sha256"]) == 5
        for pass_index, actual_pass in enumerate(result.evolution.pass_level_sets, start=1):
            expected_pass = _load_reference(
                case_dir,
                f"tgdm_pass_{pass_index:02d}.nii.gz",
            )
            java_trace_hash = hashlib.sha256(
                np.ascontiguousarray(expected_pass, dtype="<f4").tobytes()
            ).hexdigest()
            assert java_trace_hash == trace["trace_levelset_sha256"][pass_index - 1]
            np.testing.assert_allclose(actual_pass, expected_pass, rtol=0.0, atol=1e-6)

        expected_final_level_set = _load_reference(case_dir, "level_set.nii.gz")
        np.testing.assert_allclose(
            result.evolution.level_set,
            expected_final_level_set,
            rtol=0.0,
            atol=1e-6,
        )
        output_references = (
            (result.stripped_t1, "stripped_t1.nii.gz"),
            (result.corrected_white_matter, "corrected_wm.nii.gz"),
            (result.corrected_gray_matter, "corrected_gm.nii.gz"),
            (result.brain_mask_level_set, "brain_mask_level_set.nii.gz"),
            (result.brain_mask, "brain_mask.nii.gz"),
        )
        for actual, reference_name in output_references:
            expected = _load_reference(case_dir, reference_name)
            if reference_name == "brain_mask_level_set.nii.gz":
                np.testing.assert_allclose(actual.dataobj, expected, rtol=0.0, atol=1e-6)
            else:
                np.testing.assert_array_equal(actual.dataobj, expected)
            _assert_ras_output(actual, t1.affine)
