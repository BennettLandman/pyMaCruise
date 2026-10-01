from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.image_io import load_nifti_ras
from pymacruise.lesion_registration_setup_handoff import lesion_registration_setup_handoff


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data" / "synthetic" / "lesion_registration_setup_handoff"
HARNESS = ROOT / "java_reference" / "LesionToadsRegistrationSetupHandoffReferenceHarness.java"
REFERENCE = FIXTURE / "java_reference.json"
SOURCE_EXCERPTS = {
    "wrapper_defaults": FIXTURE / "wrapper_defaults_source_excerpt.java.txt",
    "wrapper_call": FIXTURE / "wrapper_registration_call_source_excerpt.java.txt",
    "membership_buffer": FIXTURE / "membership_buffer_source_excerpt.java.txt",
    "atlas_setup": FIXTURE / "atlas_registration_setup_source_excerpt.java.txt",
}
INPUT_NIFTI_SHA256 = "a9c12cae756929a8215883cb9afd7999da033ee815c98e8cdc4abcbef26f0fc4"
SOURCE_HASHES = {
    "wrapper": "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb",
    "segmentation": "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99",
    "atlas": "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e",
}


def test_registration_setup_handoff_matches_selected_wrapper_arguments() -> None:
    ref = json.loads(REFERENCE.read_text())
    assert ref["source_sha256"] == SOURCE_HASHES
    assert ref["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    for name, path in SOURCE_EXCERPTS.items():
        assert ref["source_excerpt_sha256"][name] == hashlib.sha256(path.read_bytes()).hexdigest()
    wrapper_text = SOURCE_EXCERPTS["wrapper_call"].read_text()
    wrapper_calls = [
        "segmentation.initCentroids(centroids);",
        "segmentation.setLesionCentroid(",
        "segmentation.setBlackHoleCentroid(",
        "segmentation.setGMflag();",
        "atlas.updateShapeRegistration(segmentation.getMemberships(), registrationMode, alignIter, 1, 1.0f);",
    ]
    assert [wrapper_text.index(call) for call in wrapper_calls] == sorted(
        wrapper_text.index(call) for call in wrapper_calls
    )
    assert 'registrationMode = "rigid"' in SOURCE_EXCERPTS["wrapper_defaults"].read_text()
    assert "alignIter = 1" in SOURCE_EXCERPTS["wrapper_defaults"].read_text()
    membership_text = SOURCE_EXCERPTS["membership_buffer"].read_text()
    assert "mems = new float[nx][ny][nz][classes];" in membership_text
    assert "getMemberships() { return mems; }" in membership_text
    atlas_text = SOURCE_EXCERPTS["atlas_setup"].read_text()
    atlas_assignments = [
        "mems = mems_;",
        "itMax = iter_;",
        "levels = lvl_;",
        "subsample = 2;",
        "couplingFactor = cf_;",
    ]
    assert [atlas_text.index(value) for value in atlas_assignments] == sorted(
        atlas_text.index(value) for value in atlas_assignments
    )
    input_path = FIXTURE / "initial_memberships.nii.gz"
    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == INPUT_NIFTI_SHA256
    image = load_nifti_ras(input_path)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (0.8, 1.0, 1.2), atol=1e-6)
    memberships = np.asarray(image.dataobj, dtype=np.float32)
    assert memberships.shape == tuple(ref["shape"])
    assert not np.any(memberships)
    before = memberships.copy()

    result = lesion_registration_setup_handoff(memberships)

    assert result.memberships is memberships
    assert result.transform_type == ref["mode"] == "rigid"
    assert result.iterations == ref["iterations"] == 1
    assert result.levels == ref["levels"] == 1
    assert result.subsample == ref["subsample"] == 2
    assert result.coupling_factor.view(np.uint32) == ref["coupling_factor_bits"]
    assert result.call_order == tuple(ref["call_order"])
    np.testing.assert_array_equal(memberships, before)
