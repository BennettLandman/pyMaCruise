from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_membership_pass import run_lesion_toads_membership_map_pass


FIXTURE = Path(__file__).parent / "data/synthetic/lesion_membership_map_pass"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
INPUT_SHA256 = "59c0f5976256f1488616377a0b3b08f55af71d57e59a8e0cdef0f2b57d31c791"
CLASS_TYPES = ("OBJECT", "LESION", "BLACKHOLE", "OUTMASK", "OBJECT")
RELATION_TYPES = ("direct", "prior", "prior", "prior", "direct")
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75, 0.375, 0.625], dtype=np.float32)
SHAPE_PRIORS = np.asarray([0.5, 0.75, 1.25, 0.875, 0.625], dtype=np.float32)
CLASS_SIMILARITIES = np.asarray(
    [
        [0.0, 0.2, 0.7, 0.3, 0.1],
        [0.4, 0.0, 0.3, 0.6, 0.2],
        [0.6, 0.1, 0.0, 0.2, 0.5],
        [0.2, 0.5, 0.4, 0.0, 0.3],
        [0.3, 0.4, 0.1, 0.7, 0.0],
    ],
    dtype=np.float32,
)
LESION_SIMILARITIES = np.asarray([0.25, 0.0, 0.5, 0.375, 0.625], dtype=np.float32)
BLACK_HOLE_SIMILARITIES = np.asarray([0.125, 0.5, 0.0, 0.25, 0.375], dtype=np.float32)


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def test_folded_ras_membership_maps_match_selected_java_harness() -> None:
    path = FIXTURE / "inputs_ras.nii.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == INPUT_SHA256
    image = nib.load(path)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert image.shape == (7, 8, 9, 24)
    volume = np.asarray(image.dataobj, dtype=np.float32)
    reference = json.loads((FIXTURE / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["shape"] == list(image.shape[:3])
    assert reference["class_order"] == list(CLASS_TYPES)
    assert len(reference["voxels"]) == (7 - 2) * (8 - 2) * (9 - 2)

    initial_memberships = volume[..., :5].copy()
    initial_lesion = volume[..., 5].copy()
    initial_black_hole = volume[..., 6].copy()
    initial_outlier = volume[..., 7].copy()
    original = (initial_memberships.copy(), initial_lesion.copy(), initial_black_hole.copy(), initial_outlier.copy())
    result = run_lesion_toads_membership_map_pass(
        initial_memberships,
        initial_lesion,
        initial_black_hole,
        initial_outlier,
        volume[..., 8:13],
        volume[..., 13],
        volume[..., 14],
        volume[..., 15],
        volume[..., 16:21],
        SMOOTHING_PRIORS,
        CLASS_TYPES,
        RELATION_TYPES,
        smoothing=0.2,
        shape_priors=SHAPE_PRIORS,
        class_similarities=CLASS_SIMILARITIES,
        lesion_similarities=LESION_SIMILARITIES,
        black_hole_similarities=BLACK_HOLE_SIMILARITIES,
        prior_coefficient=0.65,
        ventricle_distances=volume[..., 21].astype(np.int8),
        gray_matter_distances=volume[..., 22].astype(np.int8),
        brainstem_distances=volume[..., 23].astype(np.int8),
        max_ventricle_distance=5,
        max_gray_matter_distance=8,
        max_brainstem_distance=4,
        gm_flag=True,
        spread=2.5,
    )

    for expected in reference["voxels"]:
        coordinate = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[coordinate]) == expected["class_bits"], coordinate
        assert int(result.lesion_membership[coordinate].view(np.uint32)) == expected["lesion_bits"], coordinate
        assert int(result.black_hole_membership[coordinate].view(np.uint32)) == expected["black_hole_bits"], coordinate
        assert int(result.outlier_membership[coordinate].view(np.uint32)) == expected["outlier_bits"], coordinate

    boundary = np.ones(initial_memberships.shape[:3], dtype=bool)
    boundary[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[boundary], original[0][boundary])
    np.testing.assert_array_equal(result.lesion_membership[boundary], original[1][boundary])
    np.testing.assert_array_equal(result.black_hole_membership[boundary], original[2][boundary])
    np.testing.assert_array_equal(result.outlier_membership[boundary], original[3][boundary])
