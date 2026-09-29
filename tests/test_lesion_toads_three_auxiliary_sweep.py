from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import run_lesion_toads_blackhole_outmask_smoothing_sweep


REFERENCE_ROOT = Path(__file__).parent / "data/synthetic/lesion_toads_three_auxiliary_sweep"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
CLASS_ORDERS = {
    "blackhole-last": ("OBJECT", "LESION", "BLACKHOLE", "OUTMASK", "OBJECT"),
    "lesion-last": ("OBJECT", "BLACKHOLE", "OUTMASK", "LESION", "OBJECT"),
}
RELATIONS = ("direct", "prior", "prior", "prior", "direct")
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75, 0.375, 0.625], dtype=np.float32)
ATLAS_PRIORS = np.asarray([0.5, 0.75, 1.25, 0.875, 0.625], dtype=np.float32)
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


def _inputs(shape: tuple[int, int, int], class_types):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, 5), dtype=np.float32)
    data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(5):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    lesion = (0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    black_hole = (0.125 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    outlier = (0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    lesion_data = (0.625 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    black_hole_data = (0.5 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    outmask_data = (0.375 + 0.03125 * x + 0.015625 * y + 0.0078125 * z).astype(np.float32)
    lesion[1, 1, 1] = 0.0
    outlier[1, 1, 1] = 1.0
    if shape[0] > 3 and shape[1] > 2 and shape[2] > 4:
        outlier[3, 2, 4] = 0.0
    shapes[2, 2, 2, class_types.index("BLACKHOLE")] = 0.0
    if shape[2] > 3:
        shapes[1, 2, 3, class_types.index("LESION")] = 0.0
    ventricle = np.where(x == 1, 0, np.where(x == 2, 3, -1)).astype(np.int8)
    gray_matter = np.where(y == 1, 0, 4).astype(np.int8)
    brainstem = np.where(z == 2, 0, 2).astype(np.int8)
    return (
        memberships, lesion, black_hole, outlier, data, lesion_data, black_hole_data, outmask_data,
        shapes, ventricle, gray_matter, brainstem,
    )


def _run(inputs, class_types):
    return run_lesion_toads_blackhole_outmask_smoothing_sweep(
        *inputs[:9],
        SMOOTHING_PRIORS,
        class_types,
        RELATIONS,
        smoothing=0.2,
        shape_priors=ATLAS_PRIORS,
        class_similarities=CLASS_SIMILARITIES,
        lesion_similarities=LESION_SIMILARITIES,
        black_hole_similarities=BLACK_HOLE_SIMILARITIES,
        prior_coefficient=0.65,
        ventricle_distances=inputs[9],
        gray_matter_distances=inputs[10],
        brainstem_distances=inputs[11],
        max_ventricle_distance=5,
        max_gray_matter_distance=8,
        max_brainstem_distance=4,
        gm_flag=True,
        spread=2.5,
    )


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


@pytest.mark.parametrize(
    ("reference_name", "order_name"),
    (("java_reference_blackhole_last.json", "blackhole-last"), ("java_reference_lesion_last.json", "lesion-last")),
)
def test_three_auxiliary_sweep_matches_java_for_both_lesion_class_choices(reference_name, order_name) -> None:
    reference = json.loads((REFERENCE_ROOT / reference_name).read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["class_order"] == order_name
    class_types = CLASS_ORDERS[order_name]
    inputs = _inputs(tuple(reference["shape"]), class_types)
    originals = tuple(array.copy() for array in inputs[:4])
    result = _run(inputs, class_types)
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected["lesion_bits"]
        assert int(np.asarray(result.black_hole_membership[point], dtype=np.float32).view(np.uint32)) == expected["black_hole_bits"]
        assert int(np.asarray(result.outlier_membership[point], dtype=np.float32).view(np.uint32)) == expected["outlier_bits"]
    border = np.ones(originals[0].shape[:3], dtype=bool)
    border[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[border], originals[0][border])
    np.testing.assert_array_equal(result.lesion_membership[border], originals[1][border])
    np.testing.assert_array_equal(result.black_hole_membership[border], originals[2][border])
    np.testing.assert_array_equal(result.outlier_membership[border], originals[3][border])


def test_three_auxiliary_sweep_rejects_duplicate_special_classes() -> None:
    class_types = ("OBJECT", "LESION", "BLACKHOLE", "OUTMASK", "OUTMASK")
    inputs = _inputs((3, 3, 3), class_types)
    with pytest.raises(ValueError, match="exactly one BLACKHOLE and one OUTMASK"):
        run_lesion_toads_blackhole_outmask_smoothing_sweep(
            *inputs[:9], SMOOTHING_PRIORS, class_types, RELATIONS, smoothing=0.2,
        )
