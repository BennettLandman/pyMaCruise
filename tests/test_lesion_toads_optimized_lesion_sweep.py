from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pymacruise.lesion_membership import run_lesion_toads_lesion_smoothing_sweep


FIXTURE_ROOT = Path(__file__).parent / "data/synthetic/lesion_toads_lesion_sweep"
SOURCE_SHA256 = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
SMOOTHING_PRIORS = np.asarray([0.25, 0.5, 0.75], dtype=np.float32)


def _inputs(shape: tuple[int, int, int], lesion_class: int):
    x, y, z = np.indices(shape, dtype=np.float32)
    memberships = np.empty((*shape, 3), dtype=np.float32)
    data = np.empty_like(memberships)
    shapes = np.empty_like(memberships)
    for k in range(3):
        memberships[..., k] = 0.5 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        data[..., k] = 0.75 + 0.125 * x + 0.0625 * y + 0.03125 * z + 0.0625 * k
        shapes[..., k] = 0.25 + 0.0625 * x + 0.03125 * y + 0.015625 * z + 0.0625 * k
    lesion = (0.25 + 0.125 * x + 0.0625 * y + 0.03125 * z).astype(np.float32)
    lesion_data = (0.75 + 0.0625 * x + 0.03125 * y + 0.015625 * z).astype(np.float32)
    shapes[2, 2, 2, lesion_class] = 0.0
    return memberships, lesion, data, lesion_data, shapes


def _bits(values: np.ndarray) -> list[int]:
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


@pytest.mark.parametrize(
    ("reference_name", "class_order", "lesion_class"),
    (
        ("java_reference_optimized_middle.json", "optimized-middle", 1),
        ("java_reference_optimized_last.json", "optimized-last", 2),
    ),
)
def test_optimized_lesion_volume_sweep_matches_selected_java(
    reference_name: str, class_order: str, lesion_class: int
) -> None:
    reference = json.loads((FIXTURE_ROOT / reference_name).read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["class_order"] == class_order
    assert reference["optimized_lesion"] is True
    class_types = ["OBJECT", "OBJECT", "OBJECT"]
    class_types[lesion_class] = "OPTIMIZED_LESION"
    relation_types = ["direct", "direct", "direct"]
    relation_types[lesion_class] = "prior"
    inputs = _inputs(tuple(reference["shape"]), lesion_class)
    original_classes, original_lesion = inputs[0].copy(), inputs[1].copy()
    result = run_lesion_toads_lesion_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        class_types,
        relation_types,
        smoothing=float(np.asarray(reference["smoothing_bits"], dtype=np.uint32).view(np.float32)),
    )
    for expected in reference["voxels"]:
        point = tuple(expected["coordinate"])
        assert _bits(result.class_memberships[point]) == expected["class_bits"]
        assert int(np.asarray(result.lesion_membership[point], dtype=np.float32).view(np.uint32)) == expected[
            "lesion_bits"
        ]
    border = np.ones(original_classes.shape[:3], dtype=bool)
    border[1:-1, 1:-1, 1:-1] = False
    np.testing.assert_array_equal(result.class_memberships[border], original_classes[border])
    np.testing.assert_array_equal(result.lesion_membership[border], original_lesion[border])


def test_object_lesion_sweep_rejects_duplicate_lesion_aliases() -> None:
    inputs = _inputs((3, 3, 3), lesion_class=1)
    with pytest.raises(ValueError, match="exactly one LESION or OPTIMIZED_LESION"):
        run_lesion_toads_lesion_smoothing_sweep(
            *inputs,
            SMOOTHING_PRIORS,
            ("OBJECT", "LESION", "OPTIMIZED_LESION"),
            ("direct", "prior", "direct"),
            smoothing=0.2,
        )


def test_optimized_lesion_matches_lesion_alias_through_priors_and_distances() -> None:
    lesion_class = 2
    inputs = _inputs((5, 4, 6), lesion_class)
    shape_priors = np.asarray([0.5, 0.75, 1.25], dtype=np.float32)
    class_similarities = np.asarray(
        [[0.0, 0.2, 0.7], [0.4, 0.0, 0.3], [0.6, 0.1, 0.0]], dtype=np.float32
    )
    lesion_similarities = np.asarray([0.25, 0.0, 0.5], dtype=np.float32)
    x, y, z = np.indices((5, 4, 6))
    distances = {
        "ventricle_distances": np.where(x == 1, 0, -1).astype(np.int8),
        "gray_matter_distances": np.where(y == 1, 0, -1).astype(np.int8),
        "brainstem_distances": np.where(z == 2, 2, -1).astype(np.int8),
        "max_ventricle_distance": 5,
        "max_gray_matter_distance": 8,
        "max_brainstem_distance": 4,
        "gm_flag": True,
        "spread": 2.5,
    }
    common = {
        "smoothing": 0.2,
        "shape_priors": shape_priors,
        "class_similarities": class_similarities,
        "lesion_similarities": lesion_similarities,
        "prior_coefficient": 0.65,
        **distances,
    }
    lesion_result = run_lesion_toads_lesion_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        ("OBJECT", "OBJECT", "LESION"),
        ("direct", "direct", "prior"),
        **common,
    )
    optimized_result = run_lesion_toads_lesion_smoothing_sweep(
        *inputs,
        SMOOTHING_PRIORS,
        ("OBJECT", "OBJECT", "OPTIMIZED_LESION"),
        ("direct", "direct", "prior"),
        **common,
    )
    np.testing.assert_array_equal(optimized_result.class_memberships, lesion_result.class_memberships)
    np.testing.assert_array_equal(optimized_result.lesion_membership, lesion_result.lesion_membership)
