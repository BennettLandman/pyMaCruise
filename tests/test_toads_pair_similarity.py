from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.toads_math import toads_optimized_centroid_similarity
from pymacruise.toads_pair_similarity import toads_class_pair_similarity


REFERENCE = Path(__file__).parent / "data/synthetic/toads_pair_similarity/java_reference.json"
SOURCE_SHA256 = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
CENTERS = np.asarray([[0.17, 0.73, 0.41, 0.17], [0.81, 0.29, 0.66, 0.12]], dtype=np.float32)
STDDEV = np.asarray([[0.23, 0.31, 0.47, 0.19], [0.37, 0.28, 0.21, 0.43]], dtype=np.float32)
LESION_CENTER = np.asarray([0.84, 0.08], dtype=np.float32)
LESION_STDDEV = np.asarray([0.13, 0.34], dtype=np.float32)
WEIGHTS = np.asarray([0.63, 1.27], dtype=np.float32)
TYPES = np.asarray([1, 4, 4, 1], dtype=np.int8)  # VandyCRUISE OBJECT=1, LESION=4
PRIOR_SCALE = np.asarray(0.17, dtype=np.float32).item()


def _bits(values: np.ndarray) -> list[int]:
    output = []
    for value in np.asarray(values, dtype=np.float32).ravel():
        if np.isnan(value):
            output.append(0x7FC00000)
        else:
            output.append(int(np.asarray(value, dtype=np.float32).view(np.uint32)))
    return output


def _inputs(zero_deviation: bool) -> np.ndarray:
    deviations = STDDEV.copy()
    if zero_deviation:
        deviations[0, 0] = np.float32(0.0)
        deviations[0, 3] = np.float32(0.0)
        deviations[1, 3] = np.float32(0.0)
    return deviations


def test_selected_java_pair_dispatch_matches_float_bits() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == SOURCE_SHA256
    assert reference["source"].endswith("similarity prelude excerpt")
    assert [case["name"] for case in reference["cases"]] == [
        "asymmetric_mixed_pairs",
        "zero_deviation",
    ]

    for case in reference["cases"]:
        deviations = _inputs(case["name"] == "zero_deviation")
        result = toads_class_pair_similarity(
            CENTERS,
            deviations,
            WEIGHTS,
            LESION_CENTER,
            LESION_STDDEV,
            TYPES,
            prior_scale=PRIOR_SCALE,
        )
        assert _bits(result.similarity) == case["similarity_bits"]
        assert _bits(result.centroid_distance) == case["centroid_distance_bits"]
        assert _bits(result.lesion_distance) == case["lesion_distance_bits"]
        if case["name"] == "zero_deviation":
            assert any(value == 0x7FC00000 for value in case["similarity_bits"])


def test_mixed_lesion_distance_is_computed_but_not_used_for_coefficient() -> None:
    result = toads_class_pair_similarity(
        CENTERS, STDDEV, WEIGHTS, LESION_CENTER, LESION_STDDEV, TYPES, PRIOR_SCALE
    )
    ordinary = toads_optimized_centroid_similarity(CENTERS, STDDEV, WEIGHTS, PRIOR_SCALE)
    # OBJECT -> LESION uses lesionCentroid against object k, then retains
    # centroid[k]-centroid[m] for its similarity coefficient.
    assert result.lesion_distance[0, 1] != result.centroid_distance[0, 1]
    assert result.similarity[0, 1] == ordinary[0, 1]
    # LESION -> OBJECT dispatch is asymmetric in the lesion-distance record.
    assert result.lesion_distance[1, 0] != result.centroid_distance[1, 0]
    assert result.similarity[1, 0] == ordinary[1, 0]
    # Two lesion classes force ordinary centroid distance to zero and hence
    # coefficient one at a regular priorScale.
    assert result.centroid_distance[1, 2] == np.float32(0.0)
    assert result.similarity[1, 2] == np.float32(1.0)


def test_object_type_names_and_validation() -> None:
    named = toads_class_pair_similarity(
        CENTERS, STDDEV, WEIGHTS, LESION_CENTER, LESION_STDDEV,
        ["OBJECT", "LESION", "LESION", "OBJECT"], PRIOR_SCALE,
    )
    coded = toads_class_pair_similarity(
        CENTERS, STDDEV, WEIGHTS, LESION_CENTER, LESION_STDDEV, TYPES, PRIOR_SCALE
    )
    np.testing.assert_array_equal(named.similarity, coded.similarity)
    try:
        toads_class_pair_similarity(
            CENTERS, STDDEV, WEIGHTS, LESION_CENTER, LESION_STDDEV,
            [1, 4, 7, 1], PRIOR_SCALE,
        )
    except ValueError as exc:
        assert "unsupported object type" in str(exc)
    else:
        raise AssertionError("unsupported source class types must fail explicitly")
