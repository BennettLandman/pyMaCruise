from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_similarity import compute_lesion_toads_similarity


REFERENCE = Path(__file__).parent / "data/synthetic/lesion_toads_similarity/java_reference.json"
LESION_TYPES = {"lesion", "blackhole", "optimized_lesion"}


def _from_bits(values: list[int]) -> np.ndarray:
    return np.asarray(values, dtype=np.uint32).view(np.float32)


def _to_java_bits(values: np.ndarray) -> list[int]:
    output = []
    for value in np.asarray(values, dtype=np.float32).ravel():
        if np.isnan(value):
            output.append(0x7FC00000)  # Float.floatToIntBits canonical NaN
        else:
            output.append(int(np.asarray(value, dtype=np.float32).view(np.uint32)))
    return output


def test_centroid_similarity_matches_selected_java_float32_formula() -> None:
    reference = json.loads(REFERENCE.read_text())
    assert reference["source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert [case["name"] for case in reference["cases"]] == [
        "mixed_t1",
        "no_t1",
        "no_lesion",
        "zero_std",
    ]

    for case in reference["cases"]:
        object_types = tuple(case["types"])
        initial_black_hole_similarity = np.ones(len(object_types), dtype=np.float32)
        lesion_class = max(
            (i for i, name in enumerate(object_types) if name in LESION_TYPES),
            default=-1,
        )
        if lesion_class >= 0:
            initial_black_hole_similarity[lesion_class] = np.float32(0.0)
        result = compute_lesion_toads_similarity(
            _from_bits([v for channel in case["centroid"] for v in channel]).reshape(2, 4),
            _from_bits([v for channel in case["stddev"] for v in channel]).reshape(2, 4),
            _from_bits(case["lesion_centroid"]),
            _from_bits(case["lesion_stddev"]),
            _from_bits(case["channel_weight"]),
            object_types,
            ("Sulcal-CSF", "Cerebrum-GM", "Cerebrum-WM", "Cerebellum"),
            prior_scale=np.float32(0.15),
            t1_channel=int(case["t1_channel"]),
            black_hole_similarity=initial_black_hole_similarity,
        )
        assert _to_java_bits(result.class_similarity) == case["class_similarity"], case["name"]
        assert _to_java_bits(result.lesion_similarity) == case["lesion_similarity"], case["name"]
        assert _to_java_bits(result.black_hole_similarity) == case["black_hole_similarity"], case["name"]

    mixed = reference["cases"][0]
    no_t1 = reference["cases"][1]
    # Sulcal-CSF pairs use only T1 when available and all channels otherwise.
    assert mixed["class_similarity"] != no_t1["class_similarity"]
    assert mixed["lesion_similarity"][0] != 0
    assert mixed["lesion_similarity"][1:] == [0, 0, 0]
    assert mixed["black_hole_similarity"] == [1065353216, 1065353216, 1065353216, 0]
    assert any(value == 0x7FC00000 for value in reference["cases"][3]["class_similarity"])


def test_explicit_finite_zero_variance_mode_uses_limiting_similarity() -> None:
    centroids = np.asarray(((0.0, 0.0, 1.0), (0.0, 0.0, 1.0)), dtype=np.float32)
    deviations = np.zeros_like(centroids)
    lesion_centroid = np.asarray((0.0, 0.0), dtype=np.float32)
    lesion_deviation = np.zeros(2, dtype=np.float32)
    weights = np.ones(2, dtype=np.float32)
    types = ("object", "object", "object")
    names = ("a", "b", "c")
    black_hole = np.ones(3, dtype=np.float32)

    source = compute_lesion_toads_similarity(
        centroids, deviations, lesion_centroid, lesion_deviation, weights,
        types, names, prior_scale=0.1, black_hole_similarity=black_hole,
    )
    assert np.isnan(source.class_similarity[0, 1])
    finite = compute_lesion_toads_similarity(
        centroids, deviations, lesion_centroid, lesion_deviation, weights,
        types, names, prior_scale=0.1, black_hole_similarity=black_hole,
        finite_zero_variance=True,
    )
    assert np.isfinite(finite.class_similarity).all()
    assert np.isfinite(finite.lesion_similarity).all()
    assert finite.class_similarity[0, 1] == np.float32(1.0)
    assert finite.class_similarity[0, 2] == np.float32(0.0)
    assert finite.lesion_similarity[0] == np.float32(1.0)
    assert finite.lesion_similarity[2] == np.float32(0.0)
