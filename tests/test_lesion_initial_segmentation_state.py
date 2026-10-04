from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from pymacruise.lesion_initial_segmentation_state import (
    initialize_lesion_segmentation_state,
    run_loaded_template_first_thinning,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_initial_segmentation_state"
REFERENCE = json.loads((ROOT / "java_reference.json").read_text())
HARNESS = Path(__file__).parent / "java_reference" / "LesionConstructorInitialClassStateReferenceHarness.java"
HARNESS_SHA = "68e93d57f3d31d54a9ea70e745bfaf1666d59b03b5f74c911ba88b77563a0249"
SOURCE_SHA = "63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884"
SOURCE = (
    Path(__file__).parents[1]
    / ".."
    / "toads-cruise-scm-2026-09-26"
    / "VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java"
)
ATLAS = (
    Path(__file__).parents[1]
    / "data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"
)


def _bits(values: np.ndarray) -> list[int]:
    return [int(value) for value in values.astype(np.float32).ravel().view(np.uint32)]


def test_constructor_class_state_matches_pinned_java_harness() -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA
    template = np.fromfunction(lambda x, y, z: (x + 2 * y + z) % 3, (4, 4, 4), dtype=int)
    labels = np.asarray((1, 5, 24), dtype=np.int16)
    source_template = labels[template]
    state = initialize_lesion_segmentation_state(
        source_template, labels
    )

    expected_segmentation = np.asarray(REFERENCE["segmentation"], dtype=np.int8).reshape(4, 4, 4)
    expected_skeleton = np.asarray(REFERENCE["skeleton"], dtype=np.int8).reshape(4, 4, 4)
    expected_relation = np.asarray(REFERENCE["relation_map"], dtype=np.int8).reshape(4, 4, 4)
    np.testing.assert_array_equal(state.segmentation, expected_segmentation)
    np.testing.assert_array_equal(state.skeleton, expected_skeleton)
    np.testing.assert_array_equal(state.relation_map, expected_relation)
    assert state.available.ravel().tolist() == REFERENCE["available"]
    assert _bits(state.ordering) == REFERENCE["ordering_bits"]
    np.testing.assert_array_equal(state.lesion_membership, np.zeros((4, 4, 4), dtype=np.float32))
    assert _bits(state.lesion_membership) == REFERENCE["lesion_membership_bits"]
    assert _bits(state.outlier_membership) == REFERENCE["outlier_membership_bits"]
    assert _bits(state.class_volume) == REFERENCE["class_volume_bits"]


def test_loaded_atlas_template_state_drives_first_thinning_without_manual_state_maps() -> None:
    from pymacruise.lesion_atlas_priors import (
        load_lesion_atlas_priors,
        load_lesion_atlas_shape_priors,
    )

    priors = load_lesion_atlas_priors(str(ATLAS), ("T1", "FLAIR"))
    shapes = load_lesion_atlas_shape_priors(str(ATLAS))
    center = tuple(size // 2 for size in next(iter(shapes.shape_maps.values())).shape)
    slices = tuple(slice(value - 5, value + 6) for value in center)
    shape_patch = np.stack(
        [shapes.shape_maps[item.name][slices] for item in priors.structures], axis=-1
    ).astype(np.float32)
    dominant = np.argmax(shape_patch, axis=-1)
    images = priors.intensity_priors[:, dominant].copy()

    from pymacruise.lesion_initial_atlas_alignment import (
        initialize_lesion_atlas_alignment_from_descriptor,
    )

    alignment = initialize_lesion_atlas_alignment_from_descriptor(
        str(ATLAS),
        images,
        np.zeros(6, dtype=np.float32),
        crop_bounds=((0, 9), (0, 9), (0, 9)),
        pad_border=1,
        image_resolutions_mm=(0.9375,) * 3,
        template_resolutions_mm=(0.9375,) * 3,
        init_align_iterations=0,
        pyramid_levels=2,
    )
    transformed_template = np.asarray([item.label for item in priors.structures])[dominant]
    object_types = tuple(item.topology for item in priors.structures)
    lesion_class = next(i for i, item in enumerate(priors.structures) if item.topology == "lesion")
    result = run_loaded_template_first_thinning(
        alignment,
        transformed_template,
        shape_patch,
        np.zeros(transformed_template.shape, dtype=np.float32),
        object_types=object_types,
        lesion_class=lesion_class,
        volume_ratio=0.0,
        limit=0.0,
    )

    assert result.stage_order[:2] == (
        "atlas.computeTransformedTemplate",
        "MultipleToadSegmentation.constructor:initializeClassState",
    )
    assert result.state.segmentation.shape == (11, 11, 11)
    assert result.state.available.shape == (11, 11, 11, 13)
    assert result.first_thinning.thinning.run.segmentation.shape == (11, 11, 11)
    assert result.first_thinning.stage_order[-2:] == (
        "initThinningLabels",
        "propagateThinning",
    )


def test_constructor_state_rejects_undeclared_transformed_template_labels() -> None:
    with np.testing.assert_raises_regex(ValueError, "undeclared labels"):
        initialize_lesion_segmentation_state(
            np.asarray([[[1, 99]]], dtype=np.int16), (1, 5)
        )
