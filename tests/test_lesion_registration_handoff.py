from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_registration_handoff import register_lesion_membership_handoff
from test_lesion_statistics_membership_bridge import _run


ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_registration_handoff"
HARNESS = Path(__file__).parent / "java_reference" / "LesionRegistrationHandoffReferenceHarness.java"
JAVA_REFERENCE = ROOT / "java_reference.json"
SOURCE_ROOT = (
    Path(__file__).parents[2]
    / "toads-cruise-scm-2026-09-26"
    / "VandyCRUISE"
    / "src"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
)
ATLAS_SOURCE_SHA = "718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e"
WRAPPER_SOURCE_SHA = "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
INITIAL_TRANSFORM = (0.03, -0.02, 0.01, 0.1, -0.1, 0.05)
RESOLUTIONS = (0.9, 1.2, 1.5)


def _inputs() -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    image = nib.load(str(ROOT / "registration_handoff_ras.nii.gz"))
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert frames.shape == (6, 5, 4, 12)
    memberships_no_field = frames[..., :4]
    memberships_with_field = frames[..., 4:8]
    shapes = tuple(frames[..., 8 + index] for index in range(4))
    return np.stack((memberships_no_field, memberships_with_field)), shapes


def _assert_f32_equal(actual: float, expected: float, label: str) -> None:
    assert int(np.asarray(actual, dtype=np.float32).view(np.uint32)) == (
        int(np.asarray(expected, dtype=np.float32).view(np.uint32)) & 0xFFFFFFFF
    ), label


def test_membership_to_registration_wrapper_matches_java_trace_and_bounds() -> None:
    reference = json.loads(JAVA_REFERENCE.read_text())
    assert reference["atlas_source_sha256"] == ATLAS_SOURCE_SHA
    assert reference["wrapper_source_sha256"] == WRAPPER_SOURCE_SHA
    assert reference["harness_sha256"] == hashlib.sha256(HARNESS.read_bytes()).hexdigest()
    saved_maps, shapes = _inputs()

    for index, (name, use_field) in enumerate((("no_field", False), ("with_field", True))):
        bridge = _run(use_field)
        updated = bridge.membership.class_memberships
        np.testing.assert_array_equal(updated, saved_maps[index])
        before_memberships = updated.copy()
        before_shapes = tuple(shape.copy() for shape in shapes)
        actual = register_lesion_membership_handoff(
            bridge,
            shapes,
            INITIAL_TRANSFORM,
            image_resolutions_mm=RESOLUTIONS,
            template_resolutions_mm=RESOLUTIONS,
        )
        expected = reference["cases"][index]
        optimization = actual.optimization
        assert optimization is not None
        assert optimization.stop_reason == expected["stop_reason"]
        assert len(optimization.iterations) == len(expected["iterations"])
        np.testing.assert_array_equal(actual.transform, np.asarray(expected["transform"], dtype=np.float32))
        _assert_f32_equal(optimization.energy, expected["energy"], f"{name} final energy")
        _assert_f32_equal(optimization.lambda_value, expected["lambda"], f"{name} lambda")

        for step_index, (step, expected_step) in enumerate(
            zip(optimization.iterations, expected["iterations"])
        ):
            assert len(step.trial_energies) == len(expected_step["trials"])
            assert step.result.accepted == expected_step["accepted"]
            assert step.result.rejected_trials == expected_step["rejected_trials"]
            np.testing.assert_array_equal(
                step.result.transform,
                np.asarray(expected_step["transform"], dtype=np.float32),
                err_msg=f"{name} step {step_index} accepted transform",
            )
            _assert_f32_equal(
                step.result.energy, expected_step["energy"], f"{name} step {step_index} energy"
            )
            _assert_f32_equal(
                step.result.relative_change,
                expected_step["relative_change"],
                f"{name} step {step_index} relative change",
            )
            _assert_f32_equal(
                step.result.lambda_value,
                expected_step["lambda"],
                f"{name} step {step_index} lambda",
            )
            for trial_index, expected_trial in enumerate(expected_step["trials"]):
                np.testing.assert_array_equal(
                    step.trial_transforms[trial_index],
                    np.asarray(expected_trial["transform"], dtype=np.float32),
                    err_msg=f"{name} step {step_index} trial {trial_index} transform",
                )
                _assert_f32_equal(
                    step.trial_energies[trial_index],
                    expected_trial["energy"],
                    f"{name} step {step_index} trial {trial_index} energy",
                )

        np.testing.assert_array_equal(actual.rotation, np.asarray(expected["rotation"], dtype=np.float32))
        for class_index, (bounds, expected_bounds) in enumerate(
            zip(actual.shape_bounds, expected["bounds"])
        ):
            assert bounds.lower_xyz == tuple(expected_bounds["lower"]), class_index
            assert bounds.upper_xyz == tuple(expected_bounds["upper"]), class_index
            assert bounds.support_count == expected_bounds["support_count"], class_index
        np.testing.assert_array_equal(updated, before_memberships)
        for shape, before in zip(shapes, before_shapes):
            np.testing.assert_array_equal(shape, before)

    no_field = register_lesion_membership_handoff(
        _run(False), shapes, INITIAL_TRANSFORM,
        image_resolutions_mm=RESOLUTIONS,
        template_resolutions_mm=RESOLUTIONS,
    )
    with_field = register_lesion_membership_handoff(
        _run(True), shapes, INITIAL_TRANSFORM,
        image_resolutions_mm=RESOLUTIONS,
        template_resolutions_mm=RESOLUTIONS,
    )
    assert not np.array_equal(no_field.transform, with_field.transform)
    assert no_field.shape_bounds[2].support_count > 0  # boundary support
    assert no_field.shape_bounds[3].support_count == 0  # empty support sentinel
    assert no_field.shape_bounds[3].lower_xyz == (6, 5, 4)
    assert no_field.shape_bounds[3].upper_xyz == (0, 0, 0)


def test_disabled_alignment_preserves_seed_transform_and_refreshes_rotation_bounds() -> None:
    _, shapes = _inputs()
    bridge = _run(False)
    memberships_before = bridge.membership.class_memberships.copy()
    result = register_lesion_membership_handoff(
        bridge,
        shapes,
        INITIAL_TRANSFORM,
        image_resolutions_mm=RESOLUTIONS,
        template_resolutions_mm=RESOLUTIONS,
        align_shapes=False,
    )
    assert result.optimization is None
    np.testing.assert_array_equal(result.transform, np.asarray(INITIAL_TRANSFORM, dtype=np.float32))
    np.testing.assert_array_equal(bridge.membership.class_memberships, memberships_before)


def test_registration_handoff_pins_selected_java_sources() -> None:
    atlas_source = SOURCE_ROOT / "methods" / "LesionToadAtlas.java"
    wrapper_source = SOURCE_ROOT / "algorithms" / "AlgorithmLesionToads.java"
    assert hashlib.sha256(atlas_source.read_bytes()).hexdigest() == ATLAS_SOURCE_SHA
    assert hashlib.sha256(wrapper_source.read_bytes()).hexdigest() == WRAPPER_SOURCE_SHA
