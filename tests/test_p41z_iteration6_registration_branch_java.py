from __future__ import annotations

import json
import hashlib
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from compare_p41z_iteration6_membership_java import _stats_from_fixture
from pymacruise.lesion_next_iteration_membership import run_lesion_toads_next_gm_membership_update
from pymacruise.lesion_post_registration import lesion_toads_post_registration_branch
from pymacruise.lesion_registration_handoff import register_lesion_membership_handoff
from pymacruise.lesion_registration_handoff import _bounds_for_transform
from pymacruise.lesion_relations import lesion_compute_relations
from pymacruise.lesion_standard13_volume import STANDARD13_TYPES


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data/synthetic/p41z_iteration6_membership_java"
ATLAS = ROOT / "data/synthetic/plus_vandy_atlas_aligned_fine"
HARNESS = ROOT / "java_reference/LesionNextIterationBranchReferenceHarness.java"
REGISTRATION_HARNESS = ROOT / "java_reference/P41zIteration6RegistrationReferenceHarness.java"
OUTPUT_HARNESS = ROOT / "java_reference/LesionPostRegistrationOutputsReferenceHarness.java"
BOUNDS_HARNESS = ROOT / "java_reference/LesionShapeBoundsReferenceHarness.java"
NUMERICS = ROOT / "java_reference/upstream/edu/jhmi/rad/medic/utilities/Numerics.java"
REFERENCE_ROOT = ROOT / "java_reference"
UPSTREAM = REFERENCE_ROOT / "upstream/edu/jhmi/rad/medic"
OUTPUT_HARNESS_SHA256 = "c65f4efe46d96518c50acbd437482fa4118b88a87475561ce68f83569a2b0e76"
REGISTRATION_HARNESS_SHA256 = "da1022cdc5fd9fb1fce2c4456e9427187ebe5c546944e28985df92f112f4b366"
PREVIOUS_SUMMARY_SHA = "f03d9f0d0831022335a3296bd084e229f172ff49d86645827ec05199eb46aac0"
ATLAS_NAMES = (
    "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM", "Cerebrum-GM",
    "Caudate", "Thalamus", "Putamen", "Brainstem", "Cerebellum-WM",
    "Cerebrum-WM", "InterVentricular-WM", "Inferior-WM",
)


@pytest.mark.skipif(
    os.environ.get("PYMACRUISE_RUN_P41Z_ITER6_REGISTRATION_BRANCH") != "1",
    reason="set PYMACRUISE_RUN_P41Z_ITER6_REGISTRATION_BRANCH=1 for fixed-state registration/branch replay",
)
def test_p41z_iteration6_registration_and_wrapper_branch_match_references() -> None:
    if shutil.which("java") is None or shutil.which("javac") is None:
        pytest.skip("Java runtime and compiler are required for selected-source comparison")
    case = json.loads((FIXTURE / "case.json").read_text())
    assert hashlib.sha256(REGISTRATION_HARNESS.read_bytes()).hexdigest() == REGISTRATION_HARNESS_SHA256
    assert case["provenance"]["iteration5_summary_sha256"] == PREVIOUS_SUMMARY_SHA
    previous_objective = np.float32(case["wrapper_ratio_previous_objective"])
    previous_objective_bits = int(previous_objective.view(np.uint32))
    assert previous_objective_bits == 1164507634

    arrays, statistics, _ = _stats_from_fixture()
    update = run_lesion_toads_next_gm_membership_update(
        statistics,
        arrays["images"],
        arrays["shape_values"],
        np.ones(13, dtype=np.float32),
        np.ones(13, dtype=np.float32),
        arrays["direct_relations"],
        arrays["class_memberships"],
        arrays["lesion"],
        arrays["black"],
        arrays["outlier"],
        arrays["lesion_relation"],
        arrays["distance_maps"][..., 0],
        arrays["distance_maps"][..., 1],
        arrays["distance_maps"][..., 2],
        previous_energy=float(previous_objective),
        outlier_energy=0.1,
        lesion_power_exponent=1.2,
        spread=0.5,
        source_zero_variance_nan=True,
    )
    expected_state = case["registration_reference_iteration6"]
    assert int(np.asarray(update.membership.objective, dtype=np.float32).view(np.uint32)) == int(
        np.asarray(expected_state["transition_objective"], dtype=np.float32).view(np.uint32)
    )

    shape_priors = tuple(
        np.asarray(nib.load(str(ATLAS / f"atlas_{name}.nii.gz")).dataobj, dtype=np.float32)
        for name in ATLAS_NAMES
    )
    image_zooms = tuple(
        float(v)
        for v in nib.load(str(FIXTURE / "images_t1_flair.nii.gz")).header.get_zooms()[:3]
    )
    template_zooms = tuple(
        float(v)
        for v in nib.load(str(ATLAS / "atlas_Background.nii.gz")).header.get_zooms()[:3]
    )
    registration = register_lesion_membership_handoff(
        update,
        shape_priors,
        case["registration_initial_transform_iteration5"],
        image_resolutions_mm=image_zooms,
        template_resolutions_mm=template_zooms,
        align_shapes=True,
        init_align_iter=10,
    )

    # Run the selected Java registerAllShapes() excerpt against the exact
    # updated memberships and bundled atlas arrays. Binary layout is big
    # endian and voxel-major for memberships, matching the Java nested arrays.
    registration_input = bytearray()
    dims = (*update.membership.class_memberships.shape[:3], *shape_priors[0].shape, 13)
    registration_input.extend(struct.pack(">7i", *dims))
    registration_input.extend(struct.pack(">6f", *image_zooms, *template_zooms))
    registration_input.extend(struct.pack(">6f", *case["registration_initial_transform_iteration5"]))
    for prior in shape_priors:
        registration_input.extend(np.asarray(prior, dtype=">f4").tobytes(order="C"))
    registration_input.extend(
        np.asarray(update.membership.class_memberships, dtype=">f4").tobytes(order="C")
    )
    upstream_sources = (
        REFERENCE_ROOT / "stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java",
        UPSTREAM / "utilities/Numerics.java",
        UPSTREAM / "structures/RotationMatrix.java",
        UPSTREAM / "methods/ParametricTransform.java",
        UPSTREAM / "libraries/ImageFunctions.java",
        REGISTRATION_HARNESS,
    )
    with tempfile.TemporaryDirectory(prefix="p41z-iter6-registration-java-") as temporary:
        work = Path(temporary)
        classes = work / "classes"
        classes.mkdir()
        input_path = work / "registration.bin"
        input_path.write_bytes(registration_input)
        subprocess.run(
            ["javac", "-d", str(classes), *(str(path) for path in upstream_sources)],
            check=True, capture_output=True, text=True,
        )
        reference = json.loads(subprocess.run(
            ["java", "-Xmx2g", "-cp", str(classes),
             "P41zIteration6RegistrationReferenceHarness", str(input_path)],
            check=True, capture_output=True, text=True,
        ).stdout)

    optimization = registration.optimization
    assert optimization is not None
    assert reference["stop_reason"] == optimization.stop_reason
    np.testing.assert_array_equal(registration.transform, np.asarray(reference["transform"], dtype=np.float32))
    assert len(reference["iterations"]) == len(optimization.iterations)
    for java_step, python_step in zip(reference["iterations"], optimization.iterations, strict=True):
        assert java_step["accepted"] is python_step.result.accepted
        assert java_step["rejected_trials"] == python_step.result.rejected_trials
        np.testing.assert_array_equal(python_step.trial_transforms, np.asarray(
            [trial["transform"] for trial in java_step["trials"]], dtype=np.float32
        ))
        np.testing.assert_array_equal(python_step.trial_energies, np.asarray(
            [trial["energy"] for trial in java_step["trials"]], dtype=np.float32
        ))
        np.testing.assert_array_equal(python_step.result.transform, np.asarray(java_step["transform"], dtype=np.float32))
        assert int(np.float32(python_step.result.energy).view(np.uint32)) == int(
            np.float32(java_step["energy"]).view(np.uint32)
        )
        assert int(np.float32(python_step.result.relative_change).view(np.uint32)) == int(
            np.float32(java_step["relative_change"]).view(np.uint32)
        )
        assert int(np.float32(python_step.result.lambda_value).view(np.uint32)) == int(
            np.float32(java_step["lambda"]).view(np.uint32)
        )
    assert len(reference["bounds"]) == len(registration.shape_bounds)
    for java_bounds, python_bounds in zip(reference["bounds"], registration.shape_bounds, strict=True):
        assert python_bounds.lower_xyz == tuple(int(value) for value in java_bounds["lower_xyz"])
        assert python_bounds.upper_xyz == tuple(int(value) for value in java_bounds["upper_xyz"])
        assert python_bounds.support_count == java_bounds["support_count"]

    np.testing.assert_array_equal(
        registration.transform, np.asarray(expected_state["transform"], dtype=np.float32)
    )

    segmentation_image = nib.load(str(FIXTURE / "thinned_segmentation.nii.gz"))
    segmentation = np.asarray(segmentation_image.dataobj, dtype=np.int16)
    empty_skeleton = np.zeros(segmentation.shape, dtype=np.int16)
    empty_ordering = np.zeros(segmentation.shape, dtype=np.float32)
    branch = lesion_toads_post_registration_branch(
        update,
        registration,
        previous_energy=float(previous_objective),
        segmentation=segmentation,
        skeleton=empty_skeleton,
        ordering=empty_ordering,
        class_volume=np.ones(13, dtype=np.float32),
        object_types=STANDARD13_TYPES,
        classification=np.zeros(segmentation.shape, dtype=np.int8),
        template_labels=np.arange(13, dtype=np.int8),
        lesion_class=10,
        iteration=6,
        n_iter_max=99,
        max_distance=0.0001,
        algorithm_mode="standard",
        relation_limit=1,
        source_energy_ratio=update.energy_ratio,
        retain_growth_event_states=False,
    )

    branch_input = struct.pack(">fffii?", previous_objective, update.membership.objective,
                               np.float32(0.0001), 6, 99, False)
    with tempfile.TemporaryDirectory(prefix="p41z-iter6-branch-java-") as temporary:
        work = Path(temporary)
        classes = work / "classes"
        classes.mkdir()
        input_path = work / "control.bin"
        input_path.write_bytes(branch_input)
        subprocess.run(
            ["javac", "-d", str(classes), str(NUMERICS), str(HARNESS)],
            check=True, capture_output=True, text=True,
        )
        reference = json.loads(subprocess.run(
            ["java", "-cp", str(classes), "LesionNextIterationBranchReferenceHarness", str(input_path)],
            check=True, capture_output=True, text=True,
        ).stdout)

    assert reference["source_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert reference["ratio_bits"] == f"{int(branch.energy_ratio.view(np.uint32)):08x}"
    assert branch.stop is reference["stop"] is False
    assert branch.run_growing is reference["run_growing"] is True
    assert branch.growing is not None and branch.growing.events == ()
    assert (branch.relations is not None) is reference["refresh_relations"] is True


def test_p41z_saved_thinned_state_relation_refresh_matches_selected_java() -> None:
    """Compare computeRelations on the fixture segmentation, independent of unavailable growth state."""
    if shutil.which("java") is None or shutil.which("javac") is None:
        pytest.skip("Java runtime and compiler are required for selected-source comparison")
    assert hashlib.sha256(OUTPUT_HARNESS.read_bytes()).hexdigest() == OUTPUT_HARNESS_SHA256
    segmentation_image = nib.load(str(FIXTURE / "thinned_segmentation.nii.gz"))
    segmentation = np.asarray(segmentation_image.dataobj, dtype=np.int8)
    stored_memberships = np.asarray(
        nib.load(str(FIXTURE / "class_memberships_13ch.nii.gz")).dataobj,
        dtype=np.float32,
    )
    assert stored_memberships.shape == (*segmentation.shape, 13)
    classes, limit = 13, 3  # 14 source diffusion passes, as in the saved checkpoint.
    blank = np.zeros(segmentation.shape, dtype=np.float32)
    template_labels = np.arange(classes, dtype=np.int8)
    input_bytes = b"".join((
        struct.pack(">6i", *segmentation.shape, classes, 10, limit),
        template_labels.tobytes(),
        np.asarray(segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.zeros(segmentation.shape, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(stored_memberships, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(blank, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(blank, dtype=">f4").ravel(order="C").tobytes(),
    ))
    with tempfile.TemporaryDirectory(prefix="p41z-iter6-relations-java-") as temporary:
        work = Path(temporary)
        classes_dir = work / "classes"
        classes_dir.mkdir()
        input_path = work / "relation.bin"
        input_path.write_bytes(input_bytes)
        subprocess.run(
            ["javac", "-d", str(classes_dir), str(OUTPUT_HARNESS)],
            check=True, capture_output=True, text=True,
        )
        expected = json.loads(subprocess.run(
            ["java", "-Xmx1g", "-cp", str(classes_dir),
             "LesionPostRegistrationOutputsReferenceHarness", str(input_path)],
            check=True, capture_output=True, text=True,
        ).stdout)
    actual = lesion_compute_relations(segmentation, classes, limit)
    np.testing.assert_array_equal(
        actual.relation_map.ravel(), np.asarray(expected["relation_map"], dtype=np.int16)
    )
    np.testing.assert_array_equal(
        actual.relation_membership.ravel(),
        np.asarray(expected["relation_membership"], dtype=np.bool_),
    )
    assert actual.relation_count == expected["relation_count"]
    assert actual.diffusion_iterations == 14


def test_shape_bounds_transform_rounding_matches_java_on_compact_phantom(tmp_path: Path) -> None:
    """Keep the positive-support boundary mapper source-matched on a small case."""
    if shutil.which("java") is None or shutil.which("javac") is None:
        pytest.skip("Java runtime and compiler are required for selected-source comparison")
    image_shape = (11, 9, 7)
    template_shape = (15, 13, 11)
    image_res = (1.1, 0.9, 1.3)
    template_res = (0.8, 1.2, 0.75)
    transform = np.asarray((0.08, -0.06, 0.025, 2.2, -1.1, 0.8), dtype=np.float32)
    x, y, z = np.indices(template_shape, dtype=np.float32)
    prior = np.asarray(0.11 * x + 0.037 * y + 0.019 * z - 1.41, dtype=np.float32)
    binary = b"".join((
        struct.pack(">6i", *image_shape, *template_shape),
        struct.pack(">6f", *image_res, *template_res),
        struct.pack(">6f", *transform),
        np.asarray(prior, dtype=">f4").tobytes(order="C"),
    ))
    sources = (
        REFERENCE_ROOT / "stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java",
        UPSTREAM / "utilities/Numerics.java",
        UPSTREAM / "structures/RotationMatrix.java",
        UPSTREAM / "methods/ParametricTransform.java",
        UPSTREAM / "libraries/ImageFunctions.java",
        BOUNDS_HARNESS,
    )
    classes_dir = tmp_path / "classes"
    classes_dir.mkdir()
    input_path = tmp_path / "bounds.bin"
    input_path.write_bytes(binary)
    subprocess.run(
        ["javac", "-d", str(classes_dir), *(str(path) for path in sources)],
        check=True, capture_output=True, text=True,
    )
    expected = json.loads(subprocess.run(
        ["java", "-cp", str(classes_dir), "LesionShapeBoundsReferenceHarness", str(input_path)],
        check=True, capture_output=True, text=True,
    ).stdout)
    actual = _bounds_for_transform(prior, image_shape, image_res, template_res, transform)
    assert actual.lower_xyz == tuple(expected["lower_xyz"])
    assert actual.upper_xyz == tuple(expected["upper_xyz"])
    assert actual.support_count == expected["support_count"]
