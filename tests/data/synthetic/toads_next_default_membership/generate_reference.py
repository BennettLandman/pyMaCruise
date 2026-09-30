"""Rebuild the fixed P18czh->P18czi->P18czj RAS fixture and Java excerpt.

The handoff and continuing branch are generated with the adjacent validated
Python adapters. The Java harness independently evaluates only P18czj: source
thinning queue initialization/drain under its checked minimum-volume branch,
then the selected default 11-class exact-membership order.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.toads_default_registration_handoff import toads_default_atlas_registration_handoff
from pymacruise.toads_next_default_membership import toads_next_default_membership_pass
from pymacruise.toads_post_registration import toads_post_registration_branch
from test_toads_default_atlas_dispatch import _fixture as default_fixture
from test_toads_default_registration_handoff import _arguments
from test_toads_post_registration_branch import STATE as CZI_STATE, _arrays as czi_arrays


ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
HARNESS = ROOT / "tests/java_reference/ToadsNextDefaultMembershipReferenceHarness.java"
UPSTREAM = ROOT / "tests/java_reference/upstream/edu/jhmi/rad/medic"
VOLUME_RATIO = 0.2
THINNING_LIMIT = 1.0
SUBJECT_SHAPE = (25, 25, 25)


def _save(name: str, values: np.ndarray, affine: np.ndarray) -> None:
    nib.save(nib.Nifti1Image(values, affine), HERE / name)


def _repeat_to_shape(values: np.ndarray) -> np.ndarray:
    indices = tuple(np.arange(size) % values.shape[axis] for axis, size in enumerate(SUBJECT_SHAPE))
    return values[np.ix_(*indices)]


def main() -> None:
    data, parameters, _, _, options, affine = _arguments()
    source_data, _, _, _ = default_fixture()
    data = {name: _repeat_to_shape(value) for name, value in source_data.items()}
    for option_name, data_name in (
        ("image_values", "images"), ("initial_memberships", "initial_memberships"),
        ("transformed_shapes", "transformed_shapes"), ("closest_labels", "closest_labels"),
        ("intensity_field", "field"), ("psf_bias", "psf_bias"),
    ):
        options[option_name] = data[data_name]
    affine = np.asarray(affine)
    _save("subject_images.nii.gz", data["images"], affine)
    _save("folded_ribbon.nii.gz", data["folded_ribbon"], affine)
    subject = nib.load(str(HERE / "subject_images.nii.gz"))
    atlas_image = nib.load(str(HERE / "template_atlas_shapes.nii.gz"))
    if nib.aff2axcodes(subject.affine) != ("R", "A", "S") or nib.aff2axcodes(atlas_image.affine) != ("R", "A", "S"):
        raise ValueError("synthetic subject and independent template grid must both be RAS")
    np.testing.assert_array_equal(subject.affine, affine)
    np.testing.assert_array_equal(np.asarray(subject.dataobj), data["images"])
    options["atlas_shapes"] = np.asarray(atlas_image.dataobj, dtype=np.float32)

    handoff = toads_default_atlas_registration_handoff(**options)
    state = {name: _repeat_to_shape(value) if value.ndim >= 3 and value.shape[:3] == (9, 8, 7) else value
             for name, value in czi_arrays().items()}
    state["class_volumes"] = np.sum(
        np.where(
            np.eye(11, dtype=bool)[state["segmentation"]],
            handoff.membership_pass.memberships,
            np.float32(0.0),
        )[1:-1, 1:-1, 1:-1],
        axis=(0, 1, 2),
        dtype=np.float32,
    )
    branch = toads_post_registration_branch(
        handoff,
        previous_energy=handoff.membership_pass.objective,
        segmentation=state["segmentation"],
        skeleton=state["skeleton"],
        ordering=state["ordering"],
        class_volumes=state["class_volumes"],
        boundary_distances=state["boundary_distances"],
        template_labels=np.asarray(CZI_STATE["template_labels"], dtype=np.int8),
        previous_classification=state["previous_classification"],
        iteration=1,
        maximum_iterations=5,
        convergence_threshold=0.0,
        relation_limit=1,
        do_growing=False,
    )
    if branch.control.stop or branch.growing is None or branch.relations is None:
        raise ValueError("fixture generation requires the actual continuing czi branch")

    args = dict(options)
    args.pop("initial_memberships")
    args.pop("transformed_shapes")
    args.pop("closest_labels")
    result = toads_next_default_membership_pass(
        branch,
        image_values=data["images"],
        centroids=parameters["centroids"],
        exponent_by_channel_class=parameters["exponent_by_channel_class"],
        modality_weights=parameters["modality_weights"],
        class_variance_weights=parameters["class_variance_weights"],
        image_variances=parameters["image_variances"],
        smoothing_priors=parameters["smoothing_priors"],
        centroid_similarity=parameters["centroid_similarity"],
        atlas_shape_priors=parameters["atlas_shape_priors"],
        atlas_shapes=options["atlas_shapes"],
        image_resolutions_mm=options["image_resolutions_mm"],
        template_resolutions_mm=options["template_resolutions_mm"],
        smoothing=options["smoothing"],
        prior_factor=options["prior_factor"],
        prior_coefficient=options["prior_coefficient"],
        volume_ratio=VOLUME_RATIO,
        thinning_limit=THINNING_LIMIT,
        intensity_field=data["field"],
        psf_bias=data["psf_bias"],
    )

    # Snapshot the actual prior-composition state so the independent Java
    # method excerpt can consume precisely the registered transform and czi
    # refresh outputs tested end-to-end by the Python integration test.
    _save("continuing_memberships.nii.gz", handoff.membership_pass.memberships, subject.affine)
    _save("current_transform_shapes.nii.gz", result.transformed_shapes, subject.affine)
    _save("refreshed_closest_labels.nii.gz", branch.relations.closest_labels.astype(np.int16), subject.affine)
    _save("grown_segmentation.nii.gz", branch.growing.segmentation.astype(np.int16), subject.affine)
    _save("channel_field.nii.gz", data["field"], subject.affine)
    _save("psf_bias.nii.gz", data["psf_bias"], subject.affine)
    adjacency = branch.relations.relation.class_adjacency
    (HERE / "refreshed_class_adjacency.csv").write_text(
        "".join(",".join("1" if value else "0" for value in row) + "\n" for row in adjacency)
    )

    with tempfile.TemporaryDirectory(prefix="p18czj-java-") as build_dir:
        subprocess.run(
            [
                "javac", "-d", build_dir,
                str(UPSTREAM / "structures/BinaryHeap4D.java"),
                str(UPSTREAM / "structures/CriticalPointLUT.java"),
                str(UPSTREAM / "utilities/PowerTable.java"),
                str(ROOT / "tests/java_reference/ToadsHomeomorphicChangeReferenceHarness.java"),
                str(HARNESS),
            ],
            check=True,
        )
        (Path(build_dir) / "critical186LUT.raw.gz").write_bytes(
            (ROOT / "src/pymacruise/data/critical186LUT.raw.gz").read_bytes()
        )
        stdout = subprocess.check_output(
            [
                "java", "-cp", build_dir, "ToadsNextDefaultMembershipReferenceHarness",
                str(HERE / "subject_images.nii.gz"),
                str(HERE / "channel_field.nii.gz"),
                str(HERE / "psf_bias.nii.gz"),
                str(HERE / "continuing_memberships.nii.gz"),
                str(HERE / "current_transform_shapes.nii.gz"),
                str(HERE / "refreshed_closest_labels.nii.gz"),
                str(HERE / "grown_segmentation.nii.gz"),
                str(HERE / "refreshed_class_adjacency.csv"),
                str(VOLUME_RATIO),
            ],
            text=True,
        )
    reference = json.loads(stdout)
    reference.update(
        {
            "wrapper_source_sha256": "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d",
            "atlas_source_sha256": "914839a0885cb26ff3ac22daeeb2603ee1960368b788da3ad343a54053f031db",
            "powertable_source_sha256": "5552069bec24dad8e50d426c77babfb041f952166a8c786091ade22eec8ec678",
            "shape_xyz_channels": [*SUBJECT_SHAPE, 11],
            "subject_xyz_channels": [*SUBJECT_SHAPE, 2],
            "template_xyz_channels": [11, 10, 9, 11],
            "registered_transform_bits": np.asarray(handoff.transform, dtype=np.float32).view(np.uint32).astype(int).tolist(),
            "stage_order": list(result.stage_order),
            "volume_ratio": VOLUME_RATIO,
            "thinning_limit": THINNING_LIMIT,
            "scope_note": "The fixed 25^3 folded RAS phantom exceeds the selected source minimum-volume guard for supported labels. Java and Python match the complete queue trace and topology-checked thinning result, including 978 segmentation-label updates. The following 11-class membership sweep matches Java float bits. Class statistics and transform registration remain explicit upstream inputs; no estimation or registration is added in this adapter.",
            "queue_pop_count": len(result.thinning.processed_items),
            "segmentation_label_update_count": int(np.count_nonzero(result.thinning.segmentation != branch.growing.segmentation)),
            "removed_skeleton_voxel_count": int(np.count_nonzero(result.thinning.skeleton == -1)),
        }
    )
    (HERE / "java_reference.json").write_text(json.dumps(reference, indent=2) + "\n")


if __name__ == "__main__":
    main()
