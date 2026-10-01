from __future__ import annotations

import hashlib
import shutil
import struct
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_post_loop_outputs import lesion_toads_post_loop_outputs
from pymacruise.plus_lesion_inputs import prepare_plus_lesion_inputs
from test_lesion_nan_post_registration_branch import (
    LABELS,
    _branch,
    _p39ar_context,
)
from test_lesion_post_thinning_gm_membership import SPEC
from test_lesion_post_registration_output_routes import NAMES


ROOT = Path(__file__).parent
FIXTURE = ROOT / "data/synthetic/plus_lesion_transfer_crop"
HARNESS = ROOT / "java_reference/PlusLesionTransferCropRasIdentityReferenceHarness.java"
HARNESS_SHA256 = "a3fb7afa0365609e513d0cb480e169643637144e3010909bd0409befec266c94"
EDGE = 5


def _load(name: str) -> nib.Nifti1Image:
    image = nib.load(str(FIXTURE / name))
    assert isinstance(image, nib.Nifti1Image)
    return image


def _java_reference(tmp_path: Path, arrays: list[np.ndarray]) -> tuple[tuple[int, ...], list[np.ndarray]]:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA256
    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac is not None and java is not None, "JDK required for selected Java references"
    classes = tmp_path / "classes"
    classes.mkdir()
    subprocess.run([javac, "-d", str(classes), str(HARNESS)], check=True, capture_output=True, text=True)

    shape = arrays[0].shape
    input_path = tmp_path / "input.bin"
    output_path = tmp_path / "output.bin"
    with input_path.open("wb") as stream:
        stream.write(struct.pack(">4i", *shape, EDGE))
        for array in arrays:
            np.asarray(array, dtype=">f4").ravel(order="C").tofile(stream)
    subprocess.run(
        [java, "-cp", str(classes), "PlusLesionTransferCropRasIdentityReferenceHarness",
         str(input_path), str(output_path)],
        check=True, capture_output=True, text=True,
    )
    with output_path.open("rb") as stream:
        low_high = struct.unpack(">6i", stream.read(24))
        crop_shape = struct.unpack(">3i", stream.read(12))
        voxel_count = int(np.prod(crop_shape))
        cropped = [
            np.frombuffer(stream.read(voxel_count * 4), dtype=">f4").astype(np.float32)
            .reshape(crop_shape)
            for _ in arrays
        ]
        assert stream.read() == b""
    return (*low_high, *crop_shape), cropped


def test_p39ar_dura_products_transfer_and_crop_into_plus_ready_ras_inputs(tmp_path: Path) -> None:
    frames, update, registration, thinning = _p39ar_context()
    branch = _branch(
        update,
        registration,
        thinning,
        previous=float(registration.energy_ratio.previous_energy),
        iteration=3,
        n_iter_max=-3,
        max_distance=0.1,
        source_ratio=registration.energy_ratio,
    )
    assert branch.stop and branch.final_classification is not None
    membership = branch.membership_result.membership
    lesion_outputs = lesion_toads_post_loop_outputs(
        branch.final_classification,
        membership.class_memberships,
        membership.lesion_membership,
        membership.black_hole_membership,
        branch.segmentation,
        LABELS,
        NAMES,
        int(SPEC["lesion_class"]),
        result_mode="dura removal inputs",
    )
    assert lesion_outputs.source_order[-1] == "generateDuraRemovalOutputs"
    assert [product.name for product in lesion_outputs.products] == [
        "classification", "lesion_subtype", "dura_removal_inputs"
    ]

    source_grid = _load("lesion_output_grid.nii.gz")
    segmentation = _load("segmentation_reference.nii.gz")
    atlas_gm = _load("multi_atlas_gm.nii.gz")
    atlas_wm = _load("multi_atlas_wm.nii.gz")
    atlas_mask = _load("multi_atlas_mask.nii.gz")
    brain_mask = _load("brain_mask.nii.gz")
    before = {
        image.get_filename() or str(id(image)): np.asarray(image.dataobj).copy()
        for image in (source_grid, segmentation, atlas_gm, atlas_wm, atlas_mask, brain_mask)
    }

    result = prepare_plus_lesion_inputs(
        lesion_outputs,
        source_grid,
        segmentation,
        atlas_gm,
        atlas_wm,
        atlas_mask,
        brain_mask,
        crop_edge_width=EDGE,
    )

    channels = np.asarray(lesion_outputs.products[-1].values, dtype=np.float32)
    expected_java_geometry, expected_java = _java_reference(
        tmp_path,
        [
            np.asarray(atlas_gm.dataobj, dtype=np.float32),
            np.asarray(atlas_wm.dataobj, dtype=np.float32),
            np.asarray(atlas_mask.dataobj, dtype=np.float32),
            np.asarray(brain_mask.dataobj, dtype=np.float32),
            np.asarray(segmentation.dataobj, dtype=np.float32),
            channels[1],
            channels[2],
            channels[3],
        ],
    )
    assert result.crop_parameters.bounds == (
        (expected_java_geometry[0], expected_java_geometry[1]),
        (expected_java_geometry[2], expected_java_geometry[3]),
        (expected_java_geometry[4], expected_java_geometry[5]),
    )
    outputs = (
        result.multi_atlas_gm,
        result.multi_atlas_wm,
        result.multi_atlas_mask,
        result.brain_mask,
        result.segmentation,
        result.toads_gm,
        result.toads_wm_fill,
        result.toads_wm_mask,
    )
    assert tuple(outputs[0].shape) == expected_java_geometry[6:9]
    names = (
        "multi_atlas_gm", "multi_atlas_wm", "multi_atlas_mask", "brain_mask",
        "segmentation", "toads_gm", "toads_wm_fill", "toads_wm_mask",
    )
    for name, image, expected in zip(names, outputs, expected_java, strict=True):
        assert image.shape == outputs[0].shape
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_array_equal(
            np.asarray(image.dataobj, dtype=np.float32), expected, err_msg=name
        )
        assert np.allclose(image.affine, outputs[0].affine, rtol=0, atol=1e-6)
    assert not hasattr(result, "toads_csf")  # CSF is not returned by the selected plus wrapper.
    for image, original in zip(
        (source_grid, segmentation, atlas_gm, atlas_wm, atlas_mask, brain_mask),
        before.values(), strict=True,
    ):
        np.testing.assert_array_equal(image.dataobj, original)
