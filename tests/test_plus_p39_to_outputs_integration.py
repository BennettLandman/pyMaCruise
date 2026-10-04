from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_post_loop_outputs import (
    LesionToadsOutputProduct,
    LesionToadsPostLoopOutputs,
)
from pymacruise.lesion_dura_outputs import generate_lesion_dura_removal_outputs
from pymacruise.plus_lesion_inputs import prepare_plus_lesion_inputs
from pymacruise.plus_pipeline import PlusPipelineOptions, run_plus_post_tgdm_stages
from pymacruise.plus_topology_composition import run_plus_topology_and_tgdm


FIXTURE = Path(__file__).parent / "data/synthetic/plus_pre_tgdm_composition"
ATLAS_NAMES = (
    "Background", "Sulcal-CSF", "Ventricles", "Cerebellum-GM", "Cerebrum-GM",
    "Caudate", "Thalamus", "Putamen", "Brainstem", "Cerebellum-WM",
    "OuterVentricular-WM", "InterVentricular-WM", "Frontal-WM",
)


def _load(name: str) -> nib.Nifti1Image:
    image = nib.load(str(FIXTURE / name))
    assert isinstance(image, nib.Nifti1Image)
    return image


def test_p39_dura_products_cross_transfer_crop_and_plus_surface_tail() -> None:
    segmentation = _load("segmentation.nii.gz")
    atlas_gm = _load("multi_atlas_gm.nii.gz")
    atlas_wm = _load("multi_atlas_wm.nii.gz")
    atlas_mask = _load("multi_atlas_brain_mask.nii.gz")
    toads_gm = _load("toads_gm.nii.gz")
    toads_wm = _load("toads_wm.nii.gz")
    toads_mask = _load("toads_mask.nii.gz")
    toads_wm_boundary = _load("toads_wm_for_boundary.nii.gz")

    # Generate P39's four Dura channels from synthetic class memberships, then
    # send those public post-loop products through transfer/crop and plus.
    class_memberships = np.zeros((*segmentation.shape, len(ATLAS_NAMES)), dtype=np.float32)
    class_memberships[..., ATLAS_NAMES.index("Cerebrum-GM")] = np.asarray(
        toads_gm.dataobj, dtype=np.float32
    )
    class_memberships[..., ATLAS_NAMES.index("Brainstem")] = np.asarray(
        toads_wm.dataobj, dtype=np.float32
    )
    class_indices = np.full(segmentation.shape, ATLAS_NAMES.index("Brainstem"), dtype=np.int16)
    generated = generate_lesion_dura_removal_outputs(
        class_memberships,
        class_indices,
        ATLAS_NAMES,
        np.zeros(segmentation.shape, dtype=np.float32),
        np.zeros(segmentation.shape, dtype=np.float32),
    )
    channels = generated.as_channels()
    np.testing.assert_array_equal(channels[1], toads_gm.dataobj)
    np.testing.assert_array_equal(channels[2], toads_wm.dataobj)
    assert np.all(channels[3] == 1.0)
    labels = np.asarray(segmentation.dataobj, dtype=np.int8)
    lesion_outputs = LesionToadsPostLoopOutputs(
        classification=labels,
        membership_classification=None,
        lesion_classification=np.zeros(labels.shape, dtype=np.int16),
        result_mode="dura removal inputs",
        products=(LesionToadsOutputProduct("dura_removal_inputs", channels),),
        source_order=("generateDuraRemovalOutputs",),
    )
    source_images = (
        segmentation, atlas_gm, atlas_wm, atlas_mask, toads_gm, toads_wm,
        toads_mask, toads_wm_boundary,
    )
    originals = tuple(np.asarray(image.dataobj).copy() for image in source_images)

    prepared = prepare_plus_lesion_inputs(
        lesion_outputs,
        segmentation,
        segmentation,
        atlas_gm,
        atlas_wm,
        atlas_mask,
        atlas_mask,
        crop_edge_width=5,
    )
    cropped = run_plus_topology_and_tgdm(
        prepared.segmentation,
        prepared.toads_wm_fill,
        prepared.multi_atlas_gm,
        prepared.multi_atlas_wm,
        prepared.multi_atlas_mask,
        prepared.toads_gm,
        prepared.toads_wm_fill,
        prepared.toads_wm_mask,
        options=PlusPipelineOptions(),
    )
    post = run_plus_post_tgdm_stages(cropped.tgdm, prepared.crop_parameters)

    assert prepared.crop_parameters.original_shape == segmentation.shape
    assert cropped.tgdm_initial_membership is cropped.topology.atlas_refined
    assert cropped.tgdm_initial_membership.shape == prepared.segmentation.shape
    for image in (
        prepared.multi_atlas_gm, prepared.multi_atlas_wm, prepared.multi_atlas_mask,
        prepared.brain_mask, prepared.segmentation, prepared.toads_gm,
        prepared.toads_wm_fill, prepared.toads_wm_mask,
        cropped.tissue.gm, cropped.tissue.wm, cropped.topology.atlas_refined,
        cropped.ace.enhanced_gray_matter, cropped.gvf.vector_field,
        cropped.tgdm.inner, cropped.tgdm.central, cropped.tgdm.outer,
        post.level_sets.inner, post.level_sets.central, post.level_sets.outer,
        post.thickness.volume.image,
    ):
        assert image is not None
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        assert np.isfinite(np.asarray(image.dataobj)).all()
    assert post.level_sets.inner.shape == segmentation.shape
    assert post.thickness.volume.image.shape == segmentation.shape
    assert post.thickness.surface.mean_mm >= 0.0
    for mesh in (post.surfaces.inner, post.surfaces.central, post.surfaces.outer):
        assert mesh.coordinate_frame == "RAS"
        assert mesh.units == "mm"
        assert np.isfinite(mesh.vertices).all()

    for image, before in zip(source_images, originals, strict=True):
        np.testing.assert_array_equal(image.dataobj, before)
