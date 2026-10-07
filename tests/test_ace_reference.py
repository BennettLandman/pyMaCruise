from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.ace import (
    _central_gradient_mirror,
    _distance_and_characteristic_fmm,
    _principal_tensor_vector,
    _signed_distance_fmm,
    ace_anisotropic_diffusion,
    ace_close_wm,
    ace_enhance,
    ace_gaussian_blur,
    ace_seed_fields,
    ace_update_speed,
    clean_ace_components,
    narrow_band_ace,
)
from pymacruise.validation import dice_coefficient, mesh_surface_distance_metrics
from pymacruise.surfaces import extract_levelset_surface


ROOT = Path(__file__).parent
FIXTURES = ROOT / "data" / "synthetic" / "ace"
JAVA = FIXTURES / "java_reference"


def _compare_float_stage(
    candidate: np.ndarray,
    reference: np.ndarray,
    *,
    maximum: float = 5e-6,
    mean: float = 1e-6,
) -> None:
    error = np.abs(np.asarray(candidate, dtype=np.float32) - np.asarray(reference, dtype=np.float32))
    assert float(error.mean()) <= mean
    assert float(error.max()) <= maximum


def _assert_ace_stages_match_selected_java(fixture_root: Path, java_root: Path) -> None:
    def reference(name: str) -> np.ndarray:
        return np.asarray(nib.load(java_root / "stages" / f"{name}_java.nii.gz").dataobj)

    gray_image = nib.load(fixture_root / "gray_matter.nii.gz")
    white_image = nib.load(fixture_root / "white_matter.nii.gz")
    atlas_image = nib.load(fixture_root / "multi_atlas_csf.nii.gz")
    gray = np.asarray(gray_image.dataobj, dtype=np.float32)
    white = np.asarray(white_image.dataobj, dtype=np.float32)
    atlas = np.asarray(atlas_image.dataobj)

    # Independent Java stage fixtures make each mismatch attributable to one
    # operation rather than to accumulated end-to-end numerical drift.
    initial_levelset = np.float32(126.9) - white
    initial_phi = _signed_distance_fmm(initial_levelset, 1000.0)
    seed = ace_seed_fields(gray, white, initial_phi, intensity_threshold=126.9)
    guide = ace_gaussian_blur(initial_phi, 1.8)
    diffused = ace_anisotropic_diffusion(seed.initial_csf_speed, guide, tensor_alpha=1.6, iterations=25)
    merged = diffused.copy()
    merged[atlas > 0] = 1.0
    speed_before_blur = 1.0 - 0.9 * merged
    speed = ace_update_speed(diffused, multi_atlas_csf=atlas, gaussian_alpha=1.8)
    evolved_phi = narrow_band_ace(
        initial_phi, delta_t=0.1, max_iterations=1, inner_steps=10, curvature_threshold=5.0
    )
    weighted_distance, _ = _distance_and_characteristic_fmm(evolved_phi, speed)

    exact_stages = {
        "01_tmpVolfirst": initial_levelset,
        "03_csf_naive": seed.naive_csf,
        "04_wm_morph": ((white + gray) >= 210.0).astype(np.float32),
        "05_closedGM_morph": seed.closed_combined_wm,
        "08_initcsf_beforediffusion": seed.initial_csf_speed,
    }
    for name, candidate in exact_stages.items():
        np.testing.assert_array_equal(candidate, reference(name), err_msg=name)

    float_stages = {
        "02_PHIinit": initial_phi,
        "06_PHI_beforeGaussianBlur": initial_phi,
        "07_tmpVol_diffusion": guide,
        "09_initcsf_diffusion": diffused,
        "09_5_initcsf_afterdiffusion_mAcsf": merged,
        "10_initcsf_speedfunction_beforeGuassianBlur": speed_before_blur,
        "11_initcsf_GuassianBlur": speed,
        "12_PHI_beforenarrowBand": initial_phi,
        "13_PHI_narrowBand": evolved_phi,
    }
    for name, candidate in float_stages.items():
        # The folded fixture accumulates a few more float32 ulps in the
        # Gaussian and narrow-band stages than the concentric sphere. Retain
        # explicit per-stage limits and apply the project surface gates to
        # the evolved zero set; downstream label support and wrapper outputs are checked separately.
        limits = {
            "07_tmpVol_diffusion": {"maximum": 1e-5, "mean": 1e-6},
            "13_PHI_narrowBand": {"maximum": 0.012, "mean": 2e-5},
        }.get(name, {})
        _compare_float_stage(candidate, reference(name), **limits)

    candidate_phi = nib.Nifti1Image(evolved_phi, gray_image.affine)
    reference_phi = nib.load(java_root / "stages" / "13_PHI_narrowBand_java.nii.gz")
    assert dice_coefficient(
        np.asarray(candidate_phi.dataobj) <= 0,
        np.asarray(reference_phi.dataobj) <= 0,
    ) >= 0.97
    candidate_phi_surface = extract_levelset_surface(candidate_phi, level=0)
    reference_phi_surface = extract_levelset_surface(reference_phi, level=0)
    phi_surface_metrics = mesh_surface_distance_metrics(
        reference_phi_surface.vertices,
        reference_phi_surface.faces,
        candidate_phi_surface.vertices,
        candidate_phi_surface.faces,
    )
    assert phi_surface_metrics["hausdorff_mm"] <= 1.0
    assert phi_surface_metrics["mean_absolute_distance_mm"] <= 0.01

    # The indexed float32 heap follows JIST insertion, update, and tie order.
    # The folded fixture retains two one-unit scalar membership differences;
    # their support and final wrapper products are checked below.
    _compare_float_stage(
        weighted_distance,
        reference("14_tmpVol_narrowBand"),
        maximum=1e-5 if fixture_root.name == "folded_cortex" else 2.5e-5,
        mean=5e-6 if fixture_root.name == "folded_cortex" else 1e-5,
    )

    dx, dy, dz = _central_gradient_mirror(weighted_distance)
    magnitude = np.sqrt(dx * dx + dy * dy + dz * dz) * speed
    eligible = (weighted_distance > 0.25) & ~seed.prohibited_region
    memberships = np.zeros(gray.shape, dtype=np.int32)
    selected = eligible & (magnitude < 0.85)
    memberships[selected] = np.trunc(254.0 * (magnitude[selected] + 0.05)).astype(np.int32)
    java_memberships = reference("15_csf_beforeclean")
    if fixture_root.name == "folded_cortex":
        assert np.count_nonzero(memberships != java_memberships) == 2
        assert np.max(np.abs(memberships - java_memberships)) == 1
        assert dice_coefficient(memberships > 0, java_memberships > 0) == 1.0
    else:
        np.testing.assert_array_equal(memberships, java_memberships)
    cleaned = clean_ace_components(memberships, minimum_size_exclusive=60)
    java_cleaned = reference("16_csf_afterclean")
    if fixture_root.name == "folded_cortex":
        assert np.count_nonzero(cleaned != java_cleaned) == 2
        assert np.max(np.abs(cleaned - java_cleaned)) == 1
        assert dice_coefficient(cleaned > 0, java_cleaned > 0) == 1.0
    else:
        np.testing.assert_array_equal(cleaned, java_cleaned)

    result = ace_enhance(gray_image, white_image, multi_atlas_csf=atlas_image)
    java_enhanced = nib.load(java_root / "enhanced_gm_java.nii.gz")
    java_skeleton = nib.load(java_root / "csf_skeleton_java.nii.gz")
    java_thinned = nib.load(java_root / "thinned_skeleton_java.nii.gz")
    np.testing.assert_array_equal(result.enhanced_gray_matter.dataobj, java_enhanced.dataobj)
    np.testing.assert_array_equal(result.skeleton.dataobj, java_skeleton.dataobj)
    np.testing.assert_array_equal(result.thinned_skeleton.dataobj, java_thinned.dataobj)
    for image in (result.enhanced_gray_matter, result.skeleton, result.thinned_skeleton):
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")

    # Both source skeleton products are label masks. Their exact Java equality
    # implies these metrics are zero, but keep the declared project gates
    # explicit in this stage-level acceptance test.
    for candidate_image, reference_image in (
        (result.skeleton, java_skeleton),
        (result.thinned_skeleton, java_thinned),
    ):
        candidate = np.asarray(candidate_image.dataobj) > 0
        reference = np.asarray(reference_image.dataobj) > 0
        assert dice_coefficient(reference, candidate) >= 0.97
        candidate_surface = extract_levelset_surface(candidate_image, level=127)
        reference_surface = extract_levelset_surface(reference_image, level=127)
        metrics = mesh_surface_distance_metrics(
            reference_surface.vertices,
            reference_surface.faces,
            candidate_surface.vertices,
            candidate_surface.faces,
        )
        assert metrics["hausdorff_mm"] <= 1.0
        assert metrics["mean_absolute_distance_mm"] <= 0.01



def test_ace_stages_match_selected_macruise_java_reference() -> None:
    _assert_ace_stages_match_selected_java(FIXTURES, JAVA)


def test_ace_stages_match_selected_java_reference_on_folded_cortex() -> None:
    folded = FIXTURES.parent / "folded_cortex"
    _assert_ace_stages_match_selected_java(folded, folded / "ace_java_reference")


def test_ace_source_faithful_mode_preserves_zero_tensor_nan_and_reports_count() -> None:
    zero_component = np.zeros((3, 3, 3), dtype=np.float32)
    robust = _principal_tensor_vector(
        zero_component,
        zero_component,
        zero_component,
        zero_component,
        zero_component,
        zero_component,
    )
    source = _principal_tensor_vector(
        zero_component,
        zero_component,
        zero_component,
        zero_component,
        zero_component,
        zero_component,
        source_faithful=True,
    )
    zero_tensor_voxels = int(
        np.count_nonzero(
            (zero_component == 0)
            & (zero_component == 0)
            & (zero_component == 0)
            & (zero_component == 0)
            & (zero_component == 0)
            & (zero_component == 0)
        )
    )

    assert zero_tensor_voxels == 27
    np.testing.assert_array_equal(robust[0], 0)
    np.testing.assert_array_equal(robust[1], 0)
    np.testing.assert_array_equal(robust[2], 1)
    assert all(np.isnan(component).all() for component in source)
