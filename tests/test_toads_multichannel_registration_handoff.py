from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess

import nibabel as nib
import numpy as np
import pytest

from pymacruise.toads_default_atlas_membership import toads_default_atlas_membership_pass
from pymacruise.toads_default_exact_centroids import toads_default_exact_centroid_update
from pymacruise.toads_default_exact_variances import toads_default_exact_class_variances
from pymacruise.toads_default_membership_registration_control import continue_default_membership_through_registration_control
from pymacruise.toads_math import toads_optimized_centroid_similarity, toads_optimized_image_variance
from pymacruise.toads_multichannel_exact_field_refit import toads_multichannel_exact_field_refit
from pymacruise.toads_registration_sampling import vectorized_toads_registration_coefficients, vectorized_toads_registration_energy
from pymacruise.toads_registration import exact_registration_measure_coefficients


ROOT = Path(__file__).parent
CASE_DIR = ROOT / "data/synthetic/toads_p18czo_multichannel_exact_field"
CASE = json.loads((CASE_DIR / "case.json").read_text())
SOURCE_DIR = ROOT / "data/synthetic" / CASE["source_fixture_directory"]
JAVA_REGISTRATION = ROOT / "java_reference/generate_toads_p18czp_anisotropic_registration_reference.sh"
JAVA_CONTROL = ROOT / "java_reference/ToadsP18cyqEnergyControlReferenceHarness.java"
SOURCE_SHA256 = "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
WRAPPER_SHA256 = "b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d"


def _load(path: Path) -> nib.Nifti1Image:
    return nib.load(str(path))


def _bits(value: float) -> str:
    return str(int(np.asarray([value], dtype=np.float32).view(np.int32)[0]))


def _float_bits(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float32).ravel().view(np.int32)


def _dynamic_membership_pass():
    image_volumes = tuple(_load(SOURCE_DIR / row.split(" / ")[0]) for row in CASE["channel_order"])
    field_volumes = tuple(_load(SOURCE_DIR / row.split(" / ")[1]) for row in CASE["channel_order"])
    images = np.stack(tuple(np.asarray(v.dataobj, dtype=np.float32) for v in image_volumes), axis=-1)
    fields = np.stack(tuple(np.asarray(v.dataobj, dtype=np.float32) for v in field_volumes), axis=-1)
    memberships = np.asarray(_load(SOURCE_DIR / "exact_memberships.nii.gz").dataobj, dtype=np.float32)
    segmentation = np.asarray(_load(CASE_DIR / CASE["segmentation"]).dataobj, dtype=np.int16)
    stats = json.loads((SOURCE_DIR / "statistics_inputs.json").read_text())
    params = json.loads((SOURCE_DIR / "membership_parameters.json").read_text())
    centers_before = np.stack((
        np.asarray(stats["initial_centroids"], dtype=np.float32),
        np.asarray([(11 - k) / 12.0 for k in range(11)], dtype=np.float32),
    ))
    exponents = np.repeat(np.asarray(stats["category_exponents"], dtype=np.float32)[None, :], 2, axis=0)
    fresh_fields = toads_multichannel_exact_field_refit(
        images, memberships, segmentation, centers_before, fields
    ).fields
    centers = toads_default_exact_centroid_update(
        images, memberships, centers_before, exponents,
        intensity_field=fresh_fields, smoothness=stats["centroid_smoothness"],
    )
    deviations = toads_default_exact_class_variances(
        images, memberships, centers, segmentation, intensity_field=fresh_fields
    )
    image_variances, _ = toads_optimized_image_variance(
        images, memberships[..., 1:], centers[:, 1:], intensity_field=fresh_fields
    )
    category = np.asarray(params["category_exponents"], dtype=np.float32)
    exponents_full = np.empty((2, 11), dtype=np.float32)
    exponents_full[:, 0] = 1.0
    exponents_full[:, 1:3] = category[0]
    exponents_full[:, 3:8] = category[1]
    exponents_full[:, 8:11] = category[2]
    similarity = toads_optimized_centroid_similarity(
        centers, deviations, image_variances, params["prior_scale"]
    )
    shapes = np.asarray(_load(SOURCE_DIR / "transformed_shapes.nii.gz").dataobj, dtype=np.float32)
    closest = np.asarray(_load(SOURCE_DIR / "membership_closest_labels.nii.gz").dataobj, dtype=np.int8)
    result = toads_default_atlas_membership_pass(
        image_values=images,
        initial_memberships=memberships,
        centroids=centers,
        exponent_by_channel_class=exponents_full,
        modality_weights=np.asarray(params["modality_weights"], dtype=np.float32),
        class_variance_weights=np.asarray(params["class_variance_weights"], dtype=np.float32),
        image_variances=image_variances,
        transformed_shapes=shapes,
        closest_labels=closest,
        smoothing_priors=np.asarray(params["smoothing_priors"], dtype=np.float32),
        centroid_similarity=similarity,
        atlas_shape_priors=np.asarray(params["atlas_shape_priors"], dtype=np.float32),
        smoothing=params["smoothing"],
        prior_factor=params["prior_factor"],
        prior_coefficient=params["prior_coefficient"],
        intensity_field=fresh_fields,
    )
    return result, image_volumes


def _atlas_shapes(dimensions: tuple[int, int, int]) -> np.ndarray:
    gx, gy, gz = np.meshgrid(
        np.linspace(-1.0, 1.0, dimensions[0], dtype=np.float32),
        np.linspace(-1.0, 1.0, dimensions[1], dtype=np.float32),
        np.linspace(-1.0, 1.0, dimensions[2], dtype=np.float32),
        indexing="ij",
    )
    atlas = np.empty((*dimensions, 11), dtype=np.float32)
    for k in range(11):
        cx=np.float32(-0.45+0.09*k)
        cy=np.float32(0.3*np.sin(k))
        cz=np.float32(0.25*np.cos(k))
        radius=np.asarray((gx-cx)**2+1.3*(gy-cy)**2+0.8*(gz-cz)**2,dtype=np.float32)
        atlas[...,k]=np.asarray(np.exp(-radius/np.float32(.42))*np.float32(.75)+np.float32(.05),dtype=np.float32)
    return atlas


def test_dynamic_multichannel_membership_uses_one_shared_rigid_shape_registration(tmp_path: Path) -> None:
    membership = _dynamic_membership_pass()[0]
    dimensions = membership.memberships.shape[:3]
    assert dimensions == (14, 14, 14)
    image_volumes = tuple(_load(SOURCE_DIR / row.split(" / ")[0]) for row in CASE["channel_order"])
    for image in image_volumes:
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        np.testing.assert_allclose(nib.affines.voxel_sizes(image.affine), (0.8, 1.2, 1.7), rtol=0, atol=1e-6)
    atlas = _atlas_shapes(dimensions)
    assert np.any(atlas > 0.0)
    # Membership shape-prior maps in this case remain exactly zero; the
    # nonzero template-space atlas is a separate input to shared registration.
    zero_transformed = np.asarray(_load(SOURCE_DIR / "transformed_shapes.nii.gz").dataobj)
    assert np.count_nonzero(zero_transformed) == 0

    current_energy = np.float32(membership.objective)
    previous_energy = np.float32(current_energy * np.float32(1.01))
    image_res = (0.8, 1.2, 1.7)
    template_res = (1.0, 1.0, 1.0)
    common = dict(
        image_resolutions_mm=image_res,
        template_resolutions_mm=template_res,
        initial_transform=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        previous_energy=previous_energy,
        iteration=3,
        maximum_iterations=5,
        convergence_threshold=1e-4,
        do_growing=False,
        offset=0,
        subsample=2,
        max_registration_iterations=5,
        registration_energy_tolerance=0.0,
        minimum_lambda=1e-3,
    )
    combined = continue_default_membership_through_registration_control(membership, atlas, **common)
    registration = combined.handoff.registration
    assert registration is not None
    assert combined.handoff.membership_pass is membership
    assert combined.stage_order == (
        "computeExactMemberships", "energyRatio", "atlas.registerShapes() -> registerAllShapes()",
        "convergenceCheck", "iterationLimitCheck",
    )

    membership_path, atlas_path = tmp_path / "memberships.bin", tmp_path / "atlas.bin"
    membership_path.write_bytes(np.asarray(membership.memberships, dtype="<f4").tobytes(order="C"))
    atlas_path.write_bytes(np.asarray(atlas, dtype="<f4").tobytes(order="C"))
    completed = subprocess.run(
        ["sh", str(JAVA_REGISTRATION), str(membership_path), str(atlas_path), *(str(d) for d in dimensions),
         *(str(v) for v in image_res), *(str(v) for v in template_res)],
        check=True, capture_output=True, text=True,
    )
    lines = completed.stdout.splitlines()
    java = json.loads(lines[-1])
    assert java["image_resolution_mm"] == list(image_res)
    assert java["template_resolution_mm"] == list(template_res)
    assert java["source"].startswith("VandyCRUISE MultipleToadDeformableAtlas.updateShapeRegistration(rigid)")
    # No image channel axis enters this registration call. Both input channels
    # contribute upstream to the shared class memberships, then one atlas
    # transform is estimated for all registered classes.
    assert membership.memberships.shape == (*dimensions, 11)
    assert registration.transform.shape == (6,)

    number = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
    java_hessians = [
        np.asarray([float(v) for v in re.findall(number, line.split("vector:", 1)[1])], dtype=np.float32)
        for line in lines if line.startswith("H: vector:")
    ]
    java_gradients = [
        np.asarray([float(v) for v in re.findall(number, line.split("vector:", 1)[1])], dtype=np.float32)
        for line in lines if line.startswith("G: vector:")
    ]
    java_trials = [
        tuple(float(value) for value in match.groups())
        for line in lines if (match := re.search(r"a: ([^ (]+)\(([^)]+)\)->\(([^)]+)\)", line))
    ]
    assert len(java_hessians) == len(java_gradients) == len(java_trials) == len(registration.iterations) == 5
    low = np.zeros((11, 3), dtype=np.int32)
    high = np.broadcast_to(np.asarray(dimensions, dtype=np.int32), (11, 3)).copy()
    registered = np.asarray([False] + [True] * 10)
    transform = np.zeros(6, dtype=np.float32)
    for index, (iteration, java_hessian, java_gradient, java_trial) in enumerate(zip(
        registration.iterations, java_hessians, java_gradients, java_trials, strict=True,
    )):
        coefficients = vectorized_toads_registration_coefficients(
            membership.memberships, atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=tuple(float(v) for v in transform),
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
                registered_classes=registered,
                offset=0,
                subsample=2,
                precomputed_affine=True,
        )
        gradient, hessian = exact_registration_measure_coefficients(
            coefficients.gradient, coefficients.hessian,
            coefficients.norm_gradient, coefficients.norm_hessian, coefficients.norm,
            source_faithful=True,
        )
        np.testing.assert_allclose(java_hessian, hessian, rtol=3e-5, atol=3e-6)
        np.testing.assert_allclose(java_gradient, gradient, rtol=3e-5, atol=3e-6)
        energy_before = vectorized_toads_registration_energy(
            membership.memberships, atlas,
            image_resolutions_mm=image_res,
            template_resolutions_mm=template_res,
            transform=tuple(float(v) for v in transform),
            class_bounds_low_xyz=low,
            class_bounds_high_xyz=high,
                registered_classes=registered,
                offset=0,
                subsample=2,
                precomputed_affine=True,
        )
        assert java_trial[0] == pytest.approx(float(np.float32(1.5**index)), rel=2e-6)
        assert java_trial[1] == pytest.approx(float(energy_before), rel=3e-5, abs=3e-6)
        np.testing.assert_allclose(java_trial[2], iteration.energy, rtol=3e-5, atol=3e-6)
        transform = iteration.transform.copy()
    np.testing.assert_allclose(
        np.asarray(java["transform_bits"], dtype=np.int32).view(np.float32),
        registration.transform, rtol=3e-5, atol=3e-6,
    )
    np.testing.assert_allclose(
        np.asarray([java["energy_bits"]], dtype=np.int32).view(np.float32),
        [registration.energy], rtol=3e-5, atol=3e-6,
    )

    ratio_java_source = hashlib.sha256(JAVA_CONTROL.read_bytes()).hexdigest()
    assert ratio_java_source == "d831bae542aaa7ae79d9c741d93aa7fcc3265e9e32bc5dffc88abb30d12dfcdc"
    subprocess.run(["javac", "-d", str(tmp_path), str(JAVA_CONTROL)], check=True, capture_output=True, text=True)
    ratio = combined.control.energy_ratio
    java_control = subprocess.run(
        ["java", "-cp", str(tmp_path), "ToadsP18cyqEnergyControlReferenceHarness",
         _bits(previous_energy), _bits(current_energy), _bits(1e-4), "3", "5", "false"],
        check=True, capture_output=True, text=True,
    )
    control_reference = json.loads(java_control.stdout)
    assert control_reference["ratio_bits"] == int(ratio.view(np.int32))
    assert control_reference["stop"] is combined.control.stop is False
    assert control_reference["run_growing"] is combined.control.run_growing is True
    assert control_reference["call_order"] == [
        "computeExactMemberships", "energyRatio", "registerShapes", "convergenceCheck", "iterationLimitCheck"
    ]

    # The selected source evaluates convergence and maximum-iteration branches
    # after registration; check a max-iteration case against the same Java helper.
    limited = continue_default_membership_through_registration_control(
        membership, atlas, align_shapes=False, **{**common, "iteration": 5}
    )
    limited_java = subprocess.run(
        ["java", "-cp", str(tmp_path), "ToadsP18cyqEnergyControlReferenceHarness",
         _bits(previous_energy), _bits(current_energy), _bits(1e-4), "5", "5", "false"],
        check=True, capture_output=True, text=True,
    )
    limited_reference = json.loads(limited_java.stdout)
    assert limited_reference["stop"] is limited.control.stop is True
    assert limited_reference["iteration_limit_stop"] is True
    assert limited_reference["run_growing"] is limited.control.run_growing is False
    assert limited.stage_order == ("computeExactMemberships", "energyRatio", "convergenceCheck", "iterationLimitCheck")
    assert len(java["bounds"]) == 11
    # VandyCRUISE reports the class bounds as inclusive maximum indices.
    assert java["bounds"][0] == [0, 0, 0, *(dimension - 1 for dimension in dimensions)]
