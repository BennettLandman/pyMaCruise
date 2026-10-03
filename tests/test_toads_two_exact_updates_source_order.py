from __future__ import annotations

import json
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_control import toads_energy_ratio
from pymacruise.toads_initial_thinning_composition import (
    toads_initial_fielded_membership_to_thinning,
)
from pymacruise.toads_initial_thinning_statistics_membership import (
    continue_initial_thinning_through_memberships,
)


ROOT = Path(__file__).parent / "data" / "synthetic" / "toads_initial_thinning_startup"
JAVA_RUNNER = Path(__file__).parent.parent / "scripts" / "generate_java_toads_two_exact_updates_reference.sh"


def _load(name: str) -> np.ndarray:
    image = nib.load(str(ROOT / name))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    return np.asarray(image.dataobj)


def _bits(value: np.ndarray | float) -> int | list[int] | list[list[int]]:
    values = np.asarray(value, dtype=np.float32).view(np.uint32)
    if values.ndim == 0:
        return int(values)
    return values.tolist()


def test_two_source_ordered_exact_memberships_use_fixed_atlas_factors() -> None:
    spec = json.loads((ROOT / "spec.json").read_text())
    java = json.loads((ROOT / "two_exact_updates_java_reference.json").read_text())
    generated = subprocess.run(
        ["sh", str(JAVA_RUNNER)], check=True, capture_output=True, text=True
    )
    assert json.loads(generated.stdout) == java
    params = spec["parameters"]
    images = _load("images_ras.nii.gz").astype(np.float32)
    labels = _load("segmentation_ras.nii.gz").astype(np.uint8)
    initial = _load("initial_memberships_ras.nii.gz").astype(np.float32)
    fields = _load("fields_ras.nii.gz").astype(np.float32)
    bias = _load("psf_bias_ras.nii.gz").astype(np.float32)
    shapes = _load("transformed_shapes_ras.nii.gz").astype(np.float32)
    centers = np.asarray(params["centroids"], dtype=np.float32)
    exponents = np.asarray(params["exponent_by_channel_class"], dtype=np.float32)

    first = toads_initial_fielded_membership_to_thinning(
        image_values=images,
        segmentation=labels,
        initial_memberships=initial,
        centroids=centers,
        exponent_by_channel_class=exponents,
        modality_weights=params["modality_weights"],
        class_variance_weights=params["class_variance_weights"],
        image_variances=params["image_variances"],
        transformed_shapes=shapes,
        smoothing_priors=params["smoothing_priors"],
        centroid_similarity=params["centroid_similarity"],
        atlas_shape_priors=params["atlas_shape_priors"],
        smoothing=spec["membership_smoothing"],
        prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"],
        channel_fields=[fields[..., 0], fields[..., 1]],
        limit=spec["relation_limit"],
        volume_ratio=spec["thinning_volume_ratio"],
        thinning_limit=spec["thinning_limit"],
        psf_bias=bias,
    )
    assert _bits(first.membership_pass.objective) == java["first_objective_bits"]

    category_exponents = np.asarray(
        [exponents[:, 1], exponents[:, 3], exponents[:, 8]], dtype=np.float32
    ).T
    second = continue_initial_thinning_through_memberships(
        first,
        image_values=images,
        initial_centroids=centers,
        category_exponents=category_exponents,
        transformed_shapes=shapes,
        image_variance_weights=params["image_variances"],
        class_variance_weights=params["class_variance_weights"],
        smoothing_priors=params["smoothing_priors"],
        atlas_shape_priors=params["atlas_shape_priors"],
        exponent_by_channel_class=exponents,
        modality_weights=params["modality_weights"],
        prior_factor=spec["prior_factor"],
        prior_coefficient=spec["prior_coefficient"],
        smoothing=spec["membership_smoothing"],
        prior_scale=0.05,
        intensity_field=fields,
        psf_bias=bias,
    )

    assert java["source_order"] == [
        "first computeExactMemberships",
        "computeExactCentroids",
        "computeExactVariances",
        "computeCentroidSimilarity",
        "second computeExactMemberships",
        "energyRatio",
    ]
    assert _bits(second.membership_pass.objective) == java["second_objective_bits"]
    assert _bits(toads_energy_ratio(first.membership_pass.objective, second.membership_pass.objective)) == java["energy_ratio_bits"]
    assert _bits(params["class_variance_weights"]) == java["fixed_intensity_variance_factor_bits"]
    assert _bits(second.class_standard_deviations) == java["refreshed_standard_deviation_bits"]
    assert not np.array_equal(
        second.class_standard_deviations,
        np.asarray(params["class_variance_weights"], dtype=np.float32),
    )
