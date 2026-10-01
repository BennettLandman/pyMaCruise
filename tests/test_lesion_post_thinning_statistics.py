from __future__ import annotations

import json
import hashlib
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_initial_loop_entry import seed_lesion_toads_initial_loop
from pymacruise.lesion_initial_thinning_composition import run_lesion_initial_thinning_pass
from pymacruise.lesion_post_thinning_statistics import refresh_statistics_after_initial_thinning
from pymacruise.lesion_standard13_volume import run_lesion_toads_standard13_volume


ROOT = Path(__file__).parent / "data/synthetic/lesion_initial_thinning_composition"
SPEC = json.loads((ROOT / "fixture.json").read_text(encoding="utf-8"))
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
STATISTICS_REFERENCE = json.loads(
    (Path(__file__).parent / "data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json").read_text(encoding="utf-8")
)
STATISTICS_HARNESS = Path(__file__).parent / "java_reference/LesionPostThinningStatisticsReferenceHarness.java"


def _bits(values: np.ndarray) -> list[int] | list[list[int]]:
    return values.astype(np.float32, copy=False).view(np.uint32).astype(np.uint64).tolist()


def test_first_thinning_output_feeds_selected_statistics_refresh() -> None:
    image = nib.load(str(ROOT / "initial_thinning_inputs_ras.nii.gz"))
    frames = np.asarray(image.dataobj, dtype=np.float32)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    state = {
        "images": frames[..., 0:2],
        "shape_values": frames[..., 2:15],
        "direct_relations": frames[..., 15:28] > 0.5,
        "segmentation": frames[..., 28].astype(np.int16),
        "skeleton": frames[..., 29].astype(np.int16),
        "available": frames[..., 30:43] > 0.5,
        "ordering": frames[..., 43],
        "ventricle_distances": frames[..., 44].astype(np.int16),
        "gray_matter_distances": frames[..., 45].astype(np.int16),
        "brainstem_distances": frames[..., 46].astype(np.int16),
        "memberships": frames[..., 47:60],
        "lesion": frames[..., 60],
        "black": frames[..., 61],
        "outlier": frames[..., 62],
        "lesion_relation": frames[..., 63],
    }
    membership = run_lesion_toads_standard13_volume(
        images=state["images"],
        class_centroids=SPEC["class_centroids"],
        class_stddev=SPEC["class_stddev"],
        lesion_centroid=SPEC["lesion_centroid"],
        lesion_stddev=SPEC["lesion_stddev"],
        channel_weights=SPEC["channel_weights"],
        shape_values=state["shape_values"],
        shape_priors=SPEC["shape_priors"],
        smoothing_priors=SPEC["smoothing_priors"],
        direct_relations=state["direct_relations"],
        initial_class_memberships=state["memberships"],
        initial_lesion_membership=state["lesion"],
        initial_black_hole_membership=state["black"],
        initial_outlier_membership=state["outlier"],
        initial_lesion_relation=state["lesion_relation"],
        ventricle_distances=state["ventricle_distances"],
        gray_matter_distances=state["gray_matter_distances"],
        brainstem_distances=state["brainstem_distances"],
        outlier_energy=SPEC["outlier_energy"],
        lesion_power_exponent=SPEC["lesion_power_exponent"],
        gm_flag=False,
    )
    loop_entry = seed_lesion_toads_initial_loop(membership, iterations=1)
    thinning = run_lesion_initial_thinning_pass(
        loop_entry,
        state["segmentation"],
        SPEC["object_types"],
        state["shape_values"],
        state["skeleton"],
        state["available"],
        state["ordering"],
        lesion_class=SPEC["lesion_class"],
        volume_ratio=SPEC["volume_ratio"],
        limit=SPEC["thinning_limit"],
    )
    # This exact folded phantom's initial membership and final first-thinning
    # segmentation are pinned against the selected Java composition harness.
    assert membership.class_memberships.astype(np.float32).view(np.uint32).reshape(-1).tolist() == REFERENCE["class_membership_bits"]
    np.testing.assert_array_equal(
        thinning.thinning.run.segmentation.ravel(), REFERENCE["final_segmentation"]
    )
    assert hashlib.sha256(STATISTICS_HARNESS.read_bytes()).hexdigest() == STATISTICS_REFERENCE["harness_sha256"]
    assert STATISTICS_REFERENCE["source_sha256"] == REFERENCE["source_sha256"]
    assert STATISTICS_REFERENCE["wrapper_sha256"] == REFERENCE["wrapper_sha256"]
    assert STATISTICS_REFERENCE["power_table_file_sha256"] == "2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae"
    fixed = tuple(a.copy() for a in (state["images"], membership.class_memberships, state["lesion"], state["black"], state["outlier"]))

    refreshed = refresh_statistics_after_initial_thinning(
        thinning,
        state["images"],
        SPEC["class_centroids"],
        state["lesion"],
        state["black"],
        state["outlier"],
        SPEC["lesion_centroid"],
        (0.0, 0.0),  # explicit fixed synthetic initializer; black-hole map is zero
        SPEC["object_types"],
        centroid_options={
            "mode": "none",
            "modalities": ("T1_MPRAGE", "FLAIR"),
            "lesion_power_exponents": (0.4, 1.2),
        },
        channel_weight_options={"outlier_distance": SPEC["outlier_energy"]},
        flair_channel=1,
    )

    assert refreshed.stage_order[-3:] == (
        "computeCentroids", "computeVariances", "computeChannelWeights"
    )
    assert _bits(refreshed.centroids.class_centroids) == STATISTICS_REFERENCE["class"]
    assert refreshed.centroids.lesion_centroids is not None
    assert _bits(refreshed.centroids.lesion_centroids) == STATISTICS_REFERENCE["aux_centroids"][0]
    assert refreshed.centroids.black_hole_centroids is not None
    assert _bits(refreshed.centroids.black_hole_centroids) == STATISTICS_REFERENCE["aux_centroids"][1]
    assert _bits(refreshed.variances.class_standard_deviation.reshape(-1)) == STATISTICS_REFERENCE["class_std"]
    assert _bits(refreshed.variances.lesion_standard_deviation) == STATISTICS_REFERENCE["lesion_std"]
    assert _bits(refreshed.variances.black_hole_standard_deviation) == STATISTICS_REFERENCE["blackhole_std"]
    assert _bits(refreshed.channel_weights) == STATISTICS_REFERENCE["weights"]
    assert refreshed.centroids.class_centroids.shape == (2, 13)
    assert refreshed.centroids.lesion_centroids is not None
    assert refreshed.centroids.black_hole_centroids is not None
    np.testing.assert_array_equal(refreshed.centroids.black_hole_centroids, (0.0, 0.0))
    assert refreshed.variances.class_standard_deviation.shape == (2, 13)
    assert refreshed.variances.lesion_standard_deviation.shape == (2,)
    assert refreshed.channel_weights.shape == (2,)
    assert np.isfinite(refreshed.centroids.class_centroids).all()
    assert np.isfinite(refreshed.centroids.lesion_centroids).all()
    assert np.isfinite(refreshed.variances.class_standard_deviation).all()
    assert np.isfinite(refreshed.variances.lesion_standard_deviation).all()
    assert np.isfinite(refreshed.channel_weights).all()
    np.testing.assert_allclose(refreshed.channel_weights.sum(), np.float32(1.0), rtol=0, atol=2e-7)
    for actual, before in zip(
        (state["images"], membership.class_memberships, state["lesion"], state["black"], state["outlier"]),
        fixed,
    ):
        np.testing.assert_array_equal(actual, before)
