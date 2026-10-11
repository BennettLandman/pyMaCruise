from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from compare_p41z_next_iteration_statistics_java import compare


TYPES = ["mask", "optimized_csf", "optimized_csf", "object", "object", "object", "object", "object", "object", "object", "lesion", "object", "object"]


def _save(path: Path, values: np.ndarray) -> None:
    # A diagonal RAS affine makes the contract explicit and keeps image
    # storage order identical to the Java [x][y][z][channel] traversal.
    affine = np.diag([1.0, 1.0, 1.0, 1.0])
    nib.save(nib.Nifti1Image(values, affine), str(path))


def _write_fixture(root: Path, *, zero_energy: bool = False) -> None:
    root.mkdir(parents=True)
    shape = (4, 3, 2)
    classes = len(TYPES)
    xyz = np.indices(shape, dtype=np.float32)
    seg = ((xyz[0] * 7 + xyz[1] * 5 + xyz[2] * 3) % classes).astype(np.int16)
    if zero_energy:
        # All OBJECT-class residuals and the LESION linear residuals vanish
        # after their centroid updates; this exercises the source's INF
        # sentinel branch without treating a finite sentinel as IEEE infinity.
        image = np.zeros((*shape, 2), dtype=np.float32)
    else:
        image = np.empty((*shape, 2), dtype=np.float32)
        image[..., 0] = 0.15 + 0.017 * xyz[0] + 0.031 * xyz[1] + 0.009 * xyz[2]
        image[..., 1] = 0.23 + 0.011 * xyz[0] + 0.019 * xyz[1] + 0.027 * xyz[2]
    memberships = np.empty((*shape, classes), dtype=np.float32)
    for k in range(classes):
        memberships[..., k] = 0.12 + ((xyz[0] * 3 + xyz[1] * 5 + xyz[2] * 7 + k * 2) % 9) * 0.045
    lesion = (0.17 + ((xyz[0] + 2 * xyz[1] + xyz[2]) % 7) * 0.07).astype(np.float32)
    black = (0.09 + ((2 * xyz[0] + xyz[1] + 3 * xyz[2]) % 6) * 0.055).astype(np.float32)
    outlier = (0.03 + ((xyz[0] + xyz[1] + 2 * xyz[2]) % 5) * 0.035).astype(np.float32)
    if zero_energy:
        # Null channel-weight contributors while retaining valid
        # centroids/variance inputs. Source uses a finite 1e30 sentinel.
        memberships[..., 3:] = 0.0
        lesion[:] = 0.0
        black[:] = 0.0
    initial_centroids = np.full((2, classes), 0.5, dtype=np.float32)
    initial_centroids[:, 0] = 0.0
    initial_lesion = np.asarray([0.4, 0.4], dtype=np.float32)
    initial_black = np.asarray([0.3, 0.3], dtype=np.float32)
    tissue_power = np.ones((2, classes), dtype=np.float32)
    lesion_power = np.asarray([0.4, 1.2], dtype=np.float32)
    case = {
        "modalities": ["T1_MPRAGE", "FLAIR"],
        "object_types": TYPES,
        "initial_class_centroids": initial_centroids.tolist(),
        "initial_lesion_centroids": initial_lesion.tolist(),
        "initial_black_hole_centroids": initial_black.tolist(),
        "tissue_power_exponents": tissue_power.tolist(),
        "lesion_power_exponents": lesion_power.tolist(),
        "norm_type": "Identity",
        "outlier_distance": 0.1,
        "optimized_csf_all_modalities": True,
        "zero_energy_case": zero_energy,
    }
    (root / "case.json").write_text(json.dumps(case, indent=2) + "\n", encoding="utf-8")
    _save(root / "images.nii.gz", image)
    _save(root / "memberships.nii.gz", memberships)
    _save(root / "segmentation.nii.gz", seg)
    _save(root / "lesion_membership.nii.gz", lesion)
    _save(root / "black_hole_membership.nii.gz", black)
    _save(root / "outlier_membership.nii.gz", outlier)


@pytest.mark.parametrize("zero_energy", [False, True])
def test_p41z_next_iteration_statistics_java_io_and_source_nan_report(tmp_path: Path, zero_energy: bool) -> None:
    if shutil.which("java") is None or shutil.which("javac") is None:
        pytest.skip("Java runtime and compiler are required for this selected-source comparison")
    fixture = tmp_path / ("zero-energy" if zero_energy else "folded-synthetic")
    _write_fixture(fixture, zero_energy=zero_energy)
    report = compare(fixture)
    assert report["stage_order"] == [
        "propagateThinning", "computeCentroids", "computeVariances", "computeChannelWeights"
    ]
    assert report["norm_type"] == "Identity"
    assert report["outlier_distance"] == 0.1
    assert report["python_channel_weight_bits"] == [0x3F000000, 0x3F000000]
    if zero_energy:
        # The comparator counts NaNs explicitly; this input exercises the
        # zero-energy fallback, which remains finite for this two-channel case.
        assert report["nan_counts"] == {"java": 0, "python": 0}
        assert report["all_exact"], report
    else:
        assert report["nan_counts"] == {"java": 0, "python": 0}
        assert report["all_exact"], report
        # Regression guard: the current VandyCRUISE profile applies the
        # optimized-CSF PowerTable on FLAIR too. Keep the legacy default
        # available for older selected-source comparisons, but prove this
        # exact P41z phantom distinguishes the two policies.
        if not zero_energy:
            from pymacruise.lesion_statistics import compute_lesion_toads_centroids

            images = np.asarray(nib.load(str(fixture / "images.nii.gz")).dataobj, dtype=np.float32)
            memberships = np.asarray(nib.load(str(fixture / "memberships.nii.gz")).dataobj, dtype=np.float32)
            labels = np.asarray(nib.load(str(fixture / "segmentation.nii.gz")).dataobj)
            case = json.loads((fixture / "case.json").read_text(encoding="utf-8"))
            base_kwargs = {
                "outlier_membership": np.asarray(nib.load(str(fixture / "outlier_membership.nii.gz")).dataobj, dtype=np.float32),
                "lesion_membership": np.asarray(nib.load(str(fixture / "lesion_membership.nii.gz")).dataobj, dtype=np.float32),
                "black_hole_membership": np.asarray(nib.load(str(fixture / "black_hole_membership.nii.gz")).dataobj, dtype=np.float32),
                "lesion_centroids": np.asarray(case["initial_lesion_centroids"], dtype=np.float32),
                "black_hole_centroids": np.asarray(case["initial_black_hole_centroids"], dtype=np.float32),
                "mode": "none",
                # Give FLAIR a non-unit diagnostic exponent solely to make
                # the policy switch visible on this tiny fixture; the Java
                # and full-volume comparisons below still use P41z's exact
                # unit-exponent production profile.
                "tissue_power_exponents": np.asarray(case["tissue_power_exponents"], dtype=np.float32),
                "modalities": tuple(case["modalities"]),
                "lesion_power_exponents": np.asarray(case["lesion_power_exponents"], dtype=np.float32),
            }
            diagnostic_powers = base_kwargs["tissue_power_exponents"].copy()
            diagnostic_powers[1, (1, 2)] = np.float32(0.65)
            diagnostic_kwargs = {**base_kwargs, "tissue_power_exponents": diagnostic_powers}
            legacy = compute_lesion_toads_centroids(
                images, memberships, labels,
                np.asarray(case["initial_class_centroids"], dtype=np.float32),
                TYPES, optimized_csf_all_modalities=False, **diagnostic_kwargs,
            )
            active_profile = compute_lesion_toads_centroids(
                images, memberships, labels,
                np.asarray(case["initial_class_centroids"], dtype=np.float32),
                TYPES, optimized_csf_all_modalities=True, **diagnostic_kwargs,
            )
            flair_differences = [
                legacy.class_centroids[1, class_index].view(np.uint32)
                != active_profile.class_centroids[1, class_index].view(np.uint32)
                for class_index in (1, 2)
            ]
            assert any(flair_differences), "FLAIR optimized-CSF power-table switch had no effect on this phantom"


@pytest.mark.skipif(
    os.environ.get("PYMACRUISE_RUN_P41Z_ITER6_STATS") != "1",
    reason="set PYMACRUISE_RUN_P41Z_ITER6_STATS=1 for full-volume held-fixed comparison",
)
def test_p41z_iteration6_held_fixed_statistics_match_java() -> None:
    if shutil.which("java") is None or shutil.which("javac") is None:
        pytest.skip("Java runtime and compiler are required for this selected-source comparison")
    fixture = Path(__file__).parent / "data/synthetic/p41z_iteration6_statistics_java"
    report = compare(fixture)
    assert report["spatial_shape_xyz"] == [71, 89, 67]
    assert report["norm_type"] == "Identity"
    assert report["outlier_distance"] == 0.1
    assert report["stage_order"] == [
        "propagateThinning", "computeCentroids", "computeVariances", "computeChannelWeights"
    ]
    assert report["nan_counts"]["java"] == report["nan_counts"]["python"]
    assert report["all_exact"], report
