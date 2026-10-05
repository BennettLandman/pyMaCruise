from __future__ import annotations

import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np
from nibabel.orientations import aff2axcodes

from pymacruise.lesion_standard13_voxel import compute_lesion_toads_standard13_voxel
from pymacruise.lesion_standard13_volume import run_lesion_toads_standard13_volume
from pymacruise.lesion_similarity import compute_lesion_toads_similarity
from pymacruise.lesion_standard13_voxel import STANDARD13_NAMES, STANDARD13_TYPES

ROOT = Path(__file__).parent / "data" / "synthetic" / "lesion_standard13_volume"
HARNESS = Path(__file__).parent / "java_reference" / "LesionToadsStandard13VolumeReferenceHarness.java"
SOURCE_SHA = "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
HARNESS_SHA = "e5350b1988f65f2bb521619952ebd896242e708d503dba48c4c9913236bdc969"
SOURCE_EXCERPT_HASHES = {
    "traversal": "e4ec0a95bd65e7f1452ca80e08e8aecb3bfd2165d2db738738b84091e484b882",
    "face_reads": "b31e2b9108d4d65a3b01dfc2511edd6ae33d82a92bc3caf1adcef17c787a0847",
    "normalization_objective": "22ba167e0ed2b778d9ca1512a7551c08f0b13516a9e1d8adda7f8d413c49b033",
}
SOURCE_EXCERPTS = {
    "traversal": ROOT / "method_volume_traversal_source_excerpt.java.txt",
    "face_reads": ROOT / "method_face_reads_source_excerpt.java.txt",
    "normalization_objective": ROOT / "method_normalization_objective_source_excerpt.java.txt",
}
CLASSES = 13


def _state() -> dict[str, np.ndarray | float]:
    image = np.asarray(nib.load(str(ROOT / "folded_inputs_ras.nii.gz")).dataobj, dtype=np.float32)
    x, y, z = np.indices(image.shape[:3], dtype=np.int32)
    k = np.arange(CLASSES, dtype=np.float32)
    centers = np.stack((np.float32(0.2) + np.float32(0.04) * k,
                        np.float32(0.3) + np.float32(0.03) * k))
    deviations = np.stack((np.float32(0.3) + np.float32(0.01) * k,
                           np.float32(0.4) + np.float32(0.005) * k))
    shape = np.empty((*image.shape[:3], CLASSES), dtype=np.float32)
    direct = np.empty(shape.shape, dtype=np.bool_)
    for ix in range(5):
        for iy in range(5):
            for iz in range(5):
                for cls in range(CLASSES):
                    spatial_offset = np.float32(np.float32(0.001) * np.float32((ix + 2 * iy + 3 * iz) % 7))
                    class_offset = np.float32(np.float32(0.02) * np.float32(cls))
                    shape[ix, iy, iz, cls] = np.float32(np.float32(np.float32(0.4) + class_offset) + spatial_offset)
                    direct[ix, iy, iz, cls] = (ix + 2 * iy + 3 * iz + cls) % 2 == 0
    vent = np.where((x + y + z) % 4 == 0, -1, (x + 2 * y + z) % 3).astype(np.int16)
    gray = np.full(image.shape[:3], -1, dtype=np.int16)  # GMflag is false on the first call.
    bstem = np.where((x + z) % 5 == 0, -1, (2 * x + y + z) % 26).astype(np.int16)
    return {
        "images": image,
        "class_centroids": centers,
        "class_stddev": deviations,
        "lesion_centroid": np.asarray([0.5, 0.65], dtype=np.float32),
        "lesion_stddev": np.asarray([0.25, 0.3], dtype=np.float32),
        "channel_weights": np.asarray([1.0, 0.75], dtype=np.float32),
        "shape_values": shape,
        "shape_priors": np.float32(0.5) + np.float32(0.05) * k,
        "smoothing_priors": np.float32(0.8) + np.float32(0.025) * k,
        "direct_relations": direct,
        "initial_class_memberships": np.zeros((*image.shape[:3], CLASSES), dtype=np.float32),
        "initial_lesion_membership": np.zeros(image.shape[:3], dtype=np.float32),
        "initial_black_hole_membership": np.zeros(image.shape[:3], dtype=np.float32),
        "initial_outlier_membership": np.zeros(image.shape[:3], dtype=np.float32),
        "initial_lesion_relation": np.zeros(image.shape[:3], dtype=np.float32),
        "ventricle_distances": vent,
        "gray_matter_distances": gray,
        "brainstem_distances": bstem,
        "outlier_energy": np.float32(0.1),
        "lesion_power_exponent": np.float32(1.2),
    }


def _flat_bits(value: np.ndarray | np.float32) -> list[int] | int:
    values = np.asarray(value, dtype=np.float32)
    if values.ndim == 0:
        return int(values.view(np.uint32))
    return values.view(np.uint32).reshape(-1).astype(np.uint64).tolist()


def test_standard13_volume_matches_source_order_java_harness() -> None:
    image = nib.load(str(ROOT / "folded_inputs_ras.nii.gz"))
    assert image.shape == (5, 5, 5, 2)
    assert aff2axcodes(image.affine) == ("R", "A", "S")
    np.testing.assert_array_equal(image.affine, np.eye(4))
    values = np.asarray(image.dataobj, dtype=np.float32)
    image_hash = hashlib.sha256(values.astype("<f4", copy=False).tobytes()).hexdigest()
    reference = json.loads((ROOT / "java_reference.json").read_text())
    assert reference["source_sha256"] == SOURCE_SHA
    assert reference["image_sha256"] == image_hash
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    assert reference["source_excerpt_sha256"] == SOURCE_EXCERPT_HASHES
    for name, path in SOURCE_EXCERPTS.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == SOURCE_EXCERPT_HASHES[name]
    traversal = SOURCE_EXCERPTS["traversal"].read_text()
    assert "for (int x=1;x<nx-1;x++) for (int y=1;y<ny-1;y++) for (int z=1;z<nz-1;z++)" in traversal
    face_reads = SOURCE_EXCERPTS["face_reads"].read_text()
    face_markers = ["// case by case\t: X+", "// case by case\t: X-", "// case by case\t: Y+",
                    "// case by case\t: Y-", "// case by case\t: Z+", "// case by case\t: Z-"]
    assert [face_reads.index(marker) for marker in face_markers] == sorted(
        face_reads.index(marker) for marker in face_markers
    )
    normalization = SOURCE_EXCERPTS["normalization_objective"].read_text()
    assert "mems[x][y][z][k] = mems[x][y][z][k]/den;" in normalization
    assert "dist += energy[k]*mems[x][y][z][k]*mems[x][y][z][k];" in normalization

    # The entire RAS NIfTI must match the in-memory image formula used by the
    # independent Java harness, including boundary voxels that the update loop
    # leaves untouched.
    x, y, z = np.indices(values.shape[:3], dtype=np.int32)
    ti = (x * 7 + y * 5 + z * 3 + ((x * y) % 3) * 2) % 16
    fi = (x * 3 + y * 7 + z * 5 + ((y * z) % 4) * 3) % 15
    expected_image = np.empty_like(values)
    expected_image[..., 0] = np.float32(0.4) + np.float32(0.025) * ti.astype(np.float32)
    expected_image[..., 1] = np.float32(0.55) + np.float32(0.02) * fi.astype(np.float32)
    np.testing.assert_array_equal(values, expected_image)
    state = _state()
    state_before = {
        name: value.copy() for name, value in state.items() if isinstance(value, np.ndarray)
    }
    result = run_lesion_toads_standard13_volume(**state)
    for name, before in state_before.items():
        np.testing.assert_array_equal(state[name], before, err_msg=f"mutated input: {name}")
    actual = {
        "class_membership_bits": _flat_bits(result.class_memberships),
        "lesion_membership_bits": _flat_bits(result.lesion_membership),
        "black_hole_membership_bits": _flat_bits(result.black_hole_membership),
        "outlier_membership_bits": _flat_bits(result.outlier_membership),
        "lesion_relation_bits": _flat_bits(result.lesion_relation),
        "denominator_bits": _flat_bits(result.denominator),
        "class_energy_bits": _flat_bits(result.class_energies),
        "lesion_energy_bits": _flat_bits(result.lesion_energy),
        "black_hole_energy_bits": _flat_bits(result.black_hole_energy),
        "running_objective_bits": _flat_bits(result.running_objective),
        "objective_bits": _flat_bits(result.objective),
    }
    similarity = compute_lesion_toads_similarity(
        state["class_centroids"], state["class_stddev"], state["lesion_centroid"],
        state["lesion_stddev"], state["channel_weights"],
        tuple(kind.lower() for kind in STANDARD13_TYPES), STANDARD13_NAMES,
        prior_scale=0.1, t1_channel=0,
        black_hole_similarity=np.zeros(CLASSES, dtype=np.float32),
    )
    actual["class_similarity_bits"] = _flat_bits(similarity.class_similarity)
    actual["lesion_similarity_bits"] = _flat_bits(similarity.lesion_similarity)
    for name, bits in actual.items():
        assert bits == reference[name], name

    # The first two interior voxels are Z-face neighbors in source traversal.
    # The second receives the first voxel's updated memberships, while all
    # six source boundaries retain their supplied zeros.
    initial_faces = np.zeros((6, CLASSES), dtype=np.float32)
    independent = compute_lesion_toads_standard13_voxel(
        state["images"][1, 1, 2], state["class_centroids"], state["class_stddev"],
        state["lesion_centroid"], state["lesion_stddev"], state["channel_weights"],
        state["shape_values"][1, 1, 2], state["shape_priors"], state["smoothing_priors"],
        state["direct_relations"][1, 1, 2], np.zeros(CLASSES, dtype=np.float32), 0, 0, 0, 0,
        initial_faces, np.zeros(6, dtype=np.float32), np.zeros(6, dtype=np.float32),
        outlier_energy=state["outlier_energy"], lesion_power_exponent=state["lesion_power_exponent"],
        ventricle_distance=int(state["ventricle_distances"][1, 1, 2]),
        gray_matter_distance=-1,
        brainstem_distance=int(state["brainstem_distances"][1, 1, 2]),
    )
    assert not np.array_equal(result.class_memberships[1, 1, 2], independent.class_memberships)
    dependent_faces = np.zeros((6, CLASSES), dtype=np.float32)
    dependent_faces[5] = result.class_memberships[1, 1, 1]
    dependent_lesion_faces = np.zeros(6, dtype=np.float32)
    dependent_lesion_faces[5] = result.lesion_membership[1, 1, 1]
    dependent = compute_lesion_toads_standard13_voxel(
        state["images"][1, 1, 2], state["class_centroids"], state["class_stddev"],
        state["lesion_centroid"], state["lesion_stddev"], state["channel_weights"],
        state["shape_values"][1, 1, 2], state["shape_priors"], state["smoothing_priors"],
        state["direct_relations"][1, 1, 2], np.zeros(CLASSES, dtype=np.float32), 0, 0, 0, 0,
        dependent_faces, dependent_lesion_faces, np.zeros(6, dtype=np.float32),
        outlier_energy=state["outlier_energy"], lesion_power_exponent=state["lesion_power_exponent"],
        ventricle_distance=int(state["ventricle_distances"][1, 1, 2]),
        gray_matter_distance=-1,
        brainstem_distance=int(state["brainstem_distances"][1, 1, 2]),
    )
    np.testing.assert_array_equal(result.class_memberships[1, 1, 2], dependent.class_memberships)
    assert result.lesion_membership[1, 1, 2] == dependent.lesion_membership
    for axis in range(3):
        for edge in (0, -1):
            index = [slice(None)] * 3
            index[axis] = edge
            index_tuple = tuple(index)
            assert not np.any(result.denominator[index_tuple])
            np.testing.assert_array_equal(result.class_memberships[index_tuple], state["initial_class_memberships"][index_tuple])
            for name in (
                "lesion_membership", "black_hole_membership", "outlier_membership",
                "lesion_relation", "class_energies", "lesion_energy",
                "black_hole_energy", "running_objective",
            ):
                assert not np.any(getattr(result, name)[index_tuple]), name


def test_standard13_volume_applies_channel_correction_fields() -> None:
    state = _state()
    images = np.asarray(state["images"], dtype=np.float32)
    uncorrected = run_lesion_toads_standard13_volume(**state)
    identity_fields = np.ones_like(images, dtype=np.float32)
    identity = run_lesion_toads_standard13_volume(**state, intensity_fields=identity_fields)
    for name in (
        "class_memberships", "lesion_membership", "black_hole_membership",
        "outlier_membership", "lesion_relation", "denominator", "class_energies",
        "lesion_energy", "black_hole_energy", "running_objective",
    ):
        np.testing.assert_array_equal(getattr(identity, name), getattr(uncorrected, name), err_msg=name)
    assert identity.objective.view(np.uint32) == uncorrected.objective.view(np.uint32)

    fields = np.empty_like(images)
    fields[..., 0] = np.float32(1.1)
    fields[..., 1] = np.float32(0.9)
    corrected = run_lesion_toads_standard13_volume(**state, intensity_fields=fields)
    assert corrected.objective.view(np.uint32) != uncorrected.objective.view(np.uint32)
    assert np.count_nonzero(corrected.class_memberships != uncorrected.class_memberships) > 0


def test_standard13_volume_rejects_invalid_correction_field() -> None:
    state = _state()
    images = np.asarray(state["images"], dtype=np.float32)
    for fields in (
        np.ones((*images.shape[:3], 1), dtype=np.float32),
        np.full_like(images, np.nan),
    ):
        try:
            run_lesion_toads_standard13_volume(**state, intensity_fields=fields)
        except ValueError as error:
            assert "intensity_fields" in str(error)
        else:  # pragma: no cover - assertion branch
            raise AssertionError("invalid intensity fields were accepted")
