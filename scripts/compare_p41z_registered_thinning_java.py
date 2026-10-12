#!/usr/bin/env python3
"""Compare the next registered P41z thinning handoff with selected Java.

The registered checkpoint, transform, and support bounds are held fixed. The
script compares selected Java atlas-shape sampling and 13-class next thinning
against the Python production adapters. It does not run the full next loop.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile

import nibabel as nib
import numpy as np

from pymacruise.lesion_atlas_priors import (
    load_lesion_atlas_priors,
    load_lesion_atlas_shape_priors,
    load_lesion_atlas_topology,
)
from pymacruise.lesion_next_iteration_membership import run_lesion_toads_next_gm_membership_update
from pymacruise.lesion_next_iteration_statistics import refresh_statistics_after_next_iteration_thinning
from pymacruise.lesion_distance_refresh import refresh_lesion_toads_distance_maps
from pymacruise.lesion_registered_shapes import transform_registered_lesion_shape_priors
from pymacruise.lesion_registered_initial_membership import LesionInitialDistanceMaps
from pymacruise.lesion_standard13_volume import STANDARD13_TYPES
from pymacruise.lesion_growth import propagate_lesion_toads_growing
from pymacruise.plus_lesion_registered_loop import hydrate_plus_lesion_registered_branch
from pymacruise.lesion_next_iteration_thinning import run_lesion_toads_next_iteration_thinning


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from compare_p41z_next_iteration_statistics_java import compare as compare_statistics  # noqa: E402

JAVA_ROOT = ROOT / "tests/java_reference"
UPSTREAM = JAVA_ROOT / "upstream/edu/jhmi/rad/medic"
TYPE_CODES = {"lesion": 1, "blackhole": 2}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _checked_path(record: dict, *, name: str) -> Path:
    path = Path(record.get("path", "")).resolve()
    if not path.is_file() or path.is_symlink() or _sha256(path) != record.get("sha256"):
        raise ValueError(f"{name} is missing or fails its checkpoint hash")
    return path


def _item_bytes(item) -> bytes:
    return struct.pack(">fiiii", float(item.value), item.x, item.y, item.z, item.label)


def _python_thinning_bytes(thinning) -> bytes:
    output = bytearray()
    output.extend(np.asarray(thinning.initialization.class_volume, dtype=">f4").tobytes())
    output.extend(np.asarray(thinning.minimum_class_volume, dtype=">f4").tobytes())
    output.extend(struct.pack(">i", len(thinning.initialization.inserted_items)))
    for item in thinning.initialization.inserted_items:
        output.extend(_item_bytes(item))
    output.extend(struct.pack(">i", len(thinning.run.events)))
    for event in thinning.run.events:
        candidate = -1 if event.candidate_label is None else int(event.candidate_label)
        topology = -1 if event.topology_accepted is None else int(event.topology_accepted)
        output.extend(_item_bytes(event.popped_item))
        output.extend(struct.pack(">bibi", int(event.processed), candidate, topology, len(event.inserted_neighbors)))
        for item in event.inserted_neighbors:
            output.extend(_item_bytes(item))
    output.extend(np.asarray(thinning.run.segmentation, dtype=np.int8).ravel(order="C").tobytes())
    output.extend(np.asarray(thinning.run.skeleton, dtype=np.int8).ravel(order="C").tobytes())
    output.extend(np.asarray(thinning.run.available, dtype=np.uint8).ravel(order="C").tobytes())
    output.extend(np.asarray(thinning.run.ordering, dtype=">f4").ravel(order="C").tobytes())
    output.extend(np.asarray(thinning.run.class_volume, dtype=">f4").tobytes())
    remaining = []
    while thinning.run.queue:
        remaining.append(thinning.run.queue.pop_first())
    output.extend(struct.pack(">i", len(remaining)))
    for item in remaining:
        output.extend(_item_bytes(item))
    return bytes(output)


def _membership_input_bytes(
    spatial: tuple[int, int, int],
    image: np.ndarray,
    shapes: np.ndarray,
    direct: np.ndarray,
    memberships: np.ndarray,
    lesion: np.ndarray,
    black: np.ndarray,
    outlier: np.ndarray,
    lesion_relation: np.ndarray,
    distances: np.ndarray,
    statistics,
    previous_energy: float,
) -> bytes:
    record_type = np.dtype([
        ("image", ">f4", (2,)), ("shape", ">f4", (13,)), ("direct", "u1", (13,)),
        ("memberships", ">f4", (13,)), ("auxiliary", ">f4", (4,)), ("distances", ">i4", (3,)),
    ], align=False)
    records = np.empty(spatial, dtype=record_type)
    records["image"] = image
    records["shape"] = shapes
    records["direct"] = np.asarray(direct, dtype=np.uint8)
    records["memberships"] = memberships
    for index, values in enumerate((lesion, black, outlier, lesion_relation)):
        records["auxiliary"][..., index] = values
    records["distances"] = distances
    state = statistics.statistics
    chunks = [
        struct.pack(">3i", *spatial),
        np.asarray([previous_energy], dtype=">f4").tobytes(),
    ]
    for values in (
        state.centroids.class_centroids,
        state.variances.class_standard_deviation,
        state.centroids.lesion_centroids,
        state.variances.lesion_standard_deviation,
        state.channel_weights,
        np.ones(13, dtype=np.float32),
        np.ones(13, dtype=np.float32),
    ):
        chunks.append(np.asarray(values, dtype=">f4").tobytes(order="C"))
    chunks.append(records.tobytes(order="C"))
    return b"".join(chunks)


def _compare_membership_java(
    *,
    java: str,
    javac: str,
    temporary: Path,
    spatial: tuple[int, int, int],
    arrays: dict[str, np.ndarray],
    statistics,
    previous_energy: float,
    source_zero_variance_nan: bool,
) -> dict:
    from compare_p41z_iteration6_membership_java import _parse_java

    temporary.mkdir(parents=True, exist_ok=True)
    input_path, output_path = temporary / "membership-in.bin", temporary / "membership-java.bin"
    input_path.write_bytes(_membership_input_bytes(
        spatial,
        arrays["image"], arrays["shape"], arrays["direct"], arrays["memberships"],
        arrays["lesion"], arrays["black"], arrays["outlier"], arrays["lesion_relation"],
        arrays["distances"], statistics, previous_energy,
    ))
    classes = temporary / "membership-classes"
    classes.mkdir(parents=True, exist_ok=True)
    harness = JAVA_ROOT / "LesionNextIterationMembershipReferenceHarness.java"
    subprocess.run([javac, "-d", str(classes), str(harness)], check=True, capture_output=True, text=True)
    subprocess.run(
        [java, "-Xmx2g", "-cp", str(classes), "LesionNextIterationMembershipReferenceHarness", str(input_path), str(output_path)],
        check=True, capture_output=True, text=True,
    )
    reference = _parse_java(output_path, spatial)
    updated = run_lesion_toads_next_gm_membership_update(
        statistics,
        arrays["image"], arrays["shape"], np.ones(13, np.float32), np.ones(13, np.float32),
        arrays["direct"], arrays["memberships"], arrays["lesion"], arrays["black"],
        arrays["outlier"], arrays["lesion_relation"],
        arrays["distances"][..., 0], arrays["distances"][..., 1], arrays["distances"][..., 2],
        previous_energy=previous_energy, outlier_energy=0.1, lesion_power_exponent=1.2,
        spread=0.5, source_zero_variance_nan=source_zero_variance_nan,
    )
    actual = updated.membership
    pairs = {
        "class_memberships": actual.class_memberships,
        "class_energy": actual.class_energies,
        "lesion": actual.lesion_membership,
        "black": actual.black_hole_membership,
        "outlier": actual.outlier_membership,
        "lesion_relation": actual.lesion_relation,
        "denominator": actual.denominator,
        "lesion_energy": actual.lesion_energy,
        "black_hole_energy": actual.black_hole_energy,
        "running_objective": actual.running_objective,
    }
    exact = {}
    differing = {}
    nan_counts = {"java": int(reference["nan_count"]), "python": 0}
    for name, values in pairs.items():
        values = np.asarray(values, dtype=np.float32)
        expected = np.asarray(reference[name], dtype=np.float32)
        nan_counts["python"] += int(np.isnan(values).sum())
        difference = values.view(np.uint32) != expected.view(np.uint32)
        exact[name] = not bool(np.any(difference))
        if np.any(difference):
            differing[name] = int(np.count_nonzero(difference))
    exact["objective"] = int(np.asarray(actual.objective, np.float32).view(np.uint32)) == int(reference["objective"])
    exact["energy_ratio"] = int(np.asarray(updated.energy_ratio.energy_ratio, np.float32).view(np.uint32)) == int(
        np.asarray(reference["energy_ratio"], dtype=np.float32).view(np.uint32)
    )
    return {
        "all_exact": all(exact.values()) and nan_counts["java"] == nan_counts["python"],
        "exact_float32_bits": exact,
        "differing_value_counts": differing,
        "nan_counts": nan_counts,
        "objective_exact": exact["objective"],
        "energy_ratio_exact": exact["energy_ratio"],
        "previous_energy": float(previous_energy),
        "updated_energy": float(actual.objective),
        "energy_ratio": float(updated.energy_ratio.energy_ratio),
        "stage_order": list(updated.stage_order),
    }


def _growth_bytes(growing) -> bytes:
    output = bytearray()
    output.extend(np.asarray(growing.segmentation, dtype=np.int8).ravel(order="C").tobytes())
    output.extend(np.asarray(growing.skeleton, dtype=np.int8).ravel(order="C").tobytes())
    output.extend(np.asarray(growing.available, dtype=np.uint8).ravel(order="C").tobytes())
    output.extend(np.asarray(growing.ordering, dtype=">f4").ravel(order="C").tobytes())
    output.extend(np.asarray(growing.class_volume, dtype=">f4").tobytes())
    output.extend(struct.pack(">i", len(growing.events)))
    for event in growing.events:
        output.extend(_item_bytes(event.popped_item))
        accepted = -1 if event.topology_accepted is None else int(event.topology_accepted)
        output.extend(struct.pack(">bI", accepted, len(event.inserted_neighbors)))
        for item in event.inserted_neighbors:
            output.extend(_item_bytes(item))
    return bytes(output)


def _compare_growth_java(
    *,
    java: str,
    javac: str,
    temporary: Path,
    thinning,
    memberships,
    object_types,
    lesion_class: int,
    expected_manifest_path: Path,
) -> dict:
    manifest = json.loads(expected_manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or manifest.get("stage") != "plus-lesion-registered-loop":
        raise ValueError("expected next checkpoint must be a complete registered-loop manifest")
    expected_outputs = {
        name: _checked_path(record, name=f"iteration-20 {name}")
        for name, record in manifest.get("outputs", {}).items()
        if isinstance(record, dict)
    }
    temporary.mkdir(parents=True, exist_ok=True)
    spatial = tuple(int(value) for value in thinning.run.segmentation.shape)
    types = tuple(str(value).lower() for value in object_types)
    kind_codes = {"lesion": 1, "blackhole": 2}
    input_path, output_path = temporary / "growth-in.bin", temporary / "growth-java.bin"
    input_path.write_bytes(b"".join((
        struct.pack(">6i", *spatial, len(types), lesion_class, 0),
        bytes(kind_codes.get(name, 0) for name in types),
        np.asarray(thinning.run.segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(thinning.run.skeleton, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(memberships.class_memberships, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(thinning.run.ordering, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(thinning.run.class_volume, dtype=">f4").tobytes(),
        np.asarray(memberships.lesion_membership, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(memberships.black_hole_membership, dtype=">f4").ravel(order="C").tobytes(),
    )))
    classes = temporary / "classes"
    classes.mkdir(parents=True, exist_ok=True)
    upstream = JAVA_ROOT / "upstream/edu/jhmi/rad/medic"
    harness = JAVA_ROOT / "LesionGrowingPassDynamicReferenceHarness.java"
    sources = [
        upstream / "structures/BinaryHeap4D.java",
        upstream / "structures/CriticalPointLUT.java",
        upstream / "utilities/Numerics.java",
        harness,
    ]
    subprocess.run([javac, "-d", str(classes), *(str(path) for path in sources)], check=True, capture_output=True, text=True)
    resource = classes / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
    resource.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(ROOT / "src/pymacruise/data/critical186LUT.raw.gz", resource)
    subprocess.run(
        [java, "-Xmx2g", "-cp", str(classes), "LesionGrowingPassDynamicReferenceHarness", str(input_path), str(output_path)],
        check=True, capture_output=True, text=True,
    )
    python = propagate_lesion_toads_growing(
        thinning.run.segmentation, thinning.run.skeleton, memberships.class_memberships,
        thinning.run.ordering, thinning.run.class_volume, object_types,
        memberships.lesion_membership, memberships.black_hole_membership,
        lesion_class=lesion_class, use_lesion_max=False, retain_event_states=False,
    )
    expected = _growth_bytes(python)
    java_bytes = output_path.read_bytes()
    if java_bytes != expected:
        return {
            "all_exact": False,
            "python_bytes": len(expected),
            "java_bytes": len(java_bytes),
            "python_sha256": hashlib.sha256(expected).hexdigest(),
            "java_sha256": hashlib.sha256(java_bytes).hexdigest(),
        }

    voxel_count = int(np.prod(spatial))
    cursor = 0
    arrays_match = {}
    for name, count, dtype in (
        ("segmentation", voxel_count, np.int8),
        ("skeleton", voxel_count, np.int8),
        ("availability", voxel_count * len(types), np.uint8),
        ("ordering", voxel_count, ">f4"),
        ("class_volume", len(types), ">f4"),
    ):
        width = np.dtype(dtype).itemsize
        actual = np.frombuffer(java_bytes[cursor:cursor + count * width], dtype=dtype)
        cursor += count * width
        if name == "class_volume":
            summary = json.loads(expected_outputs["summary"].read_text(encoding="utf-8"))
            reference = np.asarray(summary["class_volume"], dtype=np.float32).ravel(order="C")
        else:
            reference = np.asarray(nib.load(str(expected_outputs[name])).dataobj)
            reference = reference.ravel(order="C")
        if np.issubdtype(actual.dtype, np.floating):
            actual_bits = actual.astype(np.float32).view(np.uint32)
            reference_bits = np.asarray(reference, dtype=np.float32).view(np.uint32)
            arrays_match[name] = np.array_equal(actual_bits, reference_bits)
        else:
            arrays_match[name] = np.array_equal(actual, np.asarray(reference, dtype=actual.dtype))
    return {
        "all_exact": True,
        "checkpoint_outputs_exact": all(arrays_match.values()),
        "checkpoint_arrays_exact": arrays_match,
        "event_count": len(python.events),
        "java_bytes": len(java_bytes),
        "java_sha256": hashlib.sha256(java_bytes).hexdigest(),
        "python_sha256": hashlib.sha256(expected).hexdigest(),
        "source_files": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in (*sources, ROOT / "src/pymacruise/data/critical186LUT.raw.gz")
        },
    }


def compare(
    manifest_path: Path,
    output_dir: Path,
    *,
    java: str = "java",
    javac: str = "javac",
    expected_next_manifest: Path | None = None,
) -> dict:
    manifest_file = manifest_path.expanduser().resolve()
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete" or manifest.get("stage") != "plus-lesion-registered-loop":
        raise ValueError("input must be a completed registered-loop checkpoint")
    if manifest.get("spatial_convention") != "RAS":
        raise ValueError("registered-loop checkpoint must be RAS")
    outputs = {
        name: _checked_path(record, name=f"registered {name}")
        for name, record in manifest.get("outputs", {}).items()
        if isinstance(record, dict)
    }
    recrop_record = manifest["inputs"]["registered_recrop_manifest"]
    recrop_path = _checked_path(recrop_record, name="registered recrop manifest")
    recrop = json.loads(recrop_path.read_text(encoding="utf-8"))
    modality_path = _checked_path(recrop["outputs"]["modalities"], name="registered modalities")
    modalities = nib.load(str(modality_path))
    if len(modalities.shape) != 4 or modalities.shape[-1] != 2 or nib.aff2axcodes(modalities.affine) != ("R", "A", "S"):
        raise ValueError("registered modalities must be two-channel RAS NIfTI")

    branch, registration, distance_state = hydrate_plus_lesion_registered_branch(manifest_file)
    atlas_path = _checked_path(manifest["inputs"]["atlas_descriptor"], name="atlas descriptor")
    atlas_classes = load_lesion_atlas_priors(atlas_path, ("T1", "FLAIR"))
    atlas_shapes = load_lesion_atlas_shape_priors(atlas_path)
    atlas_topology = load_lesion_atlas_topology(atlas_path)
    shape_priors = tuple(atlas_shapes.shape_maps[item.name] for item in atlas_classes.structures)
    image_resolution = tuple(float(value) for value in np.linalg.norm(modalities.affine[:3, :3], axis=0))
    template_resolution = tuple(
        float(value) for value in np.linalg.norm(atlas_topology.image.affine[:3, :3], axis=0)
    )
    dimensions = tuple(int(value) for value in modalities.shape[:3])
    shapes = transform_registered_lesion_shape_priors(
        registration,
        shape_priors,
        image_dimensions_xyz=dimensions,
        image_resolutions_mm=image_resolution,
        template_resolutions_mm=template_resolution,
    )
    settings = manifest["settings"]
    thinning = run_lesion_toads_next_iteration_thinning(
        branch,
        shapes,
        object_types=STANDARD13_TYPES,
        lesion_class=10,
        volume_ratio=float(settings["first_thinning_volume_ratio"]),
        limit=float(settings["thinning_limit"]),
        use_lesion_max=False,
    )

    summary = json.loads(outputs["summary"].read_text(encoding="utf-8"))
    packed = np.asarray(nib.load(str(outputs["memberships"])).dataobj, dtype=np.float32)
    statistics_input = output_dir / "statistics-input"
    statistics_input.mkdir(parents=True, exist_ok=True)
    initial_case = json.loads(
        (ROOT / "tests/data/synthetic/p41z_iteration6_statistics_java/case.json").read_text(
            encoding="utf-8"
        )
    )
    statistics_case = {
        "case": "P41z iteration-20 post-thinning statistics source comparison",
        "modalities": initial_case["modalities"],
        "object_types": initial_case["object_types"],
        "initial_class_centroids": summary["refreshed_centroids"]["class_centroids"],
        "initial_lesion_centroids": summary["refreshed_centroids"]["lesion_centroids"],
        "initial_black_hole_centroids": summary["refreshed_centroids"]["black_hole_centroids"],
        "tissue_power_exponents": initial_case["tissue_power_exponents"],
        "lesion_power_exponents": initial_case["lesion_power_exponents"],
        "norm_type": initial_case["norm_type"],
        "outlier_distance": initial_case["outlier_distance"],
        "optimized_csf_all_modalities": initial_case["optimized_csf_all_modalities"],
        "input_files": {
            "images": "images_t1_flair.nii.gz",
            "memberships": "class_memberships_13ch.nii.gz",
            "segmentation": "iteration20_thinned_segmentation.nii.gz",
            "lesion_membership": "lesion_membership.nii.gz",
            "black_hole_membership": "black_hole_membership.nii.gz",
            "outlier_membership": "outlier_membership.nii.gz",
        },
    }
    stats_arrays = {
        "images_t1_flair.nii.gz": np.asarray(modalities.dataobj, dtype=np.float32),
        "class_memberships_13ch.nii.gz": packed[..., :13],
        "iteration20_thinned_segmentation.nii.gz": np.asarray(thinning.run.segmentation, dtype=np.int16),
        "lesion_membership.nii.gz": packed[..., 13],
        "black_hole_membership.nii.gz": packed[..., 14],
        "outlier_membership.nii.gz": packed[..., 15],
    }
    for name, values in stats_arrays.items():
        nib.save(nib.Nifti1Image(values, modalities.affine), str(statistics_input / name))
    (statistics_input / "case.json").write_text(
        json.dumps(statistics_case, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    statistics_report = compare_statistics(statistics_input, java=java, javac=javac)
    summary = json.loads(outputs["summary"].read_text(encoding="utf-8"))
    image_values = np.asarray(modalities.dataobj, dtype=np.float32)
    prior_membership = branch.membership_result.membership
    relation_state = branch.relations
    direct_relations = relation_state.relation_membership[relation_state.relation_map]
    previous_distance_maps = LesionInitialDistanceMaps(
        distance_state[..., 0].astype(np.bool_),
        distance_state[..., 1].astype(np.int8),
        distance_state[..., 2].astype(np.bool_),
        distance_state[..., 3].astype(np.int8),
        distance_state[..., 4].astype(np.bool_),
        distance_state[..., 5].astype(np.int8),
    )
    refreshed_distance_state = refresh_lesion_toads_distance_maps(
        thinning.run.segmentation,
        tuple(item.name for item in atlas_classes.structures),
        maximum_ventricle_distance=2,
        maximum_gm_distance=3,
        maximum_brainstem_distance=25,
        previous_maps=previous_distance_maps,
    )
    refreshed = refreshed_distance_state.maps
    distances = np.stack(
        (refreshed.ventricle, refreshed.gray_matter, refreshed.brainstem), axis=-1
    ).astype(np.int32, copy=False)
    prior_stats = refresh_statistics_after_next_iteration_thinning(
        thinning,
        image_values,
        prior_membership.class_memberships,
        prior_membership.lesion_membership,
        prior_membership.black_hole_membership,
        prior_membership.outlier_membership,
        summary["refreshed_centroids"]["class_centroids"],
        summary["refreshed_centroids"]["lesion_centroids"],
        summary["refreshed_centroids"]["black_hole_centroids"],
        STANDARD13_TYPES,
        centroid_options={
            "mode": "none",
            "modalities": ("T1_MPRAGE", "FLAIR"),
            "optimized_csf_all_modalities": True,
            "tissue_power_exponents": np.ones((2, 13), dtype=np.float32),
            "lesion_power_exponents": np.asarray([0.4, 1.2], dtype=np.float32),
        },
        channel_weight_options={"outlier_distance": 0.1, "norm_type": "Identity"},
        use_lesion_weight=False,
        flair_channel=1,
    )
    membership_arrays = {
        "image": image_values,
        "shape": shapes,
        "direct": direct_relations,
        "memberships": prior_membership.class_memberships,
        "lesion": prior_membership.lesion_membership,
        "black": prior_membership.black_hole_membership,
        "outlier": prior_membership.outlier_membership,
        "lesion_relation": prior_membership.lesion_relation,
        "distances": distances,
    }
    membership_report = _compare_membership_java(
        java=java,
        javac=javac,
        temporary=output_dir / "membership-java-work",
        spatial=dimensions,
        arrays=membership_arrays,
        statistics=prior_stats,
        previous_energy=float(prior_membership.objective),
        source_zero_variance_nan=bool(settings["source_zero_variance_nan"]),
    )

    growth_report = None
    if expected_next_manifest is not None:
        next_membership = run_lesion_toads_next_gm_membership_update(
            prior_stats,
            membership_arrays["image"], membership_arrays["shape"],
            np.ones(13, dtype=np.float32), np.ones(13, dtype=np.float32),
            membership_arrays["direct"], membership_arrays["memberships"],
            membership_arrays["lesion"], membership_arrays["black"],
            membership_arrays["outlier"], membership_arrays["lesion_relation"],
            distances[..., 0], distances[..., 1], distances[..., 2],
            previous_energy=float(prior_membership.objective), outlier_energy=0.1,
            lesion_power_exponent=1.2, spread=0.5,
            source_zero_variance_nan=bool(settings["source_zero_variance_nan"]),
        )
        growth_report = _compare_growth_java(
            java=java,
            javac=javac,
            temporary=output_dir / "growth-java-work",
            thinning=thinning,
            memberships=next_membership.membership,
            object_types=STANDARD13_TYPES,
            lesion_class=10,
            expected_manifest_path=expected_next_manifest.expanduser().resolve(),
        )

    with tempfile.TemporaryDirectory(prefix="p41z-registered-thinning-java-") as temporary:
        work = Path(temporary)
        classes_dir = work / "classes"
        classes_dir.mkdir()
        sources = [
            JAVA_ROOT / "stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java",
            UPSTREAM / "utilities/Numerics.java",
            UPSTREAM / "structures/BinaryHeap4D.java",
            UPSTREAM / "structures/CriticalPointLUT.java",
            UPSTREAM / "structures/RotationMatrix.java",
            UPSTREAM / "methods/ParametricTransform.java",
            UPSTREAM / "libraries/ImageFunctions.java",
            JAVA_ROOT / "Lesion13ClassNextThinningReferenceHarness.java",
            JAVA_ROOT / "P41zRegisteredShapeVolumeReferenceHarness.java",
        ]
        subprocess.run([javac, "-d", str(classes_dir), *(str(path) for path in sources)], check=True, capture_output=True, text=True)
        lut_resource = classes_dir / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
        lut_resource.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / "src/pymacruise/data/critical186LUT.raw.gz", lut_resource)

        template_dims = tuple(int(value) for value in shape_priors[0].shape)
        shape_input = work / "shapes-in.bin"
        shape_output = work / "shapes-java.bin"
        shape_input.write_bytes(b"".join((
            struct.pack(">7i", *dimensions, *template_dims, len(shape_priors)),
            struct.pack(">6f", *image_resolution, *template_resolution),
            np.asarray(registration.transform, dtype=">f4").tobytes(),
            b"".join(
                struct.pack(">6i", *box.lower_xyz, *box.upper_xyz)
                for box in registration.shape_bounds
            ),
            b"".join(np.asarray(prior, dtype=">f4").ravel(order="C").tobytes() for prior in shape_priors),
        )))
        subprocess.run(
            [java, "-Xmx2g", "-cp", str(classes_dir), "P41zRegisteredShapeVolumeReferenceHarness", str(shape_input), str(shape_output)],
            check=True, capture_output=True, text=True,
        )
        python_shapes = np.asarray(shapes, dtype=">f4").ravel(order="C").tobytes()
        shape_exact = shape_output.read_bytes() == python_shapes

        segmentation = np.asarray(nib.load(str(outputs["segmentation"])).dataobj, dtype=np.int16)
        skeleton = np.asarray(nib.load(str(outputs["skeleton"])).dataobj, dtype=np.int16)
        availability = np.asarray(nib.load(str(outputs["availability"])).dataobj, dtype=bool)
        ordering = np.asarray(nib.load(str(outputs["ordering"])).dataobj, dtype=np.float32)
        thinning_input = work / "thinning-in.bin"
        thinning_output = work / "thinning-java.bin"
        type_bytes = bytes(TYPE_CODES.get(str(name).lower(), 0) for name in STANDARD13_TYPES)
        thinning_input.write_bytes(b"".join((
            struct.pack(">5i2f", *dimensions, 13, 10,
                        float(settings["first_thinning_volume_ratio"]), float(settings["thinning_limit"])),
            type_bytes,
            np.asarray(segmentation, dtype=np.int8).ravel(order="C").tobytes(),
            np.asarray(skeleton, dtype=np.int8).ravel(order="C").tobytes(),
            np.asarray(packed[..., :13], dtype=">f4").ravel(order="C").tobytes(),
            np.asarray(packed[..., 13], dtype=">f4").ravel(order="C").tobytes(),
            np.asarray(packed[..., 14], dtype=">f4").ravel(order="C").tobytes(),
            python_shapes,
            np.asarray(availability, dtype=np.uint8).ravel(order="C").tobytes(),
            np.asarray(ordering, dtype=">f4").ravel(order="C").tobytes(),
        )))
        subprocess.run(
            [java, "-Xmx2g", "-cp", str(classes_dir), "Lesion13ClassNextThinningReferenceHarness", str(thinning_input), str(thinning_output)],
            check=True, capture_output=True, text=True,
        )
        remaining_queue_count = len(thinning.run.queue)
        expected = _python_thinning_bytes(thinning)
        java_bytes = thinning_output.read_bytes()

    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "synthetic_only": True,
        "comparison_scope": "selected Java shape sampling and next-thinning methods on a hash-verified P41z iteration-19 registered state",
        "workflow_iteration": int(json.loads(outputs["summary"].read_text(encoding="utf-8"))["completed_iteration"]),
        "shape_grid_xyz": list(dimensions),
        "shape_exact": shape_exact,
        "shape_values": int(shapes.size),
        "shape_maximum_absolute_error": (
            float(np.max(np.abs(
                np.frombuffer(shape_output.read_bytes(), dtype=">f4").astype(np.float32)
                - np.asarray(shapes, dtype=np.float32).ravel(order="C")
            ))) if not shape_exact else 0.0
        ),
        "thinning_exact_serialization": java_bytes == expected,
        "post_thinning_statistics_exact_float32_bits": statistics_report["all_exact"],
        "post_thinning_statistics_nan_counts": statistics_report["nan_counts"],
        "post_thinning_statistics_exact_fields": statistics_report["exact_float32_bits"],
        "next_membership_exact_float32_bits": membership_report["all_exact"],
        "next_membership_exact_fields": membership_report["exact_float32_bits"],
        "next_membership_nan_counts": membership_report["nan_counts"],
        "next_membership_updated_energy": membership_report["updated_energy"],
        "next_membership_energy_ratio": membership_report["energy_ratio"],
        "post_registration_growth": growth_report,
        "distance_refresh_stage_order": list(refreshed_distance_state.stage_order),
        "refreshed_distance_hashes": {
            "ventricle_mask": hashlib.sha256(refreshed.ventricle_mask.tobytes(order="C")).hexdigest(),
            "ventricle": hashlib.sha256(refreshed.ventricle.tobytes(order="C")).hexdigest(),
            "gray_matter_mask": hashlib.sha256(refreshed.gray_matter_mask.tobytes(order="C")).hexdigest(),
            "gray_matter": hashlib.sha256(refreshed.gray_matter.tobytes(order="C")).hexdigest(),
            "brainstem_mask": hashlib.sha256(refreshed.brainstem_mask.tobytes(order="C")).hexdigest(),
            "brainstem": hashlib.sha256(refreshed.brainstem.tobytes(order="C")).hexdigest(),
        },
        "initial_queue_items": len(thinning.initialization.inserted_items),
        "queue_events": len(thinning.run.events),
        "topology_candidate_events": sum(event.candidate_label is not None for event in thinning.run.events),
        "remaining_queue_items": remaining_queue_count,
        "python_serialized_sha256": hashlib.sha256(expected).hexdigest(),
        "java_serialized_sha256": hashlib.sha256(java_bytes).hexdigest(),
        "java_serialized_bytes": len(java_bytes),
        "source_files": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in sources
        },
        "checkpoint_manifest": str(manifest_file),
        "checkpoint_sha256": _sha256(manifest_file),
        "claim_limit": "fixed transform and support bounds; not full registered-loop or whole-workflow equivalence",
    }
    report_path = output_dir / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    if (
        not shape_exact or java_bytes != expected or not statistics_report["all_exact"]
        or not membership_report["all_exact"]
        or (growth_report is not None and not growth_report.get("all_exact"))
        or (growth_report is not None and not growth_report.get("checkpoint_outputs_exact"))
    ):
        raise SystemExit(1)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint_manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--java", default="java")
    parser.add_argument("--javac", default="javac")
    parser.add_argument("--expected-next-checkpoint", type=Path)
    args = parser.parse_args()
    compare(
        args.checkpoint_manifest,
        args.output_dir,
        java=args.java,
        javac=args.javac,
        expected_next_manifest=args.expected_next_checkpoint,
    )


if __name__ == "__main__":
    main()
