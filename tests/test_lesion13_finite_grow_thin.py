from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
from types import SimpleNamespace

import nibabel as nib
import numpy as np

from pymacruise.lesion_growth import (
    propagate_lesion_toads_growing,
)
from pymacruise.lesion_next_iteration_thinning import run_lesion_toads_next_iteration_thinning
from pymacruise.lesion_next_iteration_statistics import refresh_statistics_after_next_iteration_thinning
from pymacruise.lesion_next_iteration_membership import run_lesion_toads_next_gm_membership_update
from pymacruise.lesion_post_registration import LesionToadsPostRegistrationBranch
from test_lesion_nan_post_registration_growth_reference import (
    HARNESS as GROW_HARNESS,
    HEAP,
    LUT_GZIP,
    LUT_JAVA,
    NUMERICS,
    _expected_bytes as expected_growth_bytes,
    _input_bytes as growth_input_bytes,
)


ROOT = Path(__file__).parent / "data/synthetic/lesion13_finite_grow_thin"
SPEC = json.loads((ROOT / "fixture.json").read_text(encoding="utf-8"))
REFERENCE = json.loads((ROOT / "java_reference.json").read_text(encoding="utf-8"))
THIN_HARNESS = Path(__file__).parent / "java_reference/Lesion13ClassNextThinningReferenceHarness.java"
STATISTICS_HARNESS = Path(__file__).parent / "java_reference/LesionPostThinningStatisticsReferenceHarness.java"
POWER_TABLE = Path(__file__).parent / "java_reference/upstream/edu/jhmi/rad/medic/utilities/PowerTable.java"
MEMBERSHIP_HARNESS = Path(__file__).parent / "java_reference/LesionPostThinningStandard13MembershipReferenceHarness.java"
RATIO_HARNESS = Path(__file__).parent / "java_reference/LesionDynamicEnergyRatioReferenceHarness.java"
TYPES = tuple(SPEC["object_types"])
CLASSES = int(SPEC["class_count"])
LESION = int(SPEC["lesion_class"])


def _item_bytes(item) -> bytes:
    return struct.pack(">fiiii", float(item.value), item.x, item.y, item.z, item.label)


def _thinning_input_bytes(seg, skel, memberships, lesion, black, shapes, available, order) -> bytes:
    nx, ny, nz = (int(value) for value in seg.shape)
    kinds = bytes(1 if value == "lesion" else 2 if value == "blackhole" else 0 for value in TYPES)
    fields = [
        struct.pack(">5i2f", nx, ny, nz, CLASSES, LESION, 0.4, 1000.0),
        kinds,
        np.asarray(seg, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(skel, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(memberships, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(lesion, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(black, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(shapes, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(available, dtype=np.uint8).ravel(order="C").tobytes(),
        np.asarray(order, dtype=">f4").ravel(order="C").tobytes(),
    ]
    return b"".join(fields)


def _expected_thinning_bytes(result) -> bytes:
    fields = [
        np.asarray(result.initialization.class_volume, dtype=">f4").tobytes(),
        np.asarray(result.minimum_class_volume, dtype=">f4").tobytes(),
        struct.pack(">i", len(result.initialization.inserted_items)),
    ]
    fields.extend(_item_bytes(item) for item in result.initialization.inserted_items)
    fields.append(struct.pack(">i", len(result.run.events)))
    for event in result.run.events:
        fields.extend((_item_bytes(event.popped_item), struct.pack(">bi", int(event.processed), len(event.inserted_neighbors))))
        fields.extend(_item_bytes(item) for item in event.inserted_neighbors)
    fields.extend((
        np.asarray(result.run.segmentation, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(result.run.skeleton, dtype=np.int8).ravel(order="C").tobytes(),
        np.asarray(result.run.available, dtype=np.uint8).ravel(order="C").tobytes(),
        np.asarray(result.run.ordering, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(result.run.class_volume, dtype=">f4").tobytes(),
    ))
    remaining = []
    while result.run.queue:
        remaining.append(result.run.queue.pop_first())
    fields.append(struct.pack(">i", len(remaining)))
    fields.extend(_item_bytes(item) for item in remaining)
    return b"".join(fields)


def _statistics_input_bytes(images, memberships, segmentation, lesion, black, centroids, lesion_centroid, black_centroid) -> bytes:
    fields = [
        np.asarray(centroids, dtype=">f4").ravel(order="C").tobytes(),
        np.asarray(lesion_centroid, dtype=">f4").tobytes(),
        np.asarray(black_centroid, dtype=">f4").tobytes(),
    ]
    for x in range(segmentation.shape[0]):
        for y in range(segmentation.shape[1]):
            for z in range(segmentation.shape[2]):
                fields.extend((
                    np.asarray(images[x, y, z], dtype=">f4").tobytes(),
                    np.asarray(memberships[x, y, z], dtype=">f4").tobytes(),
                    struct.pack(">i", int(segmentation[x, y, z])),
                    np.asarray(lesion[x, y, z], dtype=">f4").tobytes(),
                    np.asarray(black[x, y, z], dtype=">f4").tobytes(),
                ))
    return b"".join(fields)


def _float_bits(values):
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).tolist()


def _flat_float_bits(values):
    return np.asarray(values, dtype=np.float32).view(np.uint32).astype(np.uint64).ravel().tolist()


def test_finite_13class_state_growth_then_next_thinning_matches_java(tmp_path: Path) -> None:
    assert REFERENCE["lesion_toads_source_sha256"] == "545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99"
    assert REFERENCE["wrapper_source_sha256"] == "556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb"
    assert hashlib.sha256(GROW_HARNESS.read_bytes()).hexdigest() == REFERENCE["grow_harness_sha256"]
    assert hashlib.sha256(THIN_HARNESS.read_bytes()).hexdigest() == REFERENCE["thin_harness_sha256"]
    assert hashlib.sha256(HEAP.read_bytes()).hexdigest() == REFERENCE["binary_heap_sha256"]
    assert hashlib.sha256(NUMERICS.read_bytes()).hexdigest() == REFERENCE["numerics_sha256"]
    assert hashlib.sha256(STATISTICS_HARNESS.read_bytes()).hexdigest() == REFERENCE["statistics_harness_sha256"]
    assert hashlib.sha256(MEMBERSHIP_HARNESS.read_bytes()).hexdigest() == REFERENCE["membership_harness_sha256"]
    assert hashlib.sha256(RATIO_HARNESS.read_bytes()).hexdigest() == REFERENCE["energy_ratio_harness_sha256"]
    assert hashlib.sha256(POWER_TABLE.read_bytes()).hexdigest() == REFERENCE["power_table_sha256"]
    image_path = ROOT / "finite_13class_grow_thin_ras.nii.gz"
    assert hashlib.sha256(image_path.read_bytes()).hexdigest() == SPEC["sha256_nii_gz"]
    image = nib.load(str(image_path))
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    packed = np.asarray(image.dataobj, dtype=np.float32)
    assert hashlib.sha256(packed.astype("<f4", copy=False).tobytes()).hexdigest() == SPEC["sha256_uncompressed_float32le"]
    seg = packed[..., 0].astype(np.int16)
    skel = packed[..., 1].astype(np.int16)
    memberships = packed[..., 2:15]
    lesion, black = packed[..., 15], packed[..., 16]
    shapes, order = packed[..., 17:30], packed[..., 30]
    images = packed[..., 31:33]
    direct = packed[..., 33:46] > 0.5
    lesion_relation = packed[..., 46]
    ventricle = packed[..., 47].astype(np.int16)
    gray_matter = packed[..., 48].astype(np.int16)
    brainstem = packed[..., 49].astype(np.int16)
    assert memberships.shape[-1] == len(TYPES) == CLASSES == 13
    assert np.isfinite(packed).all()
    assert tuple(packed.shape[:3]) == tuple(SPEC["shape_xyz"])
    support = np.bincount(seg[1:-1, 1:-1, 1:-1].ravel(), minlength=CLASSES)
    assert np.all(support > 0)
    assert int(np.count_nonzero((seg == LESION) & (memberships[..., LESION] < lesion))) == REFERENCE["lesion_auxiliary_winner_count"]
    black_only_winners = (seg == LESION) & (memberships[..., LESION] >= lesion) & (memberships[..., LESION] < black)
    assert int(np.count_nonzero(black_only_winners)) == REFERENCE["black_hole_winner_without_lesion_eligibility_count"]
    center = tuple(SPEC["growth_center_xyz"])
    assert skel[center] == -1

    class_volume = np.zeros(CLASSES, dtype=np.float32)
    interior_seg = seg[1:-1, 1:-1, 1:-1]
    interior_memberships = memberships[1:-1, 1:-1, 1:-1]
    for k in range(CLASSES):
        values = interior_memberships[..., k][interior_seg == k]
        if k == LESION:
            values = values + lesion[1:-1, 1:-1, 1:-1][interior_seg == k]
            values = values + black[1:-1, 1:-1, 1:-1][interior_seg == k]
        class_volume[k] = np.sum(values, dtype=np.float32)
    grown = propagate_lesion_toads_growing(
        seg, skel, memberships, order, class_volume, TYPES, lesion, black,
        lesion_class=LESION,
    )
    assert len(grown.events) == 3
    assert sum(event.topology_accepted is True for event in grown.events) == 1
    assert len(grown.events) == REFERENCE["growth_event_count"]
    assert sum(event.topology_accepted is True for event in grown.events) == REFERENCE["growth_accepted_count"]
    assert grown.segmentation[center] == SPEC["growth_center_class"]
    assert grown.remaining_queue == ()

    classes_dir = tmp_path / "classes"
    classes_dir.mkdir()
    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac and java, "a JDK is required for the pinned Java method references"
    subprocess.run(
        [javac, "-d", str(classes_dir), str(HEAP), str(LUT_JAVA), str(NUMERICS),
         str(GROW_HARNESS), str(THIN_HARNESS)],
        check=True, capture_output=True, text=True,
    )
    resource = classes_dir / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
    resource.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(LUT_GZIP, resource)
    growth_input, growth_output = tmp_path / "grow-in.bin", tmp_path / "grow-out.bin"
    growth_input.write_bytes(growth_input_bytes(seg, skel, memberships, order, class_volume, lesion, black))
    subprocess.run([java, "-cp", str(classes_dir), "LesionGrowingPassDynamicReferenceHarness", str(growth_input), str(growth_output)], check=True, capture_output=True, text=True)
    assert growth_output.read_bytes() == expected_growth_bytes(grown)

    continued_branch = LesionToadsPostRegistrationBranch(
        membership_result=SimpleNamespace(
            membership=SimpleNamespace(
                class_memberships=memberships,
                lesion_membership=lesion,
                black_hole_membership=black,
            )
        ),
        registration_handoff=SimpleNamespace(),
        energy_ratio=np.float32(0.5),
        convergence_stop=False,
        iteration_limit_stop=False,
        stop=False,
        stop_reason=None,
        run_growing=True,
        growing=grown,
        segmentation=grown.segmentation,
        final_classification=None,
        relations=SimpleNamespace(),
        stage_order=("propagateGrowing", "computeRelations"),
    )
    thinned = run_lesion_toads_next_iteration_thinning(
        continued_branch,
        shapes,
        object_types=TYPES,
        lesion_class=LESION,
        volume_ratio=0.4,
        limit=1000.0,
    )
    assert len(thinned.initialization.inserted_items) == REFERENCE["thinning_initial_queue_count"]
    assert len(thinned.run.events) == REFERENCE["thinning_event_count"]
    assert sum(event.processed for event in thinned.run.events) == REFERENCE["thinning_processed_count"]
    thin_input, thin_output = tmp_path / "thin-in.bin", tmp_path / "thin-out.bin"
    thin_input.write_bytes(_thinning_input_bytes(
        continued_branch.segmentation, grown.skeleton, memberships, lesion, black, shapes,
        grown.available, grown.ordering,
    ))
    subprocess.run([java, "-cp", str(classes_dir), "Lesion13ClassNextThinningReferenceHarness", str(thin_input), str(thin_output)], check=True, capture_output=True, text=True)
    assert thin_output.read_bytes() == _expected_thinning_bytes(thinned)
    assert len(thinned.run.events) > 0

    centroid_options = {
        "mode": "none",
        "modalities": tuple(SPEC["modalities"]),
        "lesion_power_exponents": (0.4, 1.2),
    }
    initial_centroids = np.asarray(SPEC["initial_class_centroids"], dtype=np.float32)
    initial_lesion = np.asarray(SPEC["initial_lesion_centroids"], dtype=np.float32)
    initial_black = np.asarray(SPEC["initial_black_hole_centroids"], dtype=np.float32)
    outlier = np.full(seg.shape, SPEC["outlier_membership_constant"], dtype=np.float32)
    refreshed = refresh_statistics_after_next_iteration_thinning(
        thinned,
        images,
        memberships,
        lesion,
        black,
        outlier,
        initial_centroids,
        initial_lesion,
        initial_black,
        TYPES,
        centroid_options=centroid_options,
        channel_weight_options={"outlier_distance": 1.0},
        flair_channel=1,
    )
    stats_input, stats_output = tmp_path / "stats-in.bin", tmp_path / "stats-out.json"
    stats_input.write_bytes(_statistics_input_bytes(
        images,
        memberships,
        thinned.run.segmentation,
        lesion,
        black,
        initial_centroids,
        initial_lesion,
        initial_black,
    ))
    subprocess.run(
        [javac, "-d", str(classes_dir), str(POWER_TABLE), str(STATISTICS_HARNESS)],
        check=True, capture_output=True, text=True,
    )
    java_stats = json.loads(subprocess.run(
        [java, "-cp", str(classes_dir), "LesionPostThinningStatisticsReferenceHarness", str(stats_input)],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert java_stats["source_sha256"] == REFERENCE["lesion_toads_source_sha256"]
    assert java_stats["wrapper_sha256"] == REFERENCE["wrapper_source_sha256"]
    assert java_stats["class"] == _float_bits(refreshed.statistics.centroids.class_centroids)
    assert java_stats["aux_centroids"][0] == _float_bits(refreshed.statistics.centroids.lesion_centroids)
    assert java_stats["aux_centroids"][1] == _float_bits(refreshed.statistics.centroids.black_hole_centroids)
    assert java_stats["class_std"] == _float_bits(refreshed.statistics.variances.class_standard_deviation.reshape(-1))
    assert java_stats["lesion_std"] == _float_bits(refreshed.statistics.variances.lesion_standard_deviation)
    assert java_stats["blackhole_std"] == _float_bits(refreshed.statistics.variances.black_hole_standard_deviation)
    assert java_stats["weights"] == _float_bits(refreshed.statistics.channel_weights)
    assert refreshed.stage_order[-3:] == ("computeCentroids", "computeVariances", "computeChannelWeights")
    assert np.any(refreshed.statistics.centroids.lesion_centroids != 0.0)
    assert np.all(refreshed.statistics.variances.lesion_standard_deviation > 0.0)
    # Centroid routing uses max(mem, lesion, black-hole), so the black-only
    # winner contributes to the auxiliary black-hole centroid. Variance has a
    # stricter source predicate (m < lesion membership for a LESION descriptor),
    # so that same voxel remains in class variance rather than black-hole std.
    assert np.all(refreshed.statistics.centroids.black_hole_centroids > 0.0)
    assert np.array_equal(refreshed.statistics.variances.black_hole_standard_deviation, np.zeros(2, np.float32))

    update = run_lesion_toads_next_gm_membership_update(
        refreshed,
        images,
        shapes,
        SPEC["shape_priors"],
        SPEC["smoothing_priors"],
        direct,
        memberships,
        lesion,
        black,
        outlier,
        lesion_relation,
        ventricle,
        gray_matter,
        brainstem,
        previous_energy=SPEC["previous_energy"],
        outlier_energy=SPEC["outlier_energy"],
        lesion_power_exponent=SPEC["lesion_power_exponent"],
    )
    assert update.stage_order[-2:] == ("computeMemberships:GMflag=true", "energyRatio")

    membership_input = bytearray()
    for values in (
        refreshed.statistics.centroids.class_centroids,
        refreshed.statistics.variances.class_standard_deviation,
        refreshed.statistics.centroids.lesion_centroids,
        refreshed.statistics.variances.lesion_standard_deviation,
        refreshed.statistics.channel_weights,
        np.asarray(SPEC["shape_priors"], np.float32),
        np.asarray(SPEC["smoothing_priors"], np.float32),
    ):
        membership_input.extend(np.asarray(values, dtype=">f4").ravel(order="C").tobytes())
    for x in range(seg.shape[0]):
        for y in range(seg.shape[1]):
            for z in range(seg.shape[2]):
                membership_input.extend(np.asarray(images[x, y, z], dtype=">f4").tobytes())
                membership_input.extend(np.asarray(shapes[x, y, z], dtype=">f4").tobytes())
                membership_input.extend(np.asarray(direct[x, y, z], dtype=">i4").tobytes())
                membership_input.extend(np.asarray(memberships[x, y, z], dtype=">f4").tobytes())
                membership_input.extend(np.asarray((lesion[x, y, z], black[x, y, z], outlier[x, y, z], lesion_relation[x, y, z]), dtype=">f4").tobytes())
                membership_input.extend(struct.pack(">3i", int(ventricle[x, y, z]), int(gray_matter[x, y, z]), int(brainstem[x, y, z])))
    membership_path = tmp_path / "membership-in.bin"
    membership_path.write_bytes(bytes(membership_input))
    subprocess.run([javac, "-d", str(classes_dir), str(POWER_TABLE), str(MEMBERSHIP_HARNESS)], check=True, capture_output=True, text=True)
    java_membership = json.loads(subprocess.run([java, "-cp", str(classes_dir), "LesionPostThinningStandard13MembershipReferenceHarness", str(membership_path)], check=True, capture_output=True, text=True).stdout)
    actual = update.membership
    for key, values in {
        "class_membership_bits": actual.class_memberships,
        "class_energy_bits": actual.class_energies,
        "lesion_membership_bits": actual.lesion_membership,
        "black_hole_membership_bits": actual.black_hole_membership,
        "outlier_membership_bits": actual.outlier_membership,
        "lesion_relation_bits": actual.lesion_relation,
        "denominator_bits": actual.denominator,
        "lesion_energy_bits": actual.lesion_energy,
        "black_hole_energy_bits": actual.black_hole_energy,
        "running_objective_bits": actual.running_objective,
    }.items():
        assert java_membership[key] == _flat_float_bits(values), key
    assert java_membership["objective_bits"] == _float_bits(actual.objective)

    ratio_input = tmp_path / "ratio-in.bin"
    ratio_input.write_bytes(struct.pack(">2f", float(SPEC["previous_energy"]), float(actual.objective)))
    subprocess.run([javac, "-d", str(classes_dir), str(NUMERICS), str(RATIO_HARNESS)], check=True, capture_output=True, text=True)
    java_ratio = json.loads(subprocess.run([java, "-cp", str(classes_dir), "LesionDynamicEnergyRatioReferenceHarness", str(ratio_input)], check=True, capture_output=True, text=True).stdout)
    ratio = update.energy_ratio
    assert java_ratio["previous_energy_bits"] == f"{int(ratio.previous_energy.view(np.uint32)):08x}"
    assert java_ratio["current_energy_bits"] == f"{int(ratio.current_energy.view(np.uint32)):08x}"
    for field, value in {
        "difference_bits": ratio.difference,
        "numerator_bits": ratio.numerator,
        "denominator_bits": ratio.denominator,
        "quotient_bits": ratio.quotient,
        "energy_ratio_bits": ratio.energy_ratio,
    }.items():
        assert java_ratio[field] == f"{int(np.asarray(value, np.float32).view(np.uint32)):08x}"
