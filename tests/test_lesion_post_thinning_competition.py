from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_post_thinning_competition import run_post_thinning_competition_pass
from test_lesion_post_thinning_gm_membership import SPEC, _compose


ROOT = Path(__file__).parent
JAVA_ROOT = ROOT / "java_reference"
JAVA_UPSTREAM = JAVA_ROOT / "upstream/edu/jhmi/rad/medic/structures"
HARNESS = JAVA_ROOT / "LesionToadsComposedCompetitionReferenceHarness.java"
LUT = ROOT.parent / "src/pymacruise/data/critical186LUT.raw.gz"
INITIAL_REFERENCE = json.loads(
    (ROOT / "data/synthetic/lesion_initial_thinning_composition/java_reference.json").read_text()
)
COMPETITION_REFERENCE = json.loads(
    (ROOT / "data/synthetic/lesion_competition_pass/java_reference.json").read_text()
)


def _bits(values: np.ndarray | np.float32) -> list[int] | int:
    array = np.asarray(values, dtype=np.float32)
    if array.ndim == 0:
        return int(array.view(np.uint32))
    return array.view(np.uint32).reshape(-1).astype(np.uint64).tolist()


def _digest(values: np.ndarray, dtype: str) -> str:
    return hashlib.sha256(np.asarray(values, dtype=np.dtype(dtype)).tobytes(order="C")).hexdigest()


def _item(item) -> list[int]:
    return [int(item.value.view(np.uint32)), item.x, item.y, item.z, item.label]


def _run_java(tmp_path: Path, nifti: Path, volume_bits: list[int]) -> dict:
    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac and java, "a JDK is required for the selected Java competition reference"
    classes = tmp_path / "classes"
    tmp_path.mkdir(parents=True, exist_ok=True)
    classes.mkdir()
    subprocess.run(
        [javac, "-d", str(classes), str(JAVA_UPSTREAM / "BinaryHeap4D.java"),
         str(JAVA_UPSTREAM / "CriticalPointLUT.java"), str(HARNESS)],
        check=True, capture_output=True, text=True,
    )
    lut_resource = classes / "edu/jhmi/rad/medic/structures/critical186LUT.raw.gz"
    lut_resource.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(LUT, lut_resource)
    volume_argument = ",".join(f"{value:08x}" for value in volume_bits)
    completed = subprocess.run(
        [java, "-cp", str(classes), "LesionToadsComposedCompetitionReferenceHarness",
         str(nifti), volume_argument],
        check=True, capture_output=True, text=True,
    )
    return json.loads(completed.stdout)


def test_membership_output_interoperates_with_java_referenced_alternative_competition(tmp_path: Path) -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == (
        "727df8bc9ce7f4f3b3429b94d7549fd8b110bcf28273af08443c8ffe1a775772"
    )
    assert hashlib.sha256((JAVA_UPSTREAM / "BinaryHeap4D.java").read_bytes()).hexdigest() == COMPETITION_REFERENCE["heap_sha256"]
    assert hashlib.sha256((JAVA_UPSTREAM / "CriticalPointLUT.java").read_bytes()).hexdigest() == COMPETITION_REFERENCE["lut_source_sha256"]
    assert hashlib.sha256(LUT.read_bytes()).hexdigest() == COMPETITION_REFERENCE["lut_asset_sha256"]
    # _compose runs the default linear_fast_marching branch on the committed
    # folded 17^3 RAS phantom and compares its first state to selected Java.
    # The competition adapter below is from the alternative queue mode; the
    # wrapper does not run the combined sequence tested here.
    frames, state, initial, thinning = _compose()
    result = run_post_thinning_competition_pass(
        thinning,
        state["images"], SPEC["class_centroids"], state["lesion"], state["black"],
        state["outlier"], SPEC["lesion_centroid"], (0.0, 0.0),
        state["shape_values"], SPEC["shape_priors"], SPEC["smoothing_priors"],
        state["direct_relations"], frames[..., 63], state["ventricle_distances"],
        state["gray_matter_distances"], state["brainstem_distances"], SPEC["object_types"],
        centroid_options={"mode": "none", "modalities": ("T1_MPRAGE", "FLAIR"),
                          "lesion_power_exponents": (0.4, 1.2)},
        channel_weight_options={"outlier_distance": SPEC["outlier_energy"]},
        outlier_energy=SPEC["outlier_energy"],
        lesion_power_exponent=SPEC["lesion_power_exponent"],
        lesion_class=SPEC["lesion_class"],
        flair_channel=1,
    )
    update = result.membership_update
    assert update.stage_order[-4:] == (
        "computeCentroids", "computeVariances", "computeChannelWeights",
        "computeMemberships:GMflag=true",
    )
    # P39am's pinned source-faithful NaN energy cells remain computed NaNs;
    # the corresponding memberships are finite and can enter competition.
    assert np.count_nonzero(np.isnan(update.membership.class_energies)) == 37_125
    assert np.isfinite(update.membership.class_memberships).all()
    assert _bits(update.statistics.channel_weights) == json.loads(
        (ROOT / "data/synthetic/lesion_post_thinning_statistics_refresh/java_reference.json").read_text()
    )["weights"]

    java_volume = INITIAL_REFERENCE["final_volume_bits"]
    assert _bits(thinning.thinning.run.class_volume) == java_volume
    assert result.energy.previous_energy == np.float32(initial.objective)
    assert _bits(result.energy.energy_ratio) == 0xFFC00000
    assert result.stage_order[-3:] == (
        "energyRatio", "initCompetitionLabels", "propagateCompetition"
    )

    competition = result.competition
    assert competition.initialization.inserted_items
    assert competition.remaining_queue == ()
    assert len(competition.events) > 0
    java = _run_java(tmp_path, _save_composed_fixture(tmp_path, thinning, update), java_volume)
    assert _bits(competition.class_volume) == java["volume_bits"]
    assert [_item(item) for item in competition.initialization.inserted_items] == java["initial_queue"]
    events = [
        {"pop": _item(event.popped_item), "processed": event.processed,
         "candidate": event.candidate_label, "accepted": event.topology_accepted,
         "inserted": [_item(item) for item in event.inserted_neighbors]}
        for event in competition.events
    ]
    assert events == java["events"]
    assert sum(event.processed and event.topology_accepted is True for event in competition.events) == java["accepted"]
    assert sum(event.processed and event.topology_accepted is False for event in competition.events) == java["rejected"]
    assert sum(not event.processed for event in competition.events) == java["stale"]
    assert _digest(competition.segmentation, "<i2") == java["final_segmentation_sha256"]
    assert _digest(competition.skeleton, "<i2") == java["final_skeleton_sha256"]
    assert _digest(competition.available, "u1") == java["availability_sha256"]
    assert _digest(competition.ordering, "<f4") == java["ordering_sha256"]
    assert _bits(competition.class_volume) == java["volume_bits"]


def _save_composed_fixture(tmp_path: Path, thinning, update) -> Path:
    run = thinning.thinning.run
    packed = np.concatenate(
        (run.segmentation[..., None].astype(np.float32),
         update.membership.class_memberships,
         update.membership.lesion_membership[..., None],
         update.membership.black_hole_membership[..., None],
         run.ordering[..., None]),
        axis=-1,
    ).astype(np.float32)
    original = nib.load(str(ROOT / "data/synthetic/lesion_initial_thinning_composition/initial_thinning_inputs_ras.nii.gz"))
    assert nib.aff2axcodes(original.affine) == ("R", "A", "S")
    path = tmp_path / "composed_competition_ras.nii.gz"
    nib.save(nib.Nifti1Image(packed, original.affine), str(path))
    return path
