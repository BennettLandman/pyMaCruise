from __future__ import annotations

import json
import hashlib
from pathlib import Path
import subprocess

import nibabel as nib
import numpy as np

from pymacruise.toads_thinning import toads_propagate_thinning, toads_thinning_voxel_transition
from test_toads_p18czd_next_loop_thinning import _float_bits, _labels, _serialize

ROOT = Path(__file__).parent
CASE_DIR = ROOT / "data/synthetic/toads_p18czm_nonzero_shape_thinning"
CASE = json.loads((CASE_DIR / "case.json").read_text())
THINNING_SCRIPT = ROOT.parent / "scripts/run_java_toads_post_registration_next_thinning_reference.sh"


def test_selected_thinning_nonzero_shape_prior_competition_matches_java(tmp_path: Path) -> None:
    segmentation_image = nib.load(str(CASE_DIR / CASE["segmentation"]))
    memberships_image = nib.load(str(CASE_DIR / CASE["memberships"]))
    shapes_image = nib.load(str(CASE_DIR / CASE["shape_priors"]))
    for name, expected_sha in CASE["fixture_sha256"].items():
        assert hashlib.sha256((CASE_DIR / name).read_bytes()).hexdigest() == expected_sha
    for image in (segmentation_image, memberships_image, shapes_image):
        assert nib.aff2axcodes(image.affine) == tuple(CASE["source_axis_codes"])
    segmentation = np.asarray(segmentation_image.dataobj, dtype=np.int16)
    memberships = np.asarray(memberships_image.dataobj, dtype=np.float32)
    shapes = np.asarray(shapes_image.dataobj, dtype=np.float32)
    assert memberships.shape == shapes.shape == (*segmentation.shape, 11)
    assert np.isfinite(memberships).all() and np.isfinite(shapes).all()
    assert np.count_nonzero(shapes) == shapes.size

    adjacency = np.zeros((11, 11), dtype=np.bool_)
    adjacency[0, 2] = True
    adjacency[3, 2] = True
    x, y, z = CASE["processing_voxel_xyz"]
    old = int(segmentation[x, y, z])
    local = (slice(x - 1, x + 2), slice(y - 1, y + 2), slice(z - 1, z + 2))
    volumes = np.full(11, 200.0, dtype=np.float32)
    minimums = np.full(11, 100.0, dtype=np.float32)
    shaped_transition = toads_thinning_voxel_transition(
        segmentation[local], segmentation[local], memberships[x, y, z], adjacency,
        shapes[local], volumes, minimums, processing_label=old,
    )
    zero_transition = toads_thinning_voxel_transition(
        segmentation[local], segmentation[local], memberships[x, y, z], adjacency,
        np.zeros_like(shapes[local]), volumes, minimums, processing_label=old,
    )
    assert shaped_transition.candidate_label == 3
    assert zero_transition.candidate_label == 0
    # Both proposed changes are rejected by the source topology gate; the
    # distinct candidate labels prove the nonzero priors control fallback.
    assert shaped_transition.topology_accepted is False
    assert zero_transition.topology_accepted is False

    result = toads_propagate_thinning(
        segmentation, memberships, adjacency, shapes,
        volume_ratio=CASE["volume_ratio"], limit=CASE["thinning_limit"],
    )
    assert any((item.x, item.y, item.z) == (x, y, z) for item in result.processed_items)
    input_path = tmp_path / "p18czm_nonzero_shape_input.bin"
    input_path.write_bytes(_serialize(segmentation, memberships, shapes, adjacency))
    java = json.loads(subprocess.run(
        [
            "sh", str(THINNING_SCRIPT), str(input_path),
            *(str(v) for v in segmentation.shape), str(memberships.shape[-1]),
        ],
        check=True, capture_output=True, text=True,
    ).stdout)
    assert java["source"] == CASE["source"]
    assert java["source_sha256"] == "63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884"
    assert java["shape"] == [*segmentation.shape, memberships.shape[-1]]
    assert java["volume_ratio_bits"] == int(np.asarray([np.float32(CASE["volume_ratio"])]).view(np.int32)[0])
    assert java["limit_bits"] == int(np.asarray([np.float32(CASE["thinning_limit"])]).view(np.int32)[0])

    initial_volumes = np.zeros(memberships.shape[-1], dtype=np.float32)
    for i in range(1, segmentation.shape[0] - 1):
        for j in range(1, segmentation.shape[1] - 1):
            for k in range(1, segmentation.shape[2] - 1):
                klass = int(segmentation[i, j, k])
                initial_volumes[klass] = np.float32(initial_volumes[klass] + memberships[i, j, k, klass])
    assert java["initial_volumes_bits"] == _float_bits(initial_volumes)
    assert java["minimum_volumes_bits"] == _float_bits(result.minimum_class_volumes)
    assert java["boundary_counts"] == result.boundary_counts.astype(int).tolist()
    assert java["processed"] == [
        [int(np.asarray([item.value], dtype=np.float32).view(np.int32)[0]), item.x, item.y, item.z, item.label]
        for item in result.processed_items
    ]
    assert java["segmentation"] == _labels(result.segmentation)
    assert java["skeleton"] == _labels(result.skeleton, empty=True)
    assert java["available"] == "".join("1" if value else "0" for value in result.available.ravel())
    assert java["final_volumes_bits"] == _float_bits(result.class_volumes)
