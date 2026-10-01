from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest


TESTS = Path(__file__).parent
ROOT = TESTS / "data" / "synthetic" / "plus_lesion_atlas_selection"
ATLAS_DIR = TESTS.parents[0] / "data" / "atlases" / "lesionToads-atlas-2012"
HARNESS = TESTS / "java_reference" / "MaCRUISEplusLesionAtlasSelectionReferenceHarness.java"
REFERENCE_PATH = ROOT / "java_reference.json"
REFERENCE = json.loads(REFERENCE_PATH.read_text())


def _run_harness(tmp_path: Path, atlas_index: int, mprage: bool, flair: bool) -> tuple[str, str]:
    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac is not None and java is not None, "a JDK is required for the selected-Java reference"
    subprocess.run([javac, "-d", str(tmp_path), str(HARNESS)], check=True, capture_output=True, text=True)

    with_lesion = ATLAS_DIR / "cruise-atlas-12obj-lesiontoads2012.txt"
    no_lesion_t1_flair = ATLAS_DIR / "cruise-atlas-12obj-toads2012.txt"
    # The historical 10-object T1-only descriptor is packaged with the plugin,
    # but is not among the self-contained descriptors copied into this repo.
    no_lesion_t1_only = ROOT / "cruise-atlas-10obj-toads2010.txt"
    completed = subprocess.run(
        [
            java,
            "-cp",
            str(tmp_path),
            "MaCRUISEplusLesionAtlasSelectionReferenceHarness",
            str(atlas_index),
            str(mprage).lower(),
            str(flair).lower(),
            str(with_lesion),
            str(no_lesion_t1_flair),
            str(no_lesion_t1_only),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    selected_path, class_evidence = completed.stdout.rstrip("\n").split("\t", maxsplit=1)
    return Path(selected_path).name, class_evidence


def test_selected_plus_and_explicit_plugin_atlas_routes_match_java_excerpt(tmp_path: Path) -> None:
    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == REFERENCE["harness_sha256"]
    for descriptor_key, record in REFERENCE["bundled_descriptors"].items():
        if descriptor_key == "no_lesion_t1_only":
            continue
        descriptor = TESTS.parents[0] / record["path"]
        assert hashlib.sha256(descriptor.read_bytes()).hexdigest() == record["sha256"]

    expected = REFERENCE["cases"]
    selected_plus = _run_harness(tmp_path, atlas_index=0, mprage=True, flair=True)
    assert selected_plus == (
        expected["selected_plus_with_lesion"]["selected_descriptor"],
        expected["selected_plus_with_lesion"]["outer_ventricular_class"],
    )
    assert selected_plus[1] == "10|OuterVentricular-WM|25|lesion|hasLesion=true"

    plugin_no_lesion_t1_flair = _run_harness(tmp_path, atlas_index=1, mprage=True, flair=True)
    assert plugin_no_lesion_t1_flair == (
        expected["plugin_no_lesion_t1_flair"]["selected_descriptor"],
        expected["plugin_no_lesion_t1_flair"]["outer_ventricular_class"],
    )
    assert plugin_no_lesion_t1_flair[1] == "10|OuterVentricular-WM|25|obj|hasLesion=false"

    plugin_no_lesion_t1_only = _run_harness(tmp_path, atlas_index=1, mprage=True, flair=False)
    assert plugin_no_lesion_t1_only == (
        expected["plugin_no_lesion_t1_only"]["selected_descriptor"],
        "not-inspected",
    )


def test_plus_boundary_202_200_fixture_is_not_atlas_lesion_routing() -> None:
    descriptor = (ATLAS_DIR / "cruise-atlas-12obj-lesiontoads2012.txt").read_text()
    assert "OuterVentricular-WM\t25\tlesion" in descriptor
    for filename, expected_patch_value, ref_key in (
        ("segmentation_label202_overlook_patch.nii.gz", 202, "label_202_sha256"),
        ("segmentation_label200_control.nii.gz", 200, "label_200_sha256"),
    ):
        path = TESTS / "data" / "synthetic" / "plus_ace_branch" / filename
        assert hashlib.sha256(path.read_bytes()).hexdigest() == REFERENCE["boundary_fixture"][ref_key]
        image = nib.load(path)
        assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
        values = np.asarray(image.dataobj)
        assert np.any(values == expected_patch_value)
    assert 202 != 25


@pytest.mark.parametrize("atlas_index,mprage,flair", [(3, True, True), (1, False, True)])
def test_plugin_selection_keeps_source_default_and_t1_only_predicate(
    tmp_path: Path, atlas_index: int, mprage: bool, flair: bool
) -> None:
    # The plugin chooses the with-lesion path for unknown selection values;
    # No Lesion uses T1-only unless both MPRAGE and FLAIR are present.
    selected_name, class_evidence = _run_harness(tmp_path, atlas_index, mprage, flair)
    if atlas_index == 3:
        assert selected_name == REFERENCE["cases"]["selected_plus_with_lesion"]["selected_descriptor"]
        assert class_evidence.startswith("10|OuterVentricular-WM|25|lesion|")
    else:
        assert selected_name == REFERENCE["cases"]["plugin_no_lesion_t1_only"]["selected_descriptor"]
        assert class_evidence == "not-inspected"
