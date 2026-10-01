from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.lesion_initial_alignment_bounds import (
    compute_lesion_initial_alignment_bounds,
)
from pymacruise.lesion_initial_atlas_alignment import initialize_lesion_atlas_alignment


SYNTHETIC = Path(__file__).parent / "data" / "synthetic"
SOURCE_FIXTURE = SYNTHETIC / "lesion_initial_atlas_alignment"
BOUNDS_FIXTURE = SYNTHETIC / "lesion_initial_alignment_bounds"
HARNESS = Path(__file__).parent / "java_reference" / "LesionInitialAlignmentBoundsReferenceHarness.java"
HARNESS_SHA = "fc92ed5e82261ad3d24c8ac854661c1c9af8aa05a6b80d9a130fdcc6ee2aa26c"
GENERATOR_SHA = "aa965111379700a6b149fab8d03f99ef1f1986bd1c0dab96c5549e797b76dd00"
PARAMETRIC_TRANSFORM = (
    Path(__file__).parent
    / "java_reference"
    / "upstream"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
    / "methods"
    / "ParametricTransform.java"
)
ROTATION_MATRIX = (
    Path(__file__).parent
    / "java_reference"
    / "upstream"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
    / "structures"
    / "RotationMatrix.java"
)
IMAGE_FUNCTIONS = (
    Path(__file__).parent
    / "java_reference"
    / "upstream"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
    / "libraries"
    / "ImageFunctions.java"
)
NUMERICS = (
    Path(__file__).parent
    / "java_reference"
    / "upstream"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
    / "utilities"
    / "Numerics.java"
)
MEDIC_UTIL_STUB = (
    Path(__file__).parent
    / "java_reference"
    / "stubs"
    / "edu"
    / "jhmi"
    / "rad"
    / "medic"
    / "utilities"
    / "MedicUtilPublic.java"
)


def _bits(values: np.ndarray) -> list[str]:
    return [f"{int(value):08x}" for value in np.asarray(values).view(np.uint32).ravel()]


def test_initial_alignment_refreshes_bounded_shape_support_from_final_transform(
    tmp_path: Path,
) -> None:
    source_metadata = json.loads((SOURCE_FIXTURE / "inputs.json").read_text())
    fixture_metadata = json.loads((BOUNDS_FIXTURE / "inputs.json").read_text())
    reference = json.loads((BOUNDS_FIXTURE / "java_reference.json").read_text())
    initial_reference = json.loads((SOURCE_FIXTURE / "java_reference.json").read_text())

    assert hashlib.sha256(HARNESS.read_bytes()).hexdigest() == HARNESS_SHA
    generator = BOUNDS_FIXTURE / "generate_fixture.py"
    assert hashlib.sha256(generator.read_bytes()).hexdigest() == GENERATOR_SHA
    image_obj = nib.load(str(SOURCE_FIXTURE / source_metadata["input_nifti"]))
    assert nib.aff2axcodes(image_obj.affine) == ("R", "A", "S")
    channels_first = np.transpose(np.asarray(image_obj.dataobj, dtype=np.float32), (3, 0, 1, 2))
    shapes = []
    for filename in fixture_metadata["shape_niftis"]:
        path = BOUNDS_FIXTURE / filename
        assert hashlib.sha256(path.read_bytes()).hexdigest() == fixture_metadata["shape_sha256"][filename]
        shape_obj = nib.load(str(path))
        assert nib.aff2axcodes(shape_obj.affine) == ("R", "A", "S")
        shapes.append(np.asarray(shape_obj.dataobj, dtype=np.float32))

    centroids = np.asarray(source_metadata["centroids_channels_classes"], dtype=np.float32)
    preprocessing_transform = np.asarray(source_metadata["preprocessing_transform"], dtype=np.float32)
    originals = (channels_first.copy(), centroids.copy(), preprocessing_transform.copy(), [s.copy() for s in shapes])
    alignment = initialize_lesion_atlas_alignment(
        channels_first,
        centroids,
        shapes,
        preprocessing_transform,
        crop_bounds=tuple(tuple(v) for v in source_metadata["crop_bounds_source_coordinates"]),
        pad_border=source_metadata["pad_border"],
        image_resolutions_mm=tuple(source_metadata["image_resolutions_mm"]),
        template_resolutions_mm=tuple(source_metadata["template_resolutions_mm"]),
        init_align_iterations=source_metadata["init_align_iterations"],
        pyramid_levels=source_metadata["pyramid_levels"],
    )

    expected_memberships = np.asarray(initial_reference["membership_bits"], dtype=np.int32)
    expected_memberships = expected_memberships.view(np.float32).reshape(
        alignment.memberships.memberships.shape
    )
    np.testing.assert_array_equal(
        alignment.memberships.memberships.view(np.int32), expected_memberships.view(np.int32)
    )
    assert alignment.call_order == tuple(initial_reference["call_order"])
    assert [level.scale for level in alignment.registration.levels] == [8, 4, 2]
    assert _bits(alignment.registration.transform) == reference["transform_bits"]
    for before, after in zip(originals[:3], (channels_first, centroids, preprocessing_transform), strict=True):
        np.testing.assert_array_equal(after, before)
    for before, after in zip(originals[3], shapes, strict=True):
        np.testing.assert_array_equal(after, before)

    bounds = compute_lesion_initial_alignment_bounds(
        alignment,
        shapes,
        image_resolutions_mm=tuple(source_metadata["image_resolutions_mm"]),
        template_resolutions_mm=tuple(source_metadata["template_resolutions_mm"]),
    )
    assert bounds.call_order == tuple(reference["call_order"])
    assert _bits(bounds.transform) == reference["transform_bits"]
    assert [
        {
            "class": class_index,
            "lower_xyz": list(item.lower_xyz),
            "upper_xyz": list(item.upper_xyz),
            "support_count": item.support_count,
        }
        for class_index, item in enumerate(bounds.shape_bounds)
    ] == reference["bounds"]

    javac = shutil.which("javac")
    java = shutil.which("java")
    assert javac is not None and java is not None, "a JDK is required for the selected-Java reference"
    subprocess.run(
        [
            javac,
            "-d",
            str(tmp_path),
            str(MEDIC_UTIL_STUB),
            str(NUMERICS),
            str(ROTATION_MATRIX),
            str(PARAMETRIC_TRANSFORM),
            str(IMAGE_FUNCTIONS),
            str(HARNESS),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    java_result = subprocess.run(
        [java, "-cp", str(tmp_path), "LesionInitialAlignmentBoundsReferenceHarness", *_bits(bounds.transform)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert json.loads(java_result.stdout) == reference
