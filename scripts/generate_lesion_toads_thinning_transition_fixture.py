"""Write the committed synthetic 3x3x3 LesionTOADS transition input."""

from pathlib import Path

import nibabel as nib
import numpy as np


ROOT = Path(__file__).resolve().parents[1] / "tests" / "data" / "synthetic" / "lesion_toads_thinning_transition"
CASES = (
    ("direct_accept", "212221122120200120112112212", 0),
    ("direct_reject", "221001200120102002101100002", 0),
    ("fallback_accept", "222222101122202022120200001", 0),
    ("no_connected_target", "000000000000000000000000000", 0),
    ("best_already", "000000000000000000000000000", 0),
    ("critical_volume", "000000000000000000000000000", 0),
    ("lesion_volume_max", "111111111111111111111111111", 1),
)


def main() -> None:
    packed = np.zeros((3, 3, 3, len(CASES) * 11), dtype=np.float32)
    coordinates = np.indices((3, 3, 3))
    face = (
        ((coordinates[0] == 0) | (coordinates[0] == 2))
        & (coordinates[1] == 1)
        & (coordinates[2] == 1)
    ) | (
        ((coordinates[1] == 0) | (coordinates[1] == 2))
        & (coordinates[0] == 1)
        & (coordinates[2] == 1)
    ) | (
        ((coordinates[2] == 0) | (coordinates[2] == 2))
        & (coordinates[0] == 1)
        & (coordinates[1] == 1)
    )
    for case_index, (name, encoded, processing_class) in enumerate(CASES):
        base = case_index * 11
        segmentation = np.fromiter((int(char) for char in encoded), dtype=np.float32).reshape(3, 3, 3)
        packed[..., base] = segmentation
        packed[..., base + 1] = segmentation
        memberships = np.array([0.125, 0.25, 0.0625], dtype=np.float32)
        if name == "best_already":
            memberships[0] = 0.75
        elif name == "lesion_volume_max":
            memberships[:] = [0.0, 0.125, 0.0]
        for class_index, value in enumerate(memberships):
            packed[..., base + 2 + class_index] = value
        packed[..., base + 5] = 0.125
        packed[..., base + 6] = 0.0625
        if name == "lesion_volume_max":
            packed[1, 1, 1, base + 5] = 0.5
            packed[1, 1, 1, base + 6] = 0.375
        if name.startswith("direct") or name.startswith("fallback"):
            packed[..., base + 8] = np.where(face, 0.375, 0.0)
        # The last feature is a center marker used to verify queue ordering.
        packed[1, 1, 1, base + 10] = -0.625
        if name == "critical_volume":
            # NIfTI feature values include only image inputs; class volume
            # thresholds are specified in fixture_metadata.json.
            pass
        if processing_class not in range(3):
            raise AssertionError("bad processing label")

    image = nib.Nifti1Image(packed, np.eye(4))
    image.set_qform(np.eye(4), code=1)
    image.set_sform(np.eye(4), code=1)
    nib.save(image, str(ROOT / "transition_inputs.nii.gz"))


if __name__ == "__main__":
    main()
