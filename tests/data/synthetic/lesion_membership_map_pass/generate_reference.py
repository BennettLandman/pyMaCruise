"""Rebuild the selected-Java map reference from the committed NIfTI fixture."""

from pathlib import Path
import subprocess
import tempfile

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
HARNESS = Path(__file__).parents[3] / "java_reference/LesionToadsMembershipMapPassReferenceHarness.java"


def main() -> None:
    image = nib.load(HERE / "inputs_ras.nii.gz")
    values = np.asarray(image.dataobj, dtype=np.float32)
    with tempfile.TemporaryDirectory(prefix="lesion-map-pass-") as temp_name:
        temp = Path(temp_name)
        input_path = temp / "fixture.txt"
        np.savetxt(input_path, values.reshape(-1), fmt="%.9g")
        classes = temp / "classes"
        classes.mkdir()
        subprocess.run(["javac", "-d", str(classes), str(HARNESS)], check=True)
        completed = subprocess.run(
            ["java", "-cp", str(classes), "LesionToadsMembershipMapPassReferenceHarness", str(input_path)],
            check=True,
            text=True,
            capture_output=True,
        )
        (HERE / "java_reference.json").write_text(completed.stdout)


if __name__ == "__main__":
    main()
