"""Rebuild the Java objective reference from the committed RAS NIfTI fixture."""

from pathlib import Path
import subprocess
import tempfile

import nibabel as nib
import numpy as np


HERE = Path(__file__).parent
HARNESS = Path(__file__).parents[3] / "java_reference/LesionToadsObjectiveReductionReferenceHarness.java"


def _java_float(value: np.float32) -> str:
    number = float(value)
    if np.isnan(number):
        return "NaN"
    if np.isposinf(number):
        return "Infinity"
    if np.isneginf(number):
        return "-Infinity"
    return format(number, ".9g")


def main() -> None:
    volume = np.asarray(nib.load(HERE / "objective_inputs_ras.nii.gz").dataobj, dtype=np.float32)
    with tempfile.TemporaryDirectory(prefix="lesion-objective-") as temp_name:
        temp = Path(temp_name)
        input_path = temp / "fixture.txt"
        input_path.write_text("\n".join(_java_float(value) for value in volume.reshape(-1)))
        classes = temp / "classes"
        classes.mkdir()
        subprocess.run(["javac", "-d", str(classes), str(HARNESS)], check=True)
        completed = subprocess.run(
            ["java", "-cp", str(classes), "LesionToadsObjectiveReductionReferenceHarness", str(input_path)],
            check=True,
            text=True,
            capture_output=True,
        )
        (HERE / "java_reference.json").write_text(completed.stdout)


if __name__ == "__main__":
    main()
