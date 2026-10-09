"""Create a compact synthetic plus case matched to bundled LesionTOADS geometry."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).parent
ATLAS = ROOT.parents[3] / "data/atlases/lesionToads-atlas-2012"
DESCRIPTOR = ATLAS / "cruise-atlas-12obj-lesiontoads2012.txt"
TARGET = (65, 81, 62)
TOPOLOGY = "brain-atlas-cruise-12obj-lesions-topology_calc.nii.gz"


def _resample(values: np.ndarray, target: tuple[int, int, int], order: int) -> np.ndarray:
    coordinates = np.meshgrid(
        *(np.linspace(0, n - 1, m, dtype=np.float32) for n, m in zip(values.shape, target, strict=True)),
        indexing="ij",
    )
    return map_coordinates(values, coordinates, order=order, mode="nearest", prefilter=order > 1)


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    descriptor_text = DESCRIPTOR.read_text(encoding="utf-8")
    topo_image = nib.load(str(ATLAS / TOPOLOGY))
    source_topology = np.asarray(topo_image.dataobj, dtype=np.uint8)
    topology = _resample(source_topology, TARGET, order=0).astype(np.uint8)
    scale = (np.asarray(source_topology.shape, dtype=np.float64) - 1) / (np.asarray(TARGET) - 1)
    affine = topo_image.affine.copy()
    affine[:3, :3] = topo_image.affine[:3, :3] @ np.diag(scale)

    # Map source atlas regions to the explicit tissue labels used by the
    # project synthetic quality gates: 5=CSF, 44=WM, 100=GM.
    maseg = np.zeros(TARGET, dtype=np.uint16)
    maseg[topology == 5] = 5
    maseg[topology == 6] = 5
    maseg[np.isin(topology, [14, 15, 16, 17, 18, 22])] = 100
    maseg[np.isin(topology, [24, 25, 26, 27])] = 44
    brain = (topology != 1).astype(np.uint8)
    t1 = np.zeros(TARGET, dtype=np.float32)
    flair = np.full(TARGET, 80.0, dtype=np.float32)
    t1[maseg == 5], flair[maseg == 5] = 20.0, 90.0
    t1[maseg == 44], flair[maseg == 44] = 130.0, 180.0
    t1[maseg == 100], flair[maseg == 100] = 90.0, 150.0

    output_files: dict[str, np.ndarray] = {
        TOPOLOGY: topology,
        "t1.nii.gz": t1,
        "flair.nii.gz": flair,
        "brain-mask.nii.gz": brain,
        "maseg.nii.gz": maseg,
    }
    # Keep the bundled descriptor's exact names/classes/priors but place
    # compact atlas images beside the fixture and update the declared grid.
    referenced = sorted(set(re.findall(r"(?m)^File:\s*(\S+\.nii\.gz)\s*$", descriptor_text)))
    for name in referenced:
        source = np.asarray(nib.load(str(ATLAS / name)).dataobj, dtype=np.float32)
        output_files[name] = _resample(source, TARGET, order=1).astype(np.float32)
    for name, values in output_files.items():
        image = nib.Nifti1Image(values, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        image.set_data_dtype(values.dtype)
        nib.save(image, str(ROOT / name))

    spacing = np.linalg.norm(affine[:3, :3], axis=0)
    descriptor_text = descriptor_text.replace(str(ATLAS / TOPOLOGY), TOPOLOGY)
    descriptor_text = re.sub(r"(?m)^Dimensions:\s+172\s+216\s+164\s*$", "Dimensions:\t" + "\t".join(map(str, TARGET)), descriptor_text)
    descriptor_text = re.sub(r"(?m)^Resolutions:\s+0\.9375\s+0\.9375\s+0\.9375\s*$", "Resolutions:\t" + "\t".join(f"{v:.9g}" for v in spacing), descriptor_text)
    # The atlas descriptor parser resolves all NIfTI names relative to itself.
    (ROOT / DESCRIPTOR.name).write_text(descriptor_text, encoding="utf-8")
    records = {}
    for path in sorted(ROOT.glob("*.nii.gz")):
        image = nib.load(str(path))
        records[path.name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "shape": list(image.shape),
            "dtype": str(image.get_data_dtype()),
            "axis_codes": list(nib.aff2axcodes(image.affine)),
            "affine_ras_mm": image.affine.tolist(),
        }
    (ROOT / "provenance.json").write_text(json.dumps({
        "description": "compact source-atlas-aligned synthetic plus case",
        "source_atlas": "data/atlases/lesionToads-atlas-2012",
        "target_shape": list(TARGET),
        "label_mapping": {"0": "background", "5": "sulcal/ventricular CSF", "44": "white matter", "100": "gray matter"},
        "synthetic_intensity_values": {"T1": {"background": 0, "CSF": 20, "GM": 90, "WM": 130}, "FLAIR": {"background": 80, "CSF": 90, "GM": 150, "WM": 180}},
        "images": records,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
