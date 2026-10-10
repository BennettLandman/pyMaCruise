"""Create a compact synthetic plus case matched to bundled LesionTOADS geometry."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import argparse

import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).parent
ATLAS = ROOT.parents[3] / "data/atlases/lesionToads-atlas-2012"
DESCRIPTOR = ATLAS / "cruise-atlas-12obj-lesiontoads2012.txt"
TARGET = (65, 81, 62)
TOPOLOGY = "brain-atlas-cruise-12obj-lesions-topology_calc.nii.gz"
TOPOLOGY_TO_MASEG = {
    5: 5,
    6: 5,
    15: 100,
    25: 44,
    26: 44,
    27: 44,
    14: 75,
    16: 75,
    17: 75,
    18: 75,
    22: 75,
    24: 75,
}


def _resample(values: np.ndarray, target: tuple[int, int, int], order: int) -> np.ndarray:
    coordinates = np.meshgrid(
        *(np.linspace(0, n - 1, m, dtype=np.float32) for n, m in zip(values.shape, target, strict=True)),
        indexing="ij",
    )
    return map_coordinates(values, coordinates, order=order, mode="nearest", prefilter=order > 1)


def _maseg_from_topology(topology: np.ndarray) -> np.ndarray:
    """Map only cerebral cortex/white matter to the cortical MaSeg labels."""
    maseg = np.zeros(topology.shape, dtype=np.uint16)
    for topology_label, maseg_label in TOPOLOGY_TO_MASEG.items():
        maseg[topology == topology_label] = maseg_label
    return maseg


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    parser.add_argument("--shape", type=int, nargs=3, default=TARGET, metavar=("NX", "NY", "NZ"))
    args = parser.parse_args()
    output_root = args.output_dir.expanduser().resolve()
    target = tuple(args.shape)
    if any(size < 3 for size in target):
        parser.error("each target dimension must be at least 3 voxels")
    output_root.mkdir(parents=True, exist_ok=True)
    descriptor_text = DESCRIPTOR.read_text(encoding="utf-8")
    topo_image = nib.load(str(ATLAS / TOPOLOGY))
    source_topology = np.asarray(topo_image.dataobj, dtype=np.uint8)
    topology = _resample(source_topology, target, order=0).astype(np.uint8)
    scale = (np.asarray(source_topology.shape, dtype=np.float64) - 1) / (np.asarray(target) - 1)
    affine = topo_image.affine.copy()
    affine[:3, :3] = topo_image.affine[:3, :3] @ np.diag(scale)

    # Preserve the descriptor's cortical distinction: Cerebrum-GM (15) is
    # cortical GM; cerebral WM classes (25-27) are cortical WM. Cerebellar
    # and subcortical/brainstem classes remain additional brain tissue (75),
    # rather than being mislabeled as cortical GM/WM in the surface checks.
    maseg = _maseg_from_topology(topology)
    brain = (topology != 1).astype(np.uint8)
    t1 = np.zeros(target, dtype=np.float32)
    flair = np.full(target, 80.0, dtype=np.float32)
    t1[np.isin(topology, [5, 6])], flair[np.isin(topology, [5, 6])] = 20.0, 90.0
    t1[np.isin(topology, [24, 25, 26, 27])], flair[np.isin(topology, [24, 25, 26, 27])] = 130.0, 180.0
    t1[np.isin(topology, [14, 15, 16, 17, 18, 22])], flair[np.isin(topology, [14, 15, 16, 17, 18, 22])] = 90.0, 150.0

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
        if name == TOPOLOGY:
            continue
        source = np.asarray(nib.load(str(ATLAS / name)).dataobj, dtype=np.float32)
        output_files[name] = _resample(source, target, order=1).astype(np.float32)
    for name, values in output_files.items():
        image = nib.Nifti1Image(values, affine)
        image.set_qform(affine, code=1)
        image.set_sform(affine, code=1)
        image.set_data_dtype(values.dtype)
        nib.save(image, str(output_root / name))

    spacing = np.linalg.norm(affine[:3, :3], axis=0)
    descriptor_text = descriptor_text.replace(str(ATLAS / TOPOLOGY), TOPOLOGY)
    descriptor_text = re.sub(r"(?m)^Dimensions:\s+172\s+216\s+164\s*$", "Dimensions:\t" + "\t".join(map(str, target)), descriptor_text)
    descriptor_text = re.sub(r"(?m)^Resolutions:\s+0\.9375\s+0\.9375\s+0\.9375\s*$", "Resolutions:\t" + "\t".join(f"{v:.9g}" for v in spacing), descriptor_text)
    # The atlas descriptor parser resolves all NIfTI names relative to itself.
    (output_root / DESCRIPTOR.name).write_text(descriptor_text, encoding="utf-8")
    records = {}
    for path in sorted(output_root.glob("*.nii.gz")):
        image = nib.load(str(path))
        records[path.name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "shape": list(image.shape),
            "dtype": str(image.get_data_dtype()),
            "axis_codes": list(nib.aff2axcodes(image.affine)),
            "affine_ras_mm": image.affine.tolist(),
        }
    (output_root / "provenance.json").write_text(json.dumps({
        "description": "source-atlas-aligned synthetic plus case",
        "source_atlas": "data/atlases/lesionToads-atlas-2012",
        "target_shape": list(target),
        "label_mapping": {"0": "background", "5": "CSF", "44": "cerebral white matter", "75": "additional noncortical brain tissue", "100": "cerebral gray matter"},
        "topology_to_maseg": {str(key): value for key, value in TOPOLOGY_TO_MASEG.items()},
        "synthetic_intensity_values": {"T1": {"background": 0, "CSF": 20, "GM": 90, "WM": 130}, "FLAIR": {"background": 80, "CSF": 90, "GM": 150, "WM": 180}},
        "images": records,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
