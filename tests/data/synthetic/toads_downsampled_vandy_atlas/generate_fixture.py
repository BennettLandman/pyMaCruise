"""Create a coarse whole-brain folded-cortex phantom from the bundled atlas.

The generator downsamples the repository's source atlas maps to a compact grid
so a full TOADS loop remains practical. It does not read or write clinical
data. Generated voxel sizes preserve the atlas field of view.
"""

from __future__ import annotations

from pathlib import Path
import nibabel as nib
import numpy as np
from scipy.ndimage import zoom


ROOT = Path(__file__).parent
REPOSITORY = ROOT.parents[3]
SOURCE_DESCRIPTOR = REPOSITORY / "data/atlases/cruise-atlas-10obj-toads2008.txt"
TARGET_SHAPE = (36, 45, 34)


def main() -> None:
    lines = SOURCE_DESCRIPTOR.read_text(encoding="utf-8").splitlines()
    # Resolve maps from the selected source descriptor, while keeping atlas
    # names and all intensity/variance/exponent tables intact.
    structure_rows: list[tuple[str, int, str]] = []
    structure_start = next(i for i, line in enumerate(lines) if line.startswith("Structures:"))
    topology_start = next(i for i, line in enumerate(lines) if line.startswith("Topology Atlas"))
    for line in lines[structure_start + 1 : topology_start]:
        if line.strip():
            name, label, kind = line.split()
            structure_rows.append((name, int(label), kind))
    topology_file = _file_after(lines, "Topology Atlas", start=topology_start)
    topology_resolution = next(
        line for line in lines[topology_start : topology_start + 8]
        if line.startswith("Resolutions:")
    )
    source_spacing = np.asarray(
        [float(value) for value in topology_resolution.partition(":")[2].split()],
        dtype=np.float64,
    )
    shape_start = next(i for i, line in enumerate(lines) if line.startswith("Shape Atlas:"))
    shape_files: dict[str, str] = {}
    for index, line in enumerate(lines[shape_start + 1 :]):
        if line.startswith("Structure:"):
            name = line.partition(":")[2].strip()
            shape_files[name] = _file_after(lines, "Structure:", start=shape_start + 1 + index)

    topology_source = nib.as_closest_canonical(nib.load(str(SOURCE_DESCRIPTOR.parent / topology_file)))
    source_shape = np.asarray(topology_source.shape, dtype=np.int32)
    target_spacing = source_shape * source_spacing / np.asarray(TARGET_SHAPE)
    target_affine = np.diag((*target_spacing, 1.0))
    target_affine[:3, 3] = -(np.asarray(TARGET_SHAPE) - 1) / 2 * target_spacing
    zooms = np.asarray(TARGET_SHAPE, dtype=np.float64) / source_shape

    # Source atlas voxel arrays share the Java descriptor grid even though a
    # historical topology NIfTI header carries 1 mm geometry. Resample the
    # RAS arrays by index grid, not by their mutually inconsistent headers.
    topology_values = zoom(
        np.asanyarray(topology_source.dataobj), zooms, order=0, prefilter=False
    ).astype(np.int16)
    _save("topology.nii.gz", topology_values, target_affine)
    descriptor: list[str] = [lines[0], "", f"Structures:\t{len(structure_rows)}"]
    descriptor.extend(f"{name}\t{label}\t{kind}" for name, label, kind in structure_rows)
    resolution = "\t".join(f"{value:.9g}" for value in target_spacing)
    dimensions = "\t".join(map(str, TARGET_SHAPE))
    descriptor.extend(("", "Topology Atlas:", "File:\ttopology.nii.gz", f"Dimensions:\t{dimensions}", f"Resolutions:\t{resolution}", ""))

    original_shape_line = lines[shape_start].split(":", 1)[1].strip()
    sample_count = int(original_shape_line)
    descriptor.extend((f"Shape Atlas:\t{sample_count}", f"Dimensions:\t{dimensions}", f"Resolutions:\t{resolution}", "Orientations:\t1\t1\t6\t3"))
    for row, (name, _label, _kind) in enumerate(structure_rows):
        filename = f"shape-{row}.nii.gz"
        source = nib.as_closest_canonical(nib.load(str(SOURCE_DESCRIPTOR.parent / shape_files[name])))
        sampled = zoom(
            np.asanyarray(source.dataobj), zooms, order=1, prefilter=True
        ).astype(np.float32)
        _save(filename, sampled, target_affine)
        descriptor.extend((f"Structure: {name}", f"File:\t{filename}", f"Dimensions:\t{dimensions}", f"Resolutions:\t{resolution}", "Center:\t0\t0\t0"))

    # Keep selected descriptor tables, which encode the source atlas's class
    # intensity priors, variance weighting, and optimized factors.
    for heading in ("Intensity Atlas:", "Intensity Variance Factor:", "Optimized factors:"):
        start = next(i for i, line in enumerate(lines) if line.startswith(heading))
        end = next((i for i in range(start + 1, len(lines)) if not lines[i].strip()), len(lines))
        descriptor.extend(("", *lines[start:end]))
    (ROOT / "atlas.txt").write_text("\n".join(descriptor) + "\n", encoding="utf-8")

    labels = topology_values
    brain = labels != 1
    # Atlas class intensities are deterministic and follow the descriptor's
    # T1_MPRAGE class proportions with a low-amplitude smooth spatial bias.
    mprage_line = next(line for line in lines if line.startswith("T1_MPRAGE"))
    prior_values = np.asarray([float(value) for value in mprage_line.split()[1:]], dtype=np.float32)
    class_labels = np.asarray([label for _name, label, _kind in structure_rows], dtype=np.int16)
    intensities = {int(label): float(prior_values[index]) for index, label in enumerate(class_labels)}
    image = np.zeros(TARGET_SHAPE, dtype=np.float32)
    for label in np.unique(labels):
        image[labels == label] = np.float32(700.0 * intensities.get(int(label), 0.5) + 80.0)
    x, y, z = np.indices(TARGET_SHAPE, dtype=np.float32)
    image[brain] += np.float32(0.03) * x[brain] + np.float32(0.02) * y[brain] + np.float32(0.01) * z[brain]
    image[~brain] = np.float32(0.0)
    _save("t1.nii.gz", image, target_affine)
    _save("brain-mask.nii.gz", brain.astype(np.uint8), target_affine)
    maseg = np.zeros(TARGET_SHAPE, dtype=np.uint8)
    for label in (14, 15):
        maseg[labels == label] = 100  # selected GM label used by MaCRUISE cleanup
    for label in (22, 24, 25):
        maseg[labels == label] = 4  # selected WM label used by MaCRUISE cleanup
    _save("maseg.nii.gz", maseg, target_affine)
    summary = (
        "# Downsampled bundled VandyCRUISE atlas synthetic phantom\n\n"
        "This deterministic, coarse whole-brain phantom is generated from the atlas maps already bundled in this repository. "
        f"Its {TARGET_SHAPE[0]}×{TARGET_SHAPE[1]}×{TARGET_SHAPE[2]} grid preserves the source field of view and atlas labels; "
        "a synthetic T1 is synthesized from the descriptor's T1_MPRAGE class proportions with a small smooth bias. "
        "All stored images are compressed NIfTI in RAS orientation. Regenerate with `python generate_fixture.py`.\n\n"
        "This is not clinical data and does not establish full-loop source equivalence or clinical performance. The coarse "
        "spacing is intended to make a folded whole-brain convergence/runtime experiment practical.\n"
    )
    (ROOT / "README.md").write_text(summary, encoding="utf-8")


def _file_after(lines: list[str], section: str, *, start: int = 0) -> str:
    for line in lines[start : start + 10]:
        if line.startswith("File:"):
            return line.partition(":")[2].strip()
    raise ValueError(f"missing File entry after {section}")


def _save(name: str, values: np.ndarray, affine: np.ndarray) -> None:
    image = nib.Nifti1Image(values, affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    nib.save(image, str(ROOT / name))


if __name__ == "__main__":
    main()
