from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.plus_pipeline import PlusUncroppedTGDMResult, extract_plus_surfaces
from pymacruise.surfaces import extract_levelset_surface


FIXTURE = Path(__file__).parent / "data/synthetic/plus_surfaces"


def _closed_mesh_edges(mesh) -> bool:
    faces = np.asarray(mesh.faces)
    edges = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    edges.sort(axis=1)
    _unique, counts = np.unique(edges, axis=0, return_counts=True)
    return bool(np.all(counts == 2))


def test_plus_surface_adapter_routes_three_nested_ras_level_sets() -> None:
    inner = nib.load(FIXTURE / "inner.nii.gz")
    central = nib.load(FIXTURE / "central.nii.gz")
    outer = nib.load(FIXTURE / "outer.nii.gz")
    shape = inner.shape
    affine = inner.affine
    result = PlusUncroppedTGDMResult(
        initial_level_set=outer,
        tissue_sign=nib.Nifti1Image(np.zeros(shape, dtype=np.uint8), affine),
        initial_inner=inner,
        inner=inner,
        central=central,
        outer=outer,
    )

    actual = extract_plus_surfaces(result)
    for name in ("inner", "central", "outer"):
        mesh = getattr(actual, name)
        direct = extract_levelset_surface(getattr(result, name), level=0.0)
        np.testing.assert_array_equal(mesh.vertices, direct.vertices)
        np.testing.assert_array_equal(mesh.faces, direct.faces)
        assert mesh.coordinate_frame == "RAS"
        assert mesh.units == "mm"
        assert _closed_mesh_edges(mesh)

    center_ras = np.array([8.0, 12.0, 16.0])
    means = [
        np.linalg.norm(getattr(actual, name).vertices - center_ras, axis=1).mean()
        for name in ("inner", "central", "outer")
    ]
    assert means[0] < means[1] < means[2]
