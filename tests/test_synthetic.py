from __future__ import annotations

import hashlib
import json

import nibabel as nib
import numpy as np
from scipy import ndimage

from pymacruise.synthetic import make_synthetic_phantom, write_synthetic_fixtures


def test_phantom_is_deterministic_for_same_seed() -> None:
    first, first_affine, first_metadata = make_synthetic_phantom(seed=17)
    second, second_affine, second_metadata = make_synthetic_phantom(seed=17)

    assert first.keys() == second.keys()
    for key in first:
        np.testing.assert_array_equal(first[key], second[key])
    np.testing.assert_array_equal(first_affine, second_affine)
    assert first_metadata == second_metadata


def test_phantom_changes_noise_when_seed_changes() -> None:
    first, _, _ = make_synthetic_phantom(seed=17)
    second, _, _ = make_synthetic_phantom(seed=18)
    assert not np.array_equal(first["T1"], second["T1"])
    assert not np.array_equal(first["FLAIR"], second["FLAIR"])


def test_fixture_writer_is_byte_deterministic(tmp_path) -> None:
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"
    write_synthetic_fixtures(first_dir, seed=17)
    write_synthetic_fixtures(second_dir, seed=17)

    first_files = sorted(
        path.relative_to(first_dir) for path in first_dir.rglob("*") if path.is_file()
    )
    second_files = sorted(
        path.relative_to(second_dir) for path in second_dir.rglob("*") if path.is_file()
    )
    assert first_files == second_files
    for relative in first_files:
        first_hash = hashlib.sha256((first_dir / relative).read_bytes()).digest()
        second_hash = hashlib.sha256((second_dir / relative).read_bytes()).digest()
        assert first_hash == second_hash


def test_generated_nifti_fixtures_have_expected_metadata(tmp_path) -> None:
    write_synthetic_fixtures(tmp_path, seed=31)
    metadata = json.loads((tmp_path / "phantom.json").read_text())
    image = nib.load(str(tmp_path / "T1.nii.gz"))
    labels = nib.load(str(tmp_path / "MaSeg.nii.gz"))
    non_ras = nib.load(str(tmp_path / "non_ras_orientation.nii.gz"))

    assert metadata["provenance"].startswith("synthetic;")
    assert image.shape == (32, 28, 24)
    assert nib.aff2axcodes(image.affine) == ("R", "A", "S")
    assert labels.get_data_dtype() == np.dtype(np.uint16)
    assert set(np.unique(np.asanyarray(labels.dataobj))) == {0, 4, 100}
    assert nib.aff2axcodes(non_ras.affine) == ("L", "P", "S")


def test_edge_case_masks_and_ras_surface_are_present(tmp_path) -> None:
    write_synthetic_fixtures(tmp_path, seed=31)
    phantom = nib.load(str(tmp_path / "brain_mask.nii.gz"))
    cavity = nib.load(str(tmp_path / "cavity_mask.nii.gz"))
    labels = nib.load(str(tmp_path / "disconnected_labels.nii.gz"))
    surface = np.load(tmp_path / "brain_surface_ras.npz")

    brain_data = np.asanyarray(phantom.dataobj)
    cavity_data = np.asanyarray(cavity.dataobj)
    label_data = np.asanyarray(labels.dataobj)
    assert np.all(brain_data[cavity_data.astype(bool)] == 0)
    _, label_components = ndimage.label(label_data == 100)
    assert label_components == 2
    assert nib.aff2axcodes(phantom.affine) == ("R", "A", "S")

    vertices = surface["vertices_ras_mm"]
    faces = surface["faces"]
    edges = np.sort(
        np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])),
        axis=1,
    )
    unique_edges, counts = np.unique(edges, axis=0, return_counts=True)
    assert vertices.ndim == 2 and vertices.shape[1] == 3
    assert faces.ndim == 2 and faces.shape[1] == 3
    assert np.all(counts == 2), f"mesh has non-manifold/open edges: {unique_edges[counts != 2][:5]}"


def test_fixture_images_are_compressed_nifti_and_hash_manifest_matches() -> None:
    from pathlib import Path

    fixture_dir = Path(__file__).parent / "data" / "synthetic"
    manifest = fixture_dir / "SHA256SUMS"
    entries = [line.split("  ", maxsplit=1) for line in manifest.read_text().splitlines()]

    for path in fixture_dir.iterdir():
        if path.suffix == ".gz":
            assert path.name.endswith(".nii.gz")
    assert len(entries) == len([path for path in fixture_dir.iterdir() if path.is_file()]) - 1
    for digest, filename in entries:
        assert hashlib.sha256((fixture_dir / filename).read_bytes()).hexdigest() == digest
