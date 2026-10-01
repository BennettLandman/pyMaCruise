from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import zipfile


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ATLAS_CASES = {
    "cruise-atlas-10obj-toads2008.txt": ("cruise-atlas-10obj", "atlas_CerebralGM_cruise.nii.gz"),
    "lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt": (
        "lesionToads-atlas-2012", "brain-atlas-cruise-12obj-lesions-topology_calc.nii.gz"
    ),
    "lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt": (
        "lesionToads-atlas-2012", "atlas_Cerebrum-GM.nii.gz"
    ),
}


def test_sdist_and_wheel_bundle_resolvable_atlas_resources(tmp_path: Path) -> None:
    # Isolate all setuptools/build outputs from the live checkout. Hardlinks
    # avoid copying the large checked-in atlas tree when temp and source share
    # a filesystem; copy2 is the cross-device fallback.
    source = tmp_path / "source"
    source.mkdir()
    for filename in ("pyproject.toml", "setup.py", "MANIFEST.in", "readme.md"):
        _link_or_copy(PROJECT_ROOT / filename, source / filename)
    shutil.copytree(PROJECT_ROOT / "src", source / "src", copy_function=_link_or_copy)
    shutil.copytree(PROJECT_ROOT / "data" / "atlases", source / "data" / "atlases", copy_function=_link_or_copy)

    first_dist = tmp_path / "source-dist"
    first_dist.mkdir()
    _run_build(source, first_dist, "--sdist", "--wheel")
    sdist = next(first_dist.glob("*.tar.gz"))
    source_wheel = next(first_dist.glob("*.whl"))

    extracted_sdist = tmp_path / "sdist"
    extracted_sdist.mkdir()
    with tarfile.open(sdist, "r:gz") as archive:
        if hasattr(tarfile, "data_filter"):
            archive.extractall(extracted_sdist, filter="data")
        else:  # Python 3.10 and 3.11 do not expose extraction filters.
            archive.extractall(extracted_sdist)
    sdist_root = next(extracted_sdist.iterdir())
    sdist_atlases = sdist_root / "data" / "atlases"
    _assert_sdist_atlas_bundle(sdist_atlases)

    second_dist = tmp_path / "sdist-wheel"
    second_dist.mkdir()
    _run_build(sdist_root, second_dist, "--wheel")
    sdist_wheel = next(second_dist.glob("*.whl"))

    _assert_wheel_bundle(source_wheel, tmp_path / "wheel-from-source")
    _assert_wheel_bundle(sdist_wheel, tmp_path / "wheel-from-sdist")
    _assert_isolated_wheel_install(sdist_wheel, tmp_path / "installed-target")


def _link_or_copy(source: str | os.PathLike[str], destination: str | os.PathLike[str]) -> str:
    try:
        os.link(source, destination)
        return os.fspath(destination)
    except OSError:
        return shutil.copy2(source, destination)


def _run_build(source: Path, output: Path, *targets: str) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "build", "--no-isolation", "--outdir", str(output), *targets],
        cwd=source,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"build failed:\n{result.stdout}\n{result.stderr}"


def _assert_sdist_atlas_bundle(root: Path) -> None:
    assert root.is_dir()
    for descriptor, (asset_directory, asset_name) in ATLAS_CASES.items():
        descriptor_path = root / descriptor
        assert descriptor_path.is_file(), descriptor
        assert (root / asset_directory / asset_name).is_file(), asset_name
        _assert_descriptor_references_are_bundled(root, descriptor_path)


def _assert_wheel_bundle(wheel: Path, extract_to: Path) -> None:
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        atlas_root = "pymacruise/data/atlases/"
        for descriptor, (asset_directory, asset_name) in ATLAS_CASES.items():
            assert atlas_root + descriptor in names, descriptor
            assert atlas_root + asset_directory + "/" + asset_name in names, asset_name
            assert not archive.read(atlas_root + descriptor).startswith(
                b"version https://git-lfs.github.com/spec/v1\n"
            )
        archive.extractall(extract_to)

    # Exercise the package resource resolver and the public loaders from the
    # unpacked wheel tree without installing the wheel or consulting checkout
    # paths.
    script = r"""
from pathlib import Path
from pymacruise.atlas_resources import get_bundled_atlas_descriptor, get_bundled_atlas_root
from pymacruise.toads_atlas import load_toads_atlas
from pymacruise.toads_atlas_type_membership import read_toads_atlas_class_types
from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors, load_lesion_atlas_topology
root = get_bundled_atlas_root()
assert root == Path(root).resolve() and root.is_dir()
default = get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt")
assert len(load_toads_atlas(default).structures) == 11
toads2012 = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt")
assert len(read_toads_atlas_class_types(toads2012)) == 13
load_lesion_atlas_topology(toads2012)
lesion = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt")
assert len(load_lesion_atlas_priors(lesion, ["T1_MPRAGE", "FLAIR"]).structures) == 13
assert not any(part in str(root) for part in ("toads-cruise-scm", "/dev/mipav"))
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = os.fspath(extract_to)
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=extract_to.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"unpacked wheel resolution failed:\n{result.stdout}\n{result.stderr}"


def _assert_isolated_wheel_install(wheel: Path, target: Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-deps",
            "--target",
            os.fspath(target),
            os.fspath(wheel),
        ],
        cwd=target.parent,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"isolated wheel install failed:\n{result.stdout}\n{result.stderr}"

    script = r"""
from pathlib import Path
import pymacruise
from pymacruise.atlas_resources import get_bundled_atlas_descriptor, get_bundled_atlas_root
from pymacruise.toads_atlas import load_toads_atlas
from pymacruise.toads_atlas_type_membership import read_toads_atlas_class_types
from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors, load_lesion_atlas_topology
target = Path(__import__('os').environ["PYMACRUISE_INSTALL_TARGET"]).resolve()
assert Path(pymacruise.__file__).resolve().is_relative_to(target)
root = get_bundled_atlas_root()
assert root == (target / "pymacruise" / "data" / "atlases").resolve()
assert len(load_toads_atlas(get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt")).structures) == 11
descriptor = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt")
assert len(read_toads_atlas_class_types(descriptor)) == 13
load_lesion_atlas_topology(descriptor)
lesion = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt")
assert len(load_lesion_atlas_priors(lesion, ["T1_MPRAGE", "FLAIR"]).structures) == 13
assert all(Path(path).exists() for path in (descriptor, lesion))
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = os.fspath(target)
    env["PYMACRUISE_INSTALL_TARGET"] = os.fspath(target)
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=target.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"installed wheel atlas resolution failed:\n{result.stdout}\n{result.stderr}"


def _assert_descriptor_references_are_bundled(root: Path, descriptor: Path) -> None:
    for line in descriptor.read_text(encoding="utf-8").splitlines():
        if line.startswith("File:"):
            relative = line.partition(":")[2].strip()
            assert relative.endswith(".nii.gz"), relative
            relative_path = Path(relative)
            assert not relative_path.is_absolute() and ".." not in relative_path.parts, relative
            assert (descriptor.parent / relative).is_file(), relative
