from __future__ import annotations

import os
import hashlib
import json
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

    source_hashes = _source_descriptor_asset_hashes()
    script = r"""
from pathlib import Path
import hashlib
import json
import pymacruise
import nibabel as nib
from pymacruise.atlas_resources import get_bundled_atlas_descriptor, get_bundled_atlas_root
from pymacruise.toads_atlas import load_toads_atlas
from pymacruise.toads_atlas_type_membership import read_toads_atlas_class_types
from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors, load_lesion_atlas_topology
import os
target = Path(__import__('os').environ["PYMACRUISE_INSTALL_TARGET"]).resolve()
assert Path(pymacruise.__file__).resolve().is_relative_to(target)
root = get_bundled_atlas_root()
assert root == (target / "pymacruise" / "data" / "atlases").resolve()
default = load_toads_atlas(get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt"))
assert len(default.structures) == 11
assert nib.aff2axcodes(nib.load(str(default.topology.path)).affine) == ("R", "A", "S")
assert all(nib.aff2axcodes(nib.load(str(image.path)).affine) == ("R", "A", "S") for image in default.shapes.values())
descriptor = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt")
assert len(read_toads_atlas_class_types(descriptor)) == 13
toads2012_topology = load_lesion_atlas_topology(descriptor)
assert nib.aff2axcodes(toads2012_topology.image.affine) == ("R", "A", "S")
lesion = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt")
assert len(load_lesion_atlas_priors(lesion, ["T1_MPRAGE", "FLAIR"]).structures) == 13
lesion_topology = load_lesion_atlas_topology(lesion)
assert nib.aff2axcodes(lesion_topology.image.affine) == ("R", "A", "S")
expected_hashes = json.loads(os.environ["PYMACRUISE_ATLAS_HASHES"])
for relative, expected in expected_hashes.items():
    assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == expected, relative
"""
    env = os.environ.copy()
    env["PYTHONPATH"] = os.fspath(target)
    env["PYMACRUISE_INSTALL_TARGET"] = os.fspath(target)
    env["PYMACRUISE_ATLAS_HASHES"] = json.dumps(source_hashes)
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=target.parent,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"installed wheel atlas resolution failed:\n{result.stdout}\n{result.stderr}"

    _assert_installed_wheel_workflow_preflight(target)


def _assert_installed_wheel_workflow_preflight(target: Path) -> None:
    """Run the installed package's CLI against committed synthetic inputs."""
    run_root = target.parent / "preflight-run"
    config_dir = run_root / "config"
    input_dir = run_root
    fake_site = run_root / "fake-site"
    config_dir.mkdir(parents=True)
    fake_site.mkdir()

    fixture_root = PROJECT_ROOT / "tests" / "data" / "synthetic"
    config_source = fixture_root / "workflow_preflight" / "base.json"
    config_path = config_dir / "base.json"
    shutil.copy2(config_source, config_path)
    # These symlinks expose the committed fixture files at the config's
    # relative paths without copying or generating alternate image data.
    for filename in ("T1.nii.gz", "MaSeg.nii.gz", "brain_mask.nii.gz"):
        (input_dir / filename).symlink_to(fixture_root / filename)

    # Dependency preflight checks ANTsPy importability but does not call it.
    # This marker keeps the installed-wheel check independent of the optional
    # ANTsPy extra while still exercising the selected N4 dependency branch.
    (fake_site / "ants.py").write_text('"""preflight import marker"""\n', encoding="utf-8")
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [os.fspath(target.resolve()), os.fspath(fake_site.resolve())]
    )
    environment["PYTHONNOUSERSITE"] = "1"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pymacruise.cli",
            "workflow-preflight",
            "--config",
            os.fspath(config_path),
        ],
        cwd=run_root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"installed-wheel workflow-preflight failed:\n{result.stdout}\n{result.stderr}"
    )
    assert "Dry run only" in result.stdout
    assert "no N4 correction or segmentation workflow was executed" in result.stdout

    manifest_path = config_dir / "run-output" / "run.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    installed_root = (target / "pymacruise").resolve()
    atlas = manifest["preflight"]["atlas"]
    installed_atlas_root = (installed_root / "data" / "atlases").resolve()
    assert Path(atlas["descriptor"]).resolve().is_relative_to(installed_atlas_root)
    assert atlas["assets_checked"] > 1
    assert atlas["payloads_materialized"] is False
    assert manifest["dry_run"] is True
    assert manifest["preflight"]["status"] == "passed"
    assert manifest["preflight"]["workflow_executed"] is False
    assert manifest["stages"]["pipeline"] == "not_started"
    assert manifest["outputs"] == {}
    assert manifest["workflow_profile"] == "base"
    assert manifest["configuration_schema_version"] == 1
    assert manifest["source_revision"] is None
    assert manifest["dependencies"]["pymacruise"] != "not-installed"
    assert set(manifest["inputs"]) == {"t1", "maseg", "brain_mask"}
    expected_fixture_names = {
        "t1": "T1.nii.gz",
        "maseg": "MaSeg.nii.gz",
        "brain_mask": "brain_mask.nii.gz",
    }
    for role, record in manifest["inputs"].items():
        assert record["input_axis_codes"] == ["R", "A", "S"]
        assert record["output_axis_codes"] == ["R", "A", "S"]
        assert len(record["sha256"]) == 64
        assert len(record["affine_ras"]) == 4
        assert Path(record["path"]).resolve() == (
            fixture_root / expected_fixture_names[role]
        ).resolve()
    assert all(
        conversion["from"] == ["R", "A", "S"]
        and conversion["to"] == ["R", "A", "S"]
        and conversion["reoriented"] is False
        for conversion in manifest["preflight"]["orientation_conversions"].values()
    )
    assert sorted(path.name for path in manifest_path.parent.iterdir()) == ["run.json"]


def _source_descriptor_asset_hashes() -> dict[str, str]:
    atlas_root = PROJECT_ROOT / "data" / "atlases"
    hashes: dict[str, str] = {}
    for descriptor, _ in ATLAS_CASES.items():
        path = atlas_root / descriptor
        hashes[path.relative_to(atlas_root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.startswith("File:"):
                continue
            relative = Path(line.partition(":")[2].strip())
            assert not relative.is_absolute() and ".." not in relative.parts
            asset = (path.parent / relative).resolve()
            assert asset.is_relative_to(atlas_root.resolve()) and asset.is_file()
            key = asset.relative_to(atlas_root.resolve()).as_posix()
            hashes[key] = hashlib.sha256(asset.read_bytes()).hexdigest()
    return hashes


def _assert_descriptor_references_are_bundled(root: Path, descriptor: Path) -> None:
    for line in descriptor.read_text(encoding="utf-8").splitlines():
        if line.startswith("File:"):
            relative = line.partition(":")[2].strip()
            assert relative.endswith(".nii.gz"), relative
            relative_path = Path(relative)
            assert not relative_path.is_absolute() and ".." not in relative_path.parts, relative
            assert (descriptor.parent / relative).is_file(), relative
