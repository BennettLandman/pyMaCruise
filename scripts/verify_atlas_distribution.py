#!/usr/bin/env python3
"""Build and inspect PyMaCRUISE distributions and resolve atlases from a wheel.

Run from any directory with the project build requirements installed:

    python scripts/verify_atlas_distribution.py

All build products and the isolated wheel target are preserved under
``tmp/build`` in this project. The checked-in atlas tree remains the sole
repository copy.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ATLASES = ROOT / "data" / "atlases"
BUILD_ROOT = ROOT / "tmp" / "build"
EXPECTED_DESCRIPTORS = (
    "cruise-atlas-10obj-toads2008.txt",
    "lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt",
    "lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(*command: str, cwd: Path = ROOT, env: dict[str, str] | None = None) -> None:
    subprocess.run(command, cwd=cwd, env=env, check=True)


def check_archive(archive: Path, prefix: str, *, sdist: bool) -> tuple[int, int]:
    source_files = {
        path.relative_to(SOURCE_ATLASES).as_posix(): path
        for path in SOURCE_ATLASES.rglob("*")
        if path.is_file()
    }
    if sdist:
        with tarfile.open(archive, "r:gz") as tar:
            atlas_members = [
                (item.name.partition("/data/atlases/")[2], item)
                for item in tar.getmembers()
                if item.isfile() and "/data/atlases/" in item.name
            ]
            members = {relative: item for relative, item in atlas_members}
            if len(atlas_members) != len(source_files):
                raise AssertionError(f"{archive.name} contains duplicate or extra atlas entries")
            if set(members) != set(source_files):
                missing = sorted(set(source_files) - set(members))
                extra = sorted(set(members) - set(source_files))
                raise AssertionError(
                    f"{archive.name} atlas file set differs from source; "
                    f"missing={missing}, extra={extra}"
                )
            for relative, source in source_files.items():
                member = members.get(relative)
                if member is None:
                    raise AssertionError(f"{archive.name} omits atlas file {relative}")
                extracted = tar.extractfile(member)
                if extracted is None:
                    raise AssertionError(f"Cannot read {member.name} from {archive.name}")
                digest = hashlib.sha256(extracted.read()).hexdigest()
                if digest != sha256(source):
                    raise AssertionError(f"Hash mismatch in sdist for atlas file {relative}")
            for descriptor in EXPECTED_DESCRIPTORS:
                if descriptor not in members:
                    raise AssertionError(f"{archive.name} omits {descriptor}")
            return len(members), len(source_files)

    with zipfile.ZipFile(archive) as wheel:
        archive_names = wheel.namelist()
        names = set(archive_names)
        prefix = prefix.rstrip("/") + "/"
        atlas_names = [
            (name[len(prefix) :], name)
            for name in archive_names
            if name.startswith(prefix)
        ]
        members = {relative: name for relative, name in atlas_names}
        if len(atlas_names) != len(source_files):
            raise AssertionError(f"{archive.name} contains duplicate or extra atlas entries")
        if set(members) != set(source_files):
            missing = sorted(set(source_files) - set(members))
            extra = sorted(set(members) - set(source_files))
            raise AssertionError(
                f"{archive.name} atlas file set differs from source; "
                f"missing={missing}, extra={extra}"
            )
        for relative, source in source_files.items():
            name = members.get(relative)
            if name is None:
                raise AssertionError(f"{archive.name} omits atlas file {relative}")
            digest = hashlib.sha256(wheel.read(name)).hexdigest()
            if digest != sha256(source):
                raise AssertionError(f"Hash mismatch in wheel for atlas file {relative}")
        for descriptor in EXPECTED_DESCRIPTORS:
            if descriptor not in members:
                raise AssertionError(f"{archive.name} omits {descriptor}")
        return len(members), len(source_files)


def main() -> int:
    if not SOURCE_ATLASES.is_dir():
        raise SystemExit(f"Missing source atlas directory: {SOURCE_ATLASES}")
    source_nifti_count = len(list(SOURCE_ATLASES.rglob("*.nii.gz")))
    # Keep every build, pip target/cache, subprocess temp file, and inspection
    # artifact inside the project-owned scratch tree. Preserve each run so its
    # exact archives and probe result can be inspected after validation.
    BUILD_ROOT.mkdir(parents=True, exist_ok=True)
    temp_root = Path(tempfile.mkdtemp(prefix=f"atlas-distribution-{os.getpid()}-", dir=BUILD_ROOT))
    project_copy = temp_root / "source"
    project_copy.mkdir()
    for name in ("pyproject.toml", "setup.py", "MANIFEST.in", "readme.md"):
        shutil.copy2(ROOT / name, project_copy / name)
    shutil.copytree(ROOT / "src", project_copy / "src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copytree(SOURCE_ATLASES, project_copy / "data" / "atlases")

    sdist_dir = temp_root / "sdist"
    wheelhouse = temp_root / "wheelhouse"
    sdist_dir.mkdir()
    wheelhouse.mkdir()
    subprocess_tmp = temp_root / "tmp"
    pip_cache = temp_root / "pip-cache"
    subprocess_tmp.mkdir()
    pip_cache.mkdir()
    env = os.environ.copy()
    env["TMPDIR"] = str(subprocess_tmp)
    env["TEMP"] = str(subprocess_tmp)
    env["TMP"] = str(subprocess_tmp)
    env["PIP_CACHE_DIR"] = str(pip_cache)
    env["PIP_NO_CACHE_DIR"] = "1"
    run(
        sys.executable,
        "-m",
        "build",
        "--no-isolation",
        "--sdist",
        "--outdir",
        str(sdist_dir),
        str(project_copy),
        cwd=project_copy,
        env=env,
    )
    sdist = next(sdist_dir.glob("pymacruise-*.tar.gz"))
    run(
        sys.executable,
        "-m",
        "pip",
        "wheel",
        "--no-deps",
        "--no-build-isolation",
        "--wheel-dir",
        str(wheelhouse),
        str(sdist),
        cwd=temp_root,
        env=env,
    )
    wheel = next(wheelhouse.glob("pymacruise-*.whl"))

    wheel_prefix = "pymacruise/data/atlases"
    wheel_count, expected_count = check_archive(wheel, wheel_prefix, sdist=False)
    sdist_count, expected_sdist_count = check_archive(sdist, "", sdist=True)

    target = temp_root / "isolated-target"
    run(
        sys.executable,
        "-m",
        "pip",
        "install",
        "--no-deps",
        "--target",
        str(target),
        str(wheel),
        cwd=temp_root,
        env=env,
    )
    probe = r'''import json
from pathlib import Path
import nibabel as nib
from pymacruise.atlas_resources import get_bundled_atlas_descriptor, get_bundled_atlas_root
from pymacruise.toads_atlas import load_toads_atlas
from pymacruise.toads_atlas_type_membership import read_toads_atlas_class_types
from pymacruise.lesion_atlas_priors import load_lesion_atlas_priors, load_lesion_atlas_topology
root = get_bundled_atlas_root()
default = load_toads_atlas(get_bundled_atlas_descriptor("cruise-atlas-10obj-toads2008.txt"))
toads2012_path = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt")
toads2012_types = read_toads_atlas_class_types(toads2012_path)
toads2012_topology = load_lesion_atlas_topology(toads2012_path)
lesion_path = get_bundled_atlas_descriptor("lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt")
priors = load_lesion_atlas_priors(lesion_path, ["T1_MPRAGE", "FLAIR"])
topology = load_lesion_atlas_topology(lesion_path)
assert len(default.structures) == 11 and len(toads2012_types) == 13
assert len(priors.structures) == 13 and nib.aff2axcodes(topology.image.affine) == ("R", "A", "S")
assert default.topology.path.is_file() and all(v.path.is_file() for v in default.shapes.values())
assert toads2012_topology.path.is_file()
all_nifti = sorted(root.rglob("*.nii.gz"))
assert not list(root.rglob("*.nii")), "uncompressed .nii atlas image found"
nonras = [(str(path.relative_to(root)), nib.aff2axcodes(nib.load(path).affine))
          for path in all_nifti if nib.aff2axcodes(nib.load(path).affine) != ("R", "A", "S")]
assert not nonras, f"atlas NIfTI images are not RAS: {nonras}"
for descriptor in (root / "cruise-atlas-10obj-toads2008.txt",
                   root / "lesionToads-atlas-2012/cruise-atlas-12obj-toads2012.txt",
                   root / "lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt"):
    for line in descriptor.read_text(encoding="utf-8").splitlines():
        if line.startswith("File:"):
            relative = line.partition(":")[2].strip()
            assert relative.endswith(".nii.gz"), (descriptor, relative)
            assert (descriptor.parent / relative).is_file(), (descriptor, relative)
assert root.resolve() == (Path(__import__("pymacruise").__file__).resolve().parent / "data" / "atlases").resolve()
from pymacruise.asset_inventory import verify_bundled_assets
receipt = verify_bundled_assets()
assert receipt["runtime_tables_checked"] == 2
assert receipt["ibsr11objV3_referenced_payloads_checked"] == 13
print(json.dumps({"installed_asset_verification": receipt}))
print(json.dumps({"atlas_root": str(root), "default_classes": len(default.structures), "toads2012_classes": len(toads2012_types), "lesion_classes": len(priors.structures), "topology_shape": topology.image.shape, "ras_nifti_count": len(all_nifti)}))
'''
    probe_env = env.copy()
    probe_env["PYTHONPATH"] = str(target)
    probe_env["PYTHONNOUSERSITE"] = "1"
    probe_env["PYTHONDONTWRITEBYTECODE"] = "1"
    run(sys.executable, "-c", probe, cwd=temp_root, env=probe_env)

    summary = {
        "build_directory": str(temp_root),
        "atlas_source_files": expected_count,
        "wheel_atlas_files": wheel_count,
        "sdist_atlas_files": sdist_count,
        "source_atlas_files": expected_sdist_count,
        "ras_nifti_files_checked": source_nifti_count,
        "wheel_bytes": wheel.stat().st_size,
        "sdist_bytes": sdist.stat().st_size,
        "wheel_sha256": sha256(wheel),
        "sdist_sha256": sha256(sdist),
    }
    (temp_root / "verification.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
