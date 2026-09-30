#!/usr/bin/env python3
"""Build and inspect PyMaCRUISE distributions and resolve atlases from a wheel.

Run from any directory with the project build requirements installed:

    python scripts/verify_atlas_distribution.py

All build products and the isolated wheel target are created under a temporary
directory. The checked-in atlas tree remains the sole repository copy.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ATLASES = ROOT / "data" / "atlases"
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
            members = {
                item.name.partition("/data/atlases/")[2]: item
                for item in tar.getmembers()
                if item.isfile() and "/data/atlases/" in item.name
            }
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
        names = set(wheel.namelist())
        prefix = prefix.rstrip("/") + "/"
        members = {
            name[len(prefix) :]: name
            for name in names
            if name.startswith(prefix)
        }
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
    with tempfile.TemporaryDirectory(prefix="pymacruise-atlas-dist-") as temp:
        temp_root = Path(temp)
        sdist_dir = temp_root / "sdist"
        wheelhouse = temp_root / "wheelhouse"
        sdist_dir.mkdir()
        wheelhouse.mkdir()
        run(
            sys.executable,
            "-m",
            "build",
            "--no-isolation",
            "--sdist",
            "--outdir",
            str(sdist_dir),
            str(ROOT),
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
assert root.resolve() == (Path(__import__("pymacruise").__file__).resolve().parent / "data" / "atlases").resolve()
print(json.dumps({"atlas_root": str(root), "default_classes": len(default.structures), "toads2012_classes": len(toads2012_types), "lesion_classes": len(priors.structures), "topology_shape": topology.image.shape}))
'''
        env = os.environ.copy()
        env["PYTHONPATH"] = str(target)
        env["PYTHONNOUSERSITE"] = "1"
        run(sys.executable, "-c", probe, cwd=temp_root, env=env)

        summary = {
            "atlas_source_files": expected_count,
            "wheel_atlas_files": wheel_count,
            "sdist_atlas_files": sdist_count,
            "source_atlas_files": expected_sdist_count,
            "wheel_bytes": wheel.stat().st_size,
            "sdist_bytes": sdist.stat().st_size,
            "wheel_sha256": sha256(wheel),
            "sdist_sha256": sha256(sdist),
        }
        print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
