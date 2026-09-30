"""Setuptools build hooks for repository-supplied atlas assets.

Atlas files stay in one checked-in source tree at ``data/atlases``. The build
copies them into the installed package resource directory so the wheel does
not depend on the repository layout at runtime.
"""

from __future__ import annotations

from pathlib import Path
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


ROOT = Path(__file__).resolve().parent
ATLAS_SOURCE = ROOT / "data" / "atlases"
ATLAS_RESOURCE = Path("pymacruise") / "data" / "atlases"


class build_py(_build_py):
    """Add root-level atlases to the package at build time, without source copies."""

    def run(self) -> None:
        super().run()
        if not ATLAS_SOURCE.is_dir():
            raise FileNotFoundError(f"Required bundled atlas tree is missing: {ATLAS_SOURCE}")
        destination = Path(self.build_lib) / ATLAS_RESOURCE
        shutil.copytree(ATLAS_SOURCE, destination, dirs_exist_ok=True)

    def get_outputs(self, include_bytecode: int = 1) -> list[str]:
        outputs = list(super().get_outputs(include_bytecode))
        if ATLAS_SOURCE.is_dir():
            destination = Path(self.build_lib) / ATLAS_RESOURCE
            outputs.extend(
                str(destination / source.relative_to(ATLAS_SOURCE))
                for source in ATLAS_SOURCE.rglob("*")
                if source.is_file()
            )
        return outputs


setup(cmdclass={"build_py": build_py})
