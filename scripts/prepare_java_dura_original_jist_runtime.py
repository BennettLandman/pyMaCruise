#!/usr/bin/env python3
"""Stage selected Vandy classes beside original JIST numerical classes."""

from __future__ import annotations

import shutil
import sys
import zipfile
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: prepare_java_dura_original_jist_runtime.py JAR CLASSES_DIR RUNTIME_DIR")
    jar_path, classes_path, runtime_path = map(Path, sys.argv[1:])
    runtime_path.mkdir(parents=True, exist_ok=True)
    filtered_jar = runtime_path.parent / "jist-runtime-without-selected-conflicts.jar"
    omitted = (
        "edu/jhu/ece/iacl/algorithms/tgdm/",
        "edu/jhu/ece/iacl/algorithms/gvf/FastMarchingGradient",
        "edu/jhu/ece/iacl/algorithms/graphics/isosurf/",
    )
    with zipfile.ZipFile(jar_path) as source, zipfile.ZipFile(
        filtered_jar, "w", compression=zipfile.ZIP_DEFLATED
    ) as target:
        seen: set[str] = set()
        for info in source.infolist():
            if info.filename.startswith(omitted) or info.filename in seen:
                continue
            target.writestr(info, source.read(info.filename))
            seen.add(info.filename)

    keep = (
        Path("DuraTGDMFullEvolutionReferenceHarness.class"),
        Path("DuraTGDMFullEvolutionReferenceHarness$Probe.class"),
        Path("edu/jhu/ece/iacl/algorithms/tgdm"),
        Path("edu/jhu/ece/iacl/algorithms/gvf"),
        Path("edu/jhu/ece/iacl/algorithms/topology"),
        Path("edu/jhu/ece/iacl/algorithms/graphics"),
        Path("edu/jhu/ece/iacl/jist/pipeline/AbstractCalculation.class"),
        Path("edu/jhu/ece/iacl/jist/structures/geom/EmbeddedSurface.class"),
        Path("gov/nih/mipav/view/Preferences.class"),
    )
    for relative in keep:
        source = classes_path / relative
        destination = runtime_path / relative
        if source.is_dir():
            shutil.copytree(source, destination, dirs_exist_ok=True)
        elif source.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        else:
            raise FileNotFoundError(source)
    print(filtered_jar)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
