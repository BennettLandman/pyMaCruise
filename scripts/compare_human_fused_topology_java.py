#!/usr/bin/env python3
"""Compare the saved human fused-base topology field with selected Java source.

This is a held-fixed scalar-stage replay. It uses the selected
``AlgorithmTopologyCorrection`` Java source and only substitutes small
in-memory MIPAV framework shims, matching the existing synthetic source
harness. The human volumes, raw transfer files, compiled classes, and results
stay under the ignored ``build/`` directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tempfile
import time

import nibabel as nib
import numpy as np
from scipy import ndimage


REPO = Path(__file__).resolve().parents[1]
CASE = REPO / "tests/java_reference/topology_correction"
DEFAULT_CORE = REPO / "build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11"
DEFAULT_OUT = REPO / "build/invivo-base-2026-10-10-ras-grid/source-selected-fused-core-v3-t1-grid-2026-10-11/topology-java-current-human"
HARNESS = r'''import edu.jhmi.rad.medic.algorithms.AlgorithmTopologyCorrection;
import gov.nih.mipav.model.structures.ModelImage;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.BitSet;
public final class HumanTopologyCorrectionReplay {
  static float[] read(Path p, int n) throws Exception {
    byte[] b=Files.readAllBytes(p);
    if (b.length != Math.multiplyExact(n,4)) throw new IllegalArgumentException("invalid input length: "+p);
    ByteBuffer q=ByteBuffer.wrap(b).order(ByteOrder.LITTLE_ENDIAN);
    float[] a=new float[n]; for(int i=0;i<n;i++) a[i]=q.getFloat(); return a;
  }
  public static void main(String[] a) throws Exception {
    if(a.length!=5) throw new IllegalArgumentException("usage: INPUT_DIR OUTPUT NX NY NZ");
    Path in=Path.of(a[0]), out=Path.of(a[1]);
    int nx=Integer.parseInt(a[2]), ny=Integer.parseInt(a[3]), nz=Integer.parseInt(a[4]);
    int n=Math.multiplyExact(Math.multiplyExact(nx,ny),nz);
    float[] x=read(in.resolve("wmfill.f32le"),n), p=read(in.resolve("paint.f32le"),n);
    ModelImage src=new ModelImage(new int[]{nx,ny,nz},"human_wmfill");
    ModelImage dst=new ModelImage(new int[]{nx,ny,nz},"human_wmfill_topology_corrected");
    src.setData(x); BitSet paint=new BitSet(n); for(int i=0;i<n;i++) if(p[i]>0.5f) paint.set(i);
    AlgorithmTopologyCorrection op=new AlgorithmTopologyCorrection(dst,src,"scalar_image",1.0e-4f,3,6,18,
      "paint_mask",paint,0.0f,1.0f,true,"object->background");
    op.runAlgorithm(); if(!op.isCompleted()) throw new IllegalStateException("source algorithm incomplete");
    float[] y=dst.getData(); ByteBuffer b=ByteBuffer.allocate(Math.multiplyExact(n,4)).order(ByteOrder.LITTLE_ENDIAN);
    for(float v:y)b.putFloat(v); Files.write(out,b.array());
    System.out.println("completed=true dimensions="+nx+"x"+ny+"x"+nz+" paint_count="+paint.cardinality());
  }
}
'''
SOURCES = [
    "stubs/gov/nih/mipav/model/algorithms/AlgorithmBase.java",
    "stubs/gov/nih/mipav/model/structures/ModelImage.java",
    "stubs/gov/nih/mipav/view/ViewUserInterface.java",
    "stubs/gov/nih/mipav/view/ViewJProgressBar.java",
    "stubs/gov/nih/mipav/view/Preferences.java",
    "stubs/edu/jhmi/rad/medic/utilities/MedicUtilPublic.java",
    "stubs/edu/jhmi/rad/medic/libraries/ObjectProcessing.java",
    "stubs/WildMagic/LibFoundation/Mathematics/Vector3f.java",
    "upstream/edu/jhmi/rad/medic/structures/BinaryTree.java",
    "upstream/edu/jhmi/rad/medic/libraries/NeighborConnectivity.java",
    "upstream/edu/jhmi/rad/medic/libraries/TopologyPropagation.java",
    "upstream/edu/jhmi/rad/medic/algorithms/AlgorithmTopologyCorrection.java",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_f32_fortran(image_path: Path, output_path: Path) -> tuple[nib.Nifti1Image, dict]:
    image = nib.load(str(image_path))
    values = np.asarray(image.dataobj, dtype=np.float32)
    if values.ndim != 3:
        raise ValueError(f"expected 3D input: {image_path} has {values.shape}")
    if not np.isfinite(values).all():
        raise ValueError(f"non-finite voxels in {image_path}")
    # Crop products are already explicitly RAS in the selected workflow.
    orientation = nib.aff2axcodes(image.affine)
    if orientation != ("R", "A", "S"):
        raise ValueError(f"expected RAS image at {image_path}; got {orientation}")
    np.asarray(values, dtype="<f4").ravel(order="F").tofile(output_path)
    return image, {"shape_xyz": list(image.shape), "affine": image.affine.tolist(), "orientation": orientation}


def run_monitored(command: list[str], *, cwd: Path, timeout_s: int, rss_limit_gib: float) -> dict:
    started = time.monotonic()
    proc = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        stdout, _ = proc.communicate(timeout=timeout_s)
        if proc.returncode:
            raise RuntimeError(f"Java process failed ({proc.returncode}):\n{stdout}")
        child_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        # macOS reports bytes; Linux and other supported POSIX hosts report KiB.
        peak_rss = int(child_rss if sys.platform == "darwin" else child_rss * 1024)
        if peak_rss > rss_limit_gib * (1024**3):
            raise RuntimeError(f"Java max RSS {peak_rss / (1024**3):.3f} GiB exceeded {rss_limit_gib:.2f} GiB")
    except subprocess.TimeoutExpired:
        proc.kill()
        stdout, _ = proc.communicate()
        raise TimeoutError(f"Java replay exceeded {timeout_s} seconds:\n{stdout}")
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
    return {"elapsed_seconds": time.monotonic() - started, "peak_process_rss_bytes": peak_rss,
            "stdout": stdout}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-dir", type=Path, default=DEFAULT_CORE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--java-heap-mib", type=int, default=2304)
    parser.add_argument("--rss-limit-gib", type=float, default=3.8)
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    args = parser.parse_args()
    core = args.core_dir.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    source_path = core / "wmfill_cropped.nii.gz"
    paint_path = core / "wmmask_cropped.nii.gz"
    python_path = core / "wmfill_topology_corrected.nii.gz"
    input_dir = out / "input"
    input_dir.mkdir(exist_ok=True)
    source_img, source_grid = write_f32_fortran(source_path, input_dir / "wmfill.f32le")
    paint_img, paint_grid = write_f32_fortran(paint_path, input_dir / "paint.f32le")
    if source_img.shape != paint_img.shape or not np.allclose(source_img.affine, paint_img.affine, atol=1e-5, rtol=0):
        raise ValueError("source scalar and paint mask are not on the same grid")
    dimensions = source_img.shape
    class_dir = out / "classes"
    class_dir.mkdir(exist_ok=True)
    harness_path = out / "HumanTopologyCorrectionReplay.java"
    harness_path.write_text(HARNESS)
    source_files = [CASE / rel for rel in SOURCES] + [harness_path]
    javac = shutil.which("javac")
    java = shutil.which("java")
    if not javac or not java:
        raise RuntimeError("JDK javac and java are required")
    compile_cmd = [javac, "-d", str(class_dir), *map(str, source_files)]
    compile_started = time.monotonic()
    compiled = subprocess.run(compile_cmd, check=True, capture_output=True, text=True)
    compile_elapsed = time.monotonic() - compile_started
    java_out = out / "java-output.f32le"
    java_cmd = [java, f"-Xmx{args.java_heap_mib}m", "-cp", str(class_dir), "HumanTopologyCorrectionReplay",
                str(input_dir), str(java_out), *(str(int(n)) for n in dimensions)]
    java_stats = run_monitored(java_cmd, cwd=REPO, timeout_s=args.timeout_seconds, rss_limit_gib=args.rss_limit_gib)

    python_img = nib.load(str(python_path))
    java_values = np.fromfile(java_out, dtype="<f4")
    py_values = np.asarray(python_img.dataobj, dtype=np.float32)
    if tuple(dimensions) != py_values.shape or py_values.shape != java_values.reshape(dimensions, order="F").shape:
        raise ValueError("Java and Python output grids do not match")
    java_values = java_values.reshape(dimensions, order="F")
    delta = np.abs(java_values - py_values)
    thresholds = (0.0, 63.75, 127.5, 191.25, 254.0)
    level_comparisons = {}
    for threshold in thresholds:
        java_mask = java_values > threshold
        py_mask = py_values > threshold
        level_comparisons[str(threshold)] = {
            "java_foreground_voxels": int(java_mask.sum()),
            "python_foreground_voxels": int(py_mask.sum()),
            "different_voxels": int(np.count_nonzero(java_mask != py_mask)),
            "dice": float(2 * np.count_nonzero(java_mask & py_mask) / (java_mask.sum() + py_mask.sum()))
            if java_mask.sum() + py_mask.sum() else 1.0,
        }
    foreground = np.asarray(py_values > 0.0)
    structure6 = ndimage.generate_binary_structure(3, 1)
    component_count = int(ndimage.label(foreground, structure=structure6)[1])
    report = {
        "result": "PASS" if not np.count_nonzero(delta) else "DIFFERENCE",
        "comparison": "saved Python wmfill_topology_corrected vs a held-fixed replay of selected Java AlgorithmTopologyCorrection",
        "harness_semantics": "The harness invokes the actual selected AlgorithmTopologyCorrection.java, TopologyPropagation.java, NeighborConnectivity.java and BinaryTree.java sources copied into pymacruise/tests/java_reference. Only MIPAV UI/image framework classes are shims. This is not a compiled full MIPAV or full JIST workflow execution.",
        "input_shape_xyz": list(dimensions), "orientation": "RAS",
        "grid": source_grid,
        "inputs": {"scalar": {"path": str(source_path), "sha256": sha256(source_path)},
                   "paint_mask": {"path": str(paint_path), "sha256": sha256(paint_path)},
                   "python_output": {"path": str(python_path), "sha256": sha256(python_path)}},
        "java_output_sha256": sha256(java_out),
        "java_source_hashes": {str(p.relative_to(REPO)): sha256(p) for p in source_files if p.is_file()},
        "settings": {"type": "scalar_image", "regularization": 0.0001, "dimension": 3,
                     "object_connectivity": 6, "background_connectivity": 18,
                     "paint_mode": "paint_mask", "paint_threshold_strict": 0.5,
                     "lowest_level": 0.0, "highest_level": 1.0,
                     "minmax": True, "propagation": "object->background"},
        "scalar_difference": {"voxel_count": int(delta.size),
                              "different_voxels_exact": int(np.count_nonzero(delta)),
                              "max_absolute": float(delta.max(initial=0)),
                              "mean_absolute": float(delta.mean()),
                              "python_range": [float(py_values.min()), float(py_values.max())],
                              "java_range": [float(java_values.min()), float(java_values.max())]},
        "superlevel_masks": level_comparisons,
        "python_positive_6_connected_components": component_count,
        "runtime": {"javac_seconds": compile_elapsed, **java_stats,
                    "java_heap_limit_mib": args.java_heap_mib,
                    "rss_limit_gib": args.rss_limit_gib,
                    "rss_limit_exceeded": False},
        "reproduction": [
            "cd /Users/landmaba/dev/mipav/pymacruise",
            "PYTHONPATH=src .venv/bin/python scripts/compare_human_fused_topology_java.py",
        ],
    }
    report_path = out / "comparison.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["result"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
