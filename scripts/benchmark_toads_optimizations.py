"""Benchmark strict TOADS kernels in a complete selected workflow.

Default: committed folded-atlas phantom through natural convergence (cap 20).
--source-run RUN.json replays one human TOADS stage with its recorded command,
redirecting every output under --output-dir. Inputs and source settings stay
fixed. --reference selects the retained Python implementations of only the
three optimized components. Reports include per-component timings and hashes;
run each backend in a fresh process with identical thread settings.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import runpy
import sys
import time

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def compare_nifti_outputs(actual_dir, reference_dir):
    """Audit all emitted tissue products, including dtype and physical grid."""
    comparisons = {}
    for path in sorted(Path(actual_dir).glob('*.nii.gz')):
        reference = Path(reference_dir) / path.name
        actual_image = nib.load(path)
        reference_image = nib.load(reference)
        actual = np.asarray(actual_image.dataobj)
        expected = np.asarray(reference_image.dataobj)
        same_shape = actual.shape == expected.shape
        comparisons[path.name] = dict(
            same_shape=same_shape,
            same_dtype=actual.dtype == expected.dtype,
            identical_grid=np.array_equal(actual_image.affine, reference_image.affine),
            differing_values=int(np.count_nonzero(actual != expected)) if same_shape else None,
            identical_value_bytes=(same_shape and actual.dtype == expected.dtype
                                   and actual.tobytes(order='F') == expected.tobytes(order='F')),
            identical_file=hashlib.sha256(path.read_bytes()).digest() == hashlib.sha256(reference.read_bytes()).digest(),
        )
    return comparisons


def _instrument(backend, timings):
    # Import the full CLI before replacing already-bound component references.
    import pymacruise.cli
    from pymacruise.toads_inhomogeneity_exact_fit import fit_toads_exact_chebyshev_field
    from pymacruise.toads_thinning import toads_propagate_growing
    from pymacruise.toads_relations import (
        toads_closest_labels, toads_compute_relations, toads_update_boundary_distances,
    )
    targets = (fit_toads_exact_chebyshev_field, toads_propagate_growing,
               toads_closest_labels, toads_compute_relations, toads_update_boundary_distances)

    def decorate(original):
        @functools.wraps(original)
        def wrapped(*args, **kwargs):
            kwargs['backend'] = backend
            started = time.perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                name = original.__name__
                if name == 'fit_toads_exact_chebyshev_field':
                    name += f".degree{kwargs.get('degree', 1)}"
                timings.setdefault(name, []).append(time.perf_counter() - started)
        return wrapped

    replacements = {id(target): decorate(target) for target in targets}
    for name, module in list(sys.modules.items()):
        if not name.startswith('pymacruise.') or module is None:
            continue
        for key, value in list(vars(module).items()):
            replacement = replacements.get(id(value))
            if replacement is not None:
                setattr(module, key, replacement)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--source-run', type=Path)
    parser.add_argument('--stage', choices=('toads-first', 'toads-second'), default='toads-first')
    parser.add_argument('--reference', action='store_true')
    parser.add_argument('--compare', type=Path, help='classification from a matched comparison run')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    classification = output / 'classification.nii.gz'
    report_path = output / 'performance.json'
    if classification.exists() or report_path.exists():
        parser.error('refusing to overwrite an existing benchmark')
    backend = 'python' if args.reference else 'numba'
    timings = {}
    _instrument(backend, timings)
    if args.source_run:
        manifest = json.loads(args.source_run.read_text())
        command = list(manifest['stages'][args.stage]['command'])
        command[command.index('--output') + 1] = str(classification)
        if '--dura-output-dir' in command:
            command[command.index('--dura-output-dir') + 1] = str(output / 'dura-inputs')
        from pymacruise.cli import main as execute
        started = time.perf_counter()
        exit_code = execute(command)
        elapsed = time.perf_counter() - started
        baseline_output = Path(manifest['stages'][args.stage]['command'][
            manifest['stages'][args.stage]['command'].index('--output') + 1])
        reference_path = args.compare or baseline_output
    else:
        command = ['trace_synthetic_toads_loop.py', '--fixture', str(ROOT / 'tests/data/synthetic/toads_downsampled_vandy_atlas'),
                   '--atlas', str(ROOT / 'tests/data/synthetic/toads_downsampled_vandy_atlas/atlas.txt'),
                   '--maximum-iterations', '20', '--output', str(classification)]
        namespace = runpy.run_path(str(ROOT / 'scripts/trace_synthetic_toads_loop.py'))
        previous_argv = sys.argv
        sys.argv = command
        try:
            started = time.perf_counter()
            exit_code = namespace['main']()
            elapsed = time.perf_counter() - started
        finally:
            sys.argv = previous_argv
        # The Java fixture uses the full atlas and padding=10; this compact
        # benchmark uses its own atlas and padding=2. Compare matched runs.
        reference_path = args.compare
    labels = np.asarray(nib.load(classification).dataobj)
    changed = None
    if reference_path is not None:
        reference_labels = np.asarray(nib.load(reference_path).dataobj)
        if reference_labels.shape != labels.shape:
            raise ValueError('comparison classification must have the same shape')
        changed = int(np.count_nonzero(labels != reference_labels))
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        peak_rss *= 1024
    report = dict(backend=backend, source_run=str(args.source_run) if args.source_run else None,
                  command=command, elapsed_seconds=elapsed, peak_process_rss_gb=peak_rss / 1e9,
                  exit_code=exit_code, includes_jit_and_output_io=True,
                  comparison_reference=str(reference_path) if reference_path is not None else None,
                  differing_reference_label_voxels=changed,
                  classification_sha256=hashlib.sha256(classification.read_bytes()).hexdigest(),
                  components={name: dict(calls=len(values), seconds=sum(values), individual_seconds=values)
                              for name, values in timings.items()})
    report['environment'] = dict(
        python=platform.python_version(), platform=platform.platform(),
        thread_settings={key: os.environ.get(key) for key in (
            'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS')},
    )
    report['implementation_sha256'] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (ROOT / 'src/pymacruise').glob('toads_*compiled.py')
    }
    if args.source_run and '--dura-output-dir' in command:
        original_command = manifest['stages'][args.stage]['command']
        reference_dir = original_command[original_command.index('--dura-output-dir') + 1]
        report['tissue_output_comparisons'] = compare_nifti_outputs(output / 'dura-inputs', reference_dir)
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
