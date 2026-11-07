# PyMaCRUISE supervised BASE testing

> Superseded as the primary tester handoff: the requested test is a complete external-versus-Python end-to-end comparison. See [end-to-end protocol](end-to-end-external.md). This document provides supplemental Python installation, execution, and integrity checks; it does not satisfy the external comparison requirement.

Prepared October 16, 2026. Read-only access to https://github.com/BennettLandman/pyMaCruise is sufficient. You need permission to create files on your own computer, but no GitHub write access, fork, commit, or push is needed.

## What is ready

The BASE workflow is ready for an initial technical-user evaluation: installation, execution, restart, and visual review of scientific outputs. Whole-workflow equivalence to original MaCRUISE remains unverified. Completing this protocol does not establish scientific or clinical validation.

PLUS is now available for exploratory evaluation through [the full PLUS walkthrough](full-plus-walkthrough.md). AtlasEM, pre-LesionTOADS Dura, selected native/Python atlas preparation and final refinement are implemented. Full original-versus-Python PLUS agreement remains unverified.

Clone the published revision containing this guide and record its full Git commit. Keep the checkout unchanged during testing and resume. The earlier synthetic check described below used commit `16993385fde1ef0ed9a7b8e09d667ace740ac18d`; it is historical execution evidence, not a claim that every later revision has been tested.

The small synthetic command below was checked from an export of that exact revision: all eight stages completed; 68 stage-output hashes verified; 24 final-output records were present; all image outputs were finite and RAS; resume verified and skipped all eight stages. This check reused an existing dependency environment. The fresh installation below is part of the tester's evaluation. The full MRI run was not repeated when preparing this guide.

## 1. Prepare and install

Use macOS for the first evaluation, with Git and Python 3.13 installed. macOS is the platform with existing execution evidence. Other platforms can be tested separately; their installation and execution are not yet verified. Use a machine with at least 16 GB RAM; 32 GB is preferable. Reserve roughly 10 GB free disk space as a working budget. These are suggested testing resources, not measured minimum requirements. First execution can pause while scientific kernels compile.

Use a fresh folder. The commands below use Bash; start `bash` first if your terminal normally uses another shell. Authenticate to GitHub using your existing account when prompted; do not put credentials into the URL.

```bash
bash
set -o pipefail
git clone https://github.com/BennettLandman/pyMaCruise.git pymacruise
cd pymacruise
git checkout --detach HEAD

python3.13 -m venv .venv
source .venv/bin/activate
unset PYTHONPATH
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p build/tester-feedback

python -m pip install --upgrade pip
python -m pip install '.[ants]' 2>&1 | tee build/tester-feedback/install.log
python -m pip check
pymacruise --version

git rev-parse HEAD > build/tester-feedback/revision.txt
python --version > build/tester-feedback/python-version.txt
python -m pip freeze > build/tester-feedback/packages.txt
uname -a > build/tester-feedback/system.txt
python -c 'import numpy, scipy, nibabel, numba, ants; print("Scientific dependencies imported successfully")'
```

Stop and return the installation log if installation, `pip check`, or an import fails. MATLAB, MIPAV/JIST, FreeSurfer, and a GPU are not required to run this Python test. Atlas files and the example inputs are regular repository files; no separate atlas download is needed. There is no integrated graphical interface.

## 2. Run the small synthetic example and restart it

Run from the repository root, in the same activated shell. Define this function once, then invoke it as shown. The function keeps every option identical for resume.

```bash
smoke_base() {
  pymacruise base-run \
    --t1 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
    --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
    --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
    --maseg tests/data/synthetic/base_workflow/maseg.nii.gz \
    --output-dir build/tester-smoke \
    --registration-levels 1 --registration-iterations 1 \
    --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
    --toads-iterations 2 --dura-iterations 1 --crop-edge-width 2 \
    --ace-diffusion-iterations 1 --tgdm-iterations 1 "$@"
}

smoke_base 2>&1 | tee build/tester-feedback/smoke.log
# Continue only if the command above succeeds.
smoke_base --resume 2>&1 | tee build/tester-feedback/smoke-resume.log
```

Expected: ten `Stage check passed` messages on the first run; ten `Resume: verified ...; skipping` messages on the second, including native export. `build/tester-smoke/run.json` should say `status: complete`. Both TOADS passes intentionally stop at two iterations; `converged: false` is expected for this small example. Do not judge anatomy from this fixture.

Do not reuse the same output folder for a new non-resume run. An existing-run error is intentional; use `--resume` with the same options, or choose a fresh output directory.

## 3. Run the bundled human MRI BASE example

Use the three files under `source-reference-run/inputs/base-root/` together. Do not mix them with `inputs/100260/`, which is a different scan. Before starting, verify the packet's checked-in hashes:

```bash
python - <<'PY'
import hashlib
from pathlib import Path
for line in Path('source-reference-run/SHA256SUMS').read_text().splitlines():
    if not line.strip():
        continue
    expected, name = line.split(maxsplit=1)
    path = Path(name.strip())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, path
print('All packet input and atlas hashes match')
PY
```

The following explicit settings record the selected source-profile BASE recipe. From version 0.1.1 these are the default numerical settings; short synthetic runs require explicit overrides. Allow a long-running session; timing depends on hardware and first-run compilation.

```bash
mri_base() {
  pymacruise base-run \
    --t1 source-reference-run/inputs/base-root/T1.nii.gz \
    --brain-mask source-reference-run/inputs/base-root/brain_mask.nii.gz \
    --maseg source-reference-run/inputs/base-root/MaSeg.nii.gz \
    --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
    --output-dir build/tester-human-base \
    --no-regrid-inputs \
    --registration-levels 4 --registration-iterations 50 \
    --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
    --toads-iterations 50 --dura-iterations 5 --dura-membership-level 0.3 \
    --crop-edge-width 5 --ace-diffusion-iterations 25 --gvf-iterations 200 \
    --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
    --tgdm-central-iterations 7 --wm-iso-level 0.5 \
    --ace-use-updated-csf-in-ace 2 --source-faithful-refinement "$@"
}

mri_base 2>&1 | tee build/tester-feedback/human-base.log
# Continue only if the command above succeeds.
mri_base --resume 2>&1 | tee build/tester-feedback/human-base-resume.log
```

Expected: all ten stages complete, including native export, with verified skipping on resume. Previous recorded executions of this example naturally converged at TOADS iterations 20 and 18. Record the actual stop reason and iteration on your machine; completion and convergence are distinct checks. Differences should be reported rather than concealed by lowering thresholds or changing options.

`--no-regrid-inputs` is intentional for this reference-comparison profile. It uses axis permutations/flips and discards oblique directions without interpolation. Orientation/obliquity warnings are expected and should be retained. This recipe is specific to the bundled comparison case. For later testing on other inputs, choose and record the geometry profile deliberately.

If a run is interrupted, retain its output tree and log. In the same checkout and environment, repeat the matching function with `--resume`. Return the original failure even if resume later succeeds. Do not change code, inputs, or options midway through a restart test.

## 4. Verify the results and find files for review

This check verifies all recorded stage hashes, image finiteness and RAS orientation, and prints the final output paths. Run it before moving or archiving the output tree, because manifest paths are absolute. For a synthetic-only check, change `run_dir` to `build/tester-smoke`.

```bash
python - <<'PY' 2>&1 | tee build/tester-feedback/output-check.log
import hashlib
import json
from pathlib import Path
import nibabel as nib
import numpy as np

run_dir = Path('build/tester-human-base')
m = json.loads((run_dir / 'run.json').read_text())
assert m['status'] == 'complete'
assert len(m['stages']) == 9
original = nib.load(m['inputs']['t1']['path'])
images = 0
records = 0
for name, stage in m['stages'].items():
    assert stage['status'] == 'complete', name
    for item in stage['outputs']:
        p = Path(item['path'])
        h = hashlib.sha256()
        with p.open('rb') as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b''):
                h.update(chunk)
        assert h.hexdigest() == item['sha256'], p
        records += 1
        if p.name.endswith('.nii.gz'):
            image = nib.load(p)
            if name == 'native-export':
                assert image.shape == original.shape, p
                assert np.array_equal(image.affine, original.affine), p
            else:
                assert nib.aff2axcodes(image.affine) == ('R', 'A', 'S'), p
            assert np.isfinite(np.asanyarray(image.dataobj)).all(), p
            images += 1
    if 'toads_stop' in stage:
        print(name, stage['toads_stop'])
print('Verified stage records:', records, 'finite images:', images)
print('Converged:', m['converged'])
print('Interpretation:', m['result_interpretation'])
print('Final outputs:')
for item in m['final_outputs']:
    print(item['path'])
PY
```

For version 0.1.1, expect 80 stage-output records and 56 finite images. The 12 `final_outputs` records come from the native-export stage and include five NIfTIs on the original T1 grid, four RAS-mm meshes, summaries/provenance and the stage log. The historical eight-stage counts above describe the earlier revision only.

Use the final paths printed above to locate `SEG_refine.nii.gz`, `inner.nii.gz`, `central.nii.gz`, `outer.nii.gz`, the three `*_surface.vtk` meshes, `inner_surface_thickness.vtk`, `thickness_levelset.nii.gz`, and `thickness_summary.json`. Inspect the active attempts listed in `run.json`; do not assume `attempt-01` after a repair.

Use the original input T1 for overlays of the native exports listed in `final_outputs`. Working intermediates still use their recorded working T1. Final meshes use physical RAS coordinates in millimetres; account for the viewer's coordinate convention.

Record:

- Whether the final refined labels and level-set zero contours align with the original T1 in axial, coronal, and sagittal views.
- Regions with missing cortex, misplaced tissue boundaries, clipping, or unusual disconnected surface pieces; include screenshots and voxel/world coordinates.
- Whether surface and thickness files open and contain usable values.
- Any confusing instruction, unexplained warning, silent pause, or recovery problem.

These anatomical observations are diagnostics. The port's acceptance objective is faithful reproduction of the selected original workflow, including its limitations. Record problems without editing outputs or assuming that improving anatomy establishes source equivalence. A controlled original-workflow comparison is a separate task requiring access beyond this repository.

## 5. Return feedback without GitHub write access

Provide OS/architecture, installed RAM, Python version, elapsed times, pass/fail for installation/synthetic/MRI/resume, and screenshots/observations. Preserve all error text. Include `build/tester-feedback/`, the run manifests, and, when practical, complete run output trees. The complete output archive lets a maintainer inspect the stage logs and scientific results.

```bash
tar -czf build/pymacruise-tester-results.tgz \
  build/tester-feedback build/tester-smoke build/tester-human-base
```

If MRI testing was not attempted, omit `build/tester-human-base` from the archive command. Return the archive through the project owner's agreed transfer channel. A GitHub Issue can carry the summary if your account has permission to create one; do not depend on that permission or on uploading a large output archive to an Issue. No commit or push is needed.

## Troubleshoot a silent stop at N4

`base-run` executes N4 before TOADS and the remaining BASE stages. Input
regridding is preparation, not a completed segmentation workflow. In older
packages, a missing or unimportable ANTsPy dependency can leave an empty
`stages/n4/attempt-01/stage.log`: the nested CLI exits with status 2 while its
error text remains in an unwritten buffer. This symptom alone does not prove
that ANTsPy is the cause on a tester's machine.

From the same activated environment and working directory, expose the error:

```bash
python -c 'import sys, ants; print(sys.executable); print(ants.__version__)'
python -m pip check
python -m pymacruise.cli n4 \
  --input test/outputs/resources/regridded-inputs/T1.nii.gz \
  --output test/n4-diagnostic.nii.gz
```

Use the actual working T1 path if the run did not regrid. The diagnostic
output sits outside the stage checkpoint; do not copy it into an attempt
directory or edit `run.json`. If ANTsPy is missing, install the existing
local wheel with its `[ants]` extra, or, from a source checkout, run
`python -m pip install '.[ants]'`. When repairing an existing installation
without either artifact available, the matching dependency constraints are:

```bash
python -m pip install 'antspyx>=0.6,<0.7' 'numpy>=1.23,<2.3' 'scipy>=1.10,<1.15'
python -m pip check
```

This may download Python dependencies; no separate ANTs executable, original
MaCRUISE installation, or additional atlas is needed. After fixing the
reported problem, repeat the original `base-run` command with `--resume`.
Keep the existing output directory so partial attempts remain available.
If it still fails, return the direct N4 output, dependency checks, and
`run.json`. The source correction prepared October 28, 2026 surfaces nested
CLI exit errors and records both the N4 attempt and run as failed; updating
the code is not required to run these diagnostic commands.

## Earlier execution evidence

- [CLI contract](https://github.com/BennettLandman/pyMaCruise/blob/16993385fde1ef0ed9a7b8e09d667ace740ac18d/docs/cli-help.md)
- [Synthetic run evidence](https://github.com/BennettLandman/pyMaCruise/blob/16993385fde1ef0ed9a7b8e09d667ace740ac18d/docs/reports/synthetic-base-cli-run-2026-10-15.md)
- [Dependency-install evidence](https://github.com/BennettLandman/pyMaCruise/blob/16993385fde1ef0ed9a7b8e09d667ace740ac18d/docs/reports/clean-dependency-install-attempt-2026-10-15.md)
- [Prior full-MRI BASE recipe and limitations](https://github.com/BennettLandman/pyMaCruise/blob/16993385fde1ef0ed9a7b8e09d667ace740ac18d/docs/reports/in-vivo-base-v3-no-mask-n4-source-profile-2026-10-14.md)
- [Source-equivalence acceptance policy](https://github.com/BennettLandman/pyMaCruise/blob/16993385fde1ef0ed9a7b8e09d667ace740ac18d/docs/validation-tolerances.md)
