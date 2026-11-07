# PyMaCRUISE

PyMaCRUISE project code is licensed under the [Apache License 2.0](LICENSE). Third-party materials retain their applicable terms. See the [license and attribution inventory](license_info.md) for upstream notices and the status of bundled atlas data.

**Distribution scope:** The tester bundle is Python-only, with required data and documentation. Original MaCRUISE, MATLAB reference drivers and Java validation harnesses are supplied separately; the full development-source archive is not the tester distribution.

**PLUS testing update — 2026-10-27:** Exploratory end-to-end testing instructions are ready for both original and Python scientific pipelines. Start with the [PLUS walkthrough](docs/testing/full-plus-walkthrough.md), [individual tools](docs/testing/individual-tools.md), and [complete CLI reference](docs/tools-reference.md). The kit prepares the selected `ibsr11objV3` atlas from bundled files. Full original-runtime execution and numerical equivalence remain unverified.

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

For human testing against the complete external MaCRUISE workflow, use the [full BASE walkthrough](docs/testing/full-base-walkthrough.md), including Linux prerequisites, virtual environment, a pinned external container, both full runs, comparisons, and return files. External execution is not yet verified on the authoring host. The [end-to-end comparison protocol](docs/testing/end-to-end-external.md) defines the acceptance boundary and remaining gaps; [Python BASE checks](docs/testing/python-base-checks.md) provide supplemental preflight checks.

## Beta tester user manual

This guide covers installing PyMaCRUISE, running its supported BASE workflow,
interpreting common messages, finding the results, and reviewing images and
surfaces. The command line is the user interface; PyMaCRUISE does not include a
graphical viewer. BASE and PLUS are available for supervised exploratory testing. PLUS atlas
loading and conversion are implemented; independent original-versus-Python
end-to-end agreement remains unverified. Use the PLUS walkthrough for its
explicit assumptions, complete run commands and stage-by-stage inspection.

### BASE and PLUS: scope and intent

**MaCRUISE BASE** is the T1-only cortical reconstruction workflow. Its intended
inputs are a T1-weighted MRI, a corresponding MaSeg label image, and a binary
brain mask on the T1 grid. It runs tissue classification and Dura handling,
builds the inner, central, and outer cortical surfaces, computes thickness
products, and refines the MaSeg labels. Use BASE when evaluating the core
cortical surface and thickness workflow. The Python runner requires the
brain mask from the user; it does not perform the historical BET/SPECTRE
skull-stripping step itself. BASE applies the selected mask-then-windowed-sinc
axial handoff before TOADS and runs Dura in the same LPS processing order.
Final outputs return to the original grid.

**MaCRUISE PLUS** is the lesion-aware extension intended for cases with a
T1-weighted MRI, a corresponding FLAIR image, MaSeg labels, and a T1-grid brain
mask. Its workflow registers FLAIR to T1, uses both modalities in the
LesionTOADS pathway, handles Dura, and produces cortical surfaces, thickness,
and refined labels. The implementation also has an optional AtlasEM/Dura
handoff. PLUS is ready for supervised exploratory evaluation using the prepared
testing kit, including AtlasEM/Dura and final refinement. The native-to-Python
atlas conversion and selected loading are implemented. Full original-runtime
execution, matched selected-atlas validation and complete original-versus-Python
agreement remain unverified. Results are research outputs.

Both workflows are research software for technical evaluation of the selected
MaCRUISE behavior. A completed run means the configured stages finished; it
does not establish source equivalence, anatomical correctness, or clinical
validity. Current implementation evidence and remaining gaps are summarized in
the project status section below.

### Install

PyMaCRUISE supports Python 3.10–3.13. Python 3.13 on macOS is the platform with
the most direct tester execution evidence; other operating systems and clean
dependency installation need tester confirmation. For a BASE run with N4 bias
correction, install the optional ANTs dependency too. A GitHub read-only clone
is sufficient; you do not need to commit or push changes.

```bash
git clone https://github.com/BennettLandman/pyMaCruise.git pymacruise
cd pymacruise
git checkout --detach HEAD

python3.13 -m venv .venv
source .venv/bin/activate
unset PYTHONPATH
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
python -m pip install --upgrade pip
python -m pip install '.[ants]'
python -m pip check
pymacruise --version
pymacruise --help
```

If your Python 3.13 command is named `python3`, use that in place of
`python3.13`. On Windows, create the environment with `py -3.13 -m venv .venv`
and activate it with `.venv\Scripts\Activate.ps1` in PowerShell.
Keep the same checkout, environment, input files, and command options when
resuming a run. First execution may pause while Numba compiles scientific
kernels. MATLAB, MIPAV/JIST, FreeSurfer, and a GPU are not required. See
[supplemental BASE checks](docs/testing/python-base-checks.md) for a complete
installation log, a small smoke run, and the bundled human example.

### Run BASE

The main command is `pymacruise base-run`. It requires a T1-weighted NIfTI,
MaSeg label image, T1-grid binary brain mask, atlas descriptor, and a new
output folder. For the selected source-profile settings and bundled human
example, follow the [Python BASE checks](docs/testing/python-base-checks.md).
That guide includes input hash checks and the complete command. BASE defaults
use full numerical settings and legacy working geometry; final volumes return
to the original T1 grid and meshes use physical RAS+ millimetres. Version 0.1.6
retains mask-before-axial interpolation, LPS Dura processing, v3.2 hemisphere
separation and upward topology correction, with the optional TGDM bound unset.
N4 defaults to the installed backend's current unmasked policy; explicit
`--n4-estimation-mask legacy-otsu` selects historical automatic Otsu estimation.
The policy is recorded and checked during resume. The tester bundle adds a
portable launcher, dependency constraints, original-output comparison and
feedback collection. See the [0.1.6 tester release report](docs/reports/tester-release-0.1.6-2026-10-29.md).

For a quick installation check, run the small synthetic example from the
repository root (its deliberately short iteration limits are not suitable for
judging anatomy):

```bash
mkdir -p build
set -o pipefail
pymacruise base-run \
  --t1 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
  --maseg tests/data/synthetic/base_workflow/maseg.nii.gz \
  --output-dir build/tester-smoke \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 2 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1 \
  2>&1 | tee build/tester-smoke.log
```

The exit status should be zero and all ten stages, including native export, should report completion.
To check restart behavior, repeat the command with `--resume` added before
`2>&1 | tee`, and use a different log filename. Resume verifies completed
stage outputs and skips valid checkpoints. A fresh run needs a fresh output folder. Use
`pymacruise base-run --help` for all supported options and
[the CLI reference](docs/cli-help.md) for input, geometry, and output details.

### Read messages and recover from errors

- **RAS/reorientation warning:** The input was not stored in RAS axis order and
  PyMaCRUISE normalized it. This warning is informative, not automatically a
  failure. Keep it in the log. Working-stage intermediates use their recorded
  working T1; the final native exports restore the original input grid and
  oblique directions, so use the original T1 for their overlays.
- **`Run already exists` / output folder is not empty:** Use `--resume` only to
  continue the same run with the same inputs, settings, and checkout. Choose a
  new output folder for an independent run. Do not delete or overwrite a run
  before saving its logs and results.
- **Resume signature or checkpoint verification error:** Inputs, options,
  assets, code, or a saved output changed. Return to the original checkout and
  use the original command. If you intended a different experiment, start it in
  a new output folder instead of editing the existing manifest.
- **Missing ANTs/`ants` dependency:** Activate the environment and install with
  `python -m pip install '.[ants]'`, then retry. If installation or import
  still fails, return the full installation log and Python/OS versions.
- **BASE stops at N4 with an empty stage log:** Older packages hide the error
  when ANTsPy is missing or cannot import. Run N4 directly to expose it, check
  the active environment, then resume the original command. See the
  [N4 troubleshooting commands](docs/testing/python-base-checks.md#troubleshoot-a-silent-stop-at-n4).
- **Grid/geometry mismatch, invalid mask, or missing-file error:** Check that
  T1, MaSeg, and brain mask are the matching case, that the mask is binary, and
  that the images share the required T1 grid. Do not silence the check by
  changing headers or resampling labels manually; include the exact command
  and error text in the report.
- **TOADS reached its iteration limit / `converged: false`:** This describes
  the algorithm's stopping condition; by itself it does not mean the program
  crashed. Record the stage, stop reason, iteration count, and final `status`
  from `run.json`. A completed workflow and a converged TOADS stage are separate
  facts. The synthetic example is expected to hit its deliberately low limit.
- **Any traceback, nonzero exit, missing output, or failed stage check:** Keep
  the console log and complete output folder. If the run is restartable, retry
  with `--resume` only after saving the original failure. Report both the
  original error and whether resume succeeded.

### Find and understand the outputs

The selected output folder contains a top-level `run.json` manifest, shared
derived inputs under `resources/`, and ordered stage directories under
`stages/`. Each stage has an `attempt-NN/` directory with its products, a
`stage.log`, and often `provenance.json`. A repaired stage may have more than
one attempt: use the active attempt and paths recorded in `run.json`, not an
assumption that `attempt-01` is current. Keep the whole folder together because
the manifest records paths, settings, hashes, and stage status.

| Result | What it contains / how to use it |
| --- | --- |
| `run.json` | Overall status, inputs and hashes, geometry handling, stage status and stop reasons, output paths, and final-output index. `status: complete` means all workflow stages completed; check `converged` separately. |
| `resources/` | Shared prepared or reoriented inputs used by stages. For orientation-only runs, `resources/orientation-only-inputs/T1.nii.gz` is the matching T1 for visual overlays. |
| `stages/*/attempt-NN/` | Intermediate NIfTI images, stage logs, provenance, and stage-specific JSON summaries. These are retained for restart and review. |
| `native-export/.../SEG_refine.nii.gz` | Final integer MaSeg label refinement on the original T1 grid. |
| `native-export/.../{inner,central,outer}.nii.gz` | Level sets; their zero contours define the inner, central, and outer cortical surfaces. |
| `native-export/.../*_surface.vtk` | Scientific surface meshes in RAS+ millimetres; load as model/surface data in a viewer that supports legacy VTK. |
| `native-export/.../thickness_levelset.nii.gz` | Voxelwise thickness product. |
| `native-export/.../inner_surface_thickness.vtk` and `thickness_summary.json` | Per-vertex thickness on the inner mesh and summary measurements. JSON `null` with `available: false` means no valid measurement was available, not a measured zero. |

Saved image products use compressed NIfTI (`.nii.gz`). Final native-export
images preserve the original T1 storage axes and affine; working images use
RAS axes. Surface coordinates are physical RAS+ millimetres. The exact paths for your run are in
`run.json` under `final_outputs` and each stage's `outputs`. Do not move files
individually or rename them before you run a hash/provenance check.

### Visual review

Use a medical-image viewer that can load NIfTI volumes and legacy VTK surfaces.
Open the original input T1 first, then overlay the final native-export `SEG_refine.nii.gz`
with label coloring and inspect axial, coronal, and sagittal slices. Add the
inner, central, and outer level sets as contours or load the corresponding VTK
files as surfaces. Review whether the contours and meshes align to the cortex,
whether there are missing or disconnected regions, and whether the result is
clipped. Load `thickness_levelset.nii.gz` as a scalar image to inspect spatial
variation; the inner mesh's thickness values are in `inner_surface_thickness.vtk`.

For viewers that use GIFTI, convert a workflow surface with
`pymacruise surface-convert --input inner_surface_thickness.vtk --output inner.surf.gii`.
Vertex measurements are included; `--data-output thickness.shape.gii` also
writes them separately. Extraction and surface thickness commands accept GIFTI,
and `surface-info` inspects geometry and vertex fields. See the
[GIFTI guide](docs/gifti.md) for metadata, labels, and coordinate/precision details.

Use the original T1 for products listed in `final_outputs`. Intermediate
working-stage products require their matching working T1 from `resources/`.
VTK surfaces use physical RAS+ millimetres; confirm the viewer uses that frame.
PyMaCRUISE writes scientific data products but
does not render figures or prescribe a viewer. Save screenshots with the
three-plane view, identify the file and stage attempt, and record any concern
with slice or world coordinates. Visual review is diagnostic; it does not by
itself establish equivalence to original MaCRUISE or clinical validity.

### Send beta-test feedback

Include the PyMaCRUISE revision, operating system/architecture, Python version,
available RAM, commands run, elapsed times, and pass/fail for installation,
run, resume, and visual review. Preserve full console logs and error text;
include `run.json`, screenshots, and the complete run folder when practical.
The [supplemental BASE check guide](docs/testing/python-base-checks.md) has a
feedback checklist and archive example. No GitHub write access is required;
return results through the project owner's agreed channel.

## Current status — 2026-10-29

- **BASE:** The installed 0.1.6 candidate passes 128 targeted checks; the tester tools pass eight additional checks. Full defaults, ten stages and original-grid final exports remain enabled. N4 backend selection is explicit, with the current unmasked policy as default. Source-held fusion, topology and ACE comparisons match Adam’s original component products; these remain separate from independent reconstruction agreement. See the [tester release report](docs/reports/tester-release-0.1.6-2026-10-29.md) for final-case measurements and limitations.
- **TOADS performance:** Source-ordered compiled growing, relation, and bias-field kernels reduce the compact full run from 36.5 to 8.3 seconds (4.4×). A fresh human first pass takes 13.5 minutes versus its saved 50.2-minute baseline (3.7×), with all five tissue outputs byte-identical and the same stopping iteration. See [the performance evaluation](docs/reports/toads-compiled-performance-2026-10-17.md).
- **PLUS:** The selected LesionTOADS path is implemented and the eleven-stage AtlasEM-enabled Python workflow has human execution evidence. The testing kit preserves oblique FLAIR geometry, uses the owner-selected atlas, and supports explicit stage checkpoints. Original-run helpers include native-atlas selection, all final products and original thickness calls. Use the [end-to-end walkthrough](docs/testing/full-plus-walkthrough.md) and [readiness receipt](docs/reports/plus-testing-readiness-2026-10-27.md); independent source agreement remains unverified.
- **Synthetic PLUS checkpoint:** The labeled folded-cortex case exercises source-selected surface evolution and explicit ribbon checks. Its current long TOADS profile reaches the selected iteration cap without convergence; held-fixed Java transitions and terminal TGDM comparisons pass their stated gates. This is a nonconvergent synthetic input profile, not evidence that the workflow cannot run. See the [iteration 90–99 report](docs/reports/p41z-iterations90-99-selected-java-2026-10-15.md).
- **Current validation:** See the [October 27 receipt](docs/reports/plus-testing-readiness-2026-10-27.md) for current tests and human-run evidence. Older reports retain their historical results and scope.
- **External references:** The complete original MATLAB/MIPAV/JIST runtime and `reg_aladin` remain unavailable here. Native `ibsr11objV3` assets, Java descriptor preparation, Python conversion, source-run recording and comparison commands are now supplied. Atlas format differences are reported explicitly; they do not establish atlas or whole-workflow equivalence.

## Scope and target

- **Base workflow:** T1-weighted MRI, MaSeg labels, and a user-provided brain mask on the T1 grid.
- **Plus workflow:** T1, MaSeg, matching brain mask, and FLAIR, following the selected `run_maCRUISEplus_v1` profile in [docs/plus-baseline.md](docs/plus-baseline.md).
- **Reference source:** MATLAB MaCRUISE V3 and the supplied VandyCRUISE algorithms. Numbered historical runners are bug-fix versions; versions before V3 contain known bugs.
- **Goal:** run the selected workflows without MATLAB or a MIPAV/JIST JVM runtime. ANTsPy/ANTs is the selected external backend for ANTs operations.
- **Excluded:** visualization porting and FreeSurfer requirements. A post-port discussion will select modern Python visualization tools.

## Implementation principles

1. Use existing Python libraries for general image I/O, array arithmetic, registration, filtering, morphology, and mesh operations when their behavior fits the source; port source-specific behavior with direct evidence where it does not.
2. Preserve useful MaCRUISE outputs, labels, data types, transforms, and source parameters. Treat image calculation as a generic geometry-aware N-dimensional array operator.
3. Save image inputs, fixtures, intermediates, and outputs as compressed NIfTI (`.nii.gz`). Report spatial outputs in RAS; warn and record reorientation of non-RAS inputs.
4. Build and validate on deterministic synthetic data before using the supplied human scan. Keep synthetic fixtures and generators in `tests/data/synthetic/`.
5. Report source-equivalence evidence, Python integration regressions, and real-data validation as distinct evidence categories. Do not infer whole-workflow equivalence from component matches.

The runtime stack is NiBabel, NumPy, Numba, and SciPy; pytest, build, and scikit-image are development dependencies, and `antspyx` is optional for ANTs-backed operations. The license/attribution inventory and follow-up TTO work are in [license_info.md](license_info.md). Atlas assets are regular Git files. The owner-selected PLUS `ibsr11objV3` files are preserved in native form; runtime descriptor preparation and Python conversion are implemented; original-runtime execution and matched validation remain open.

## Inputs and outputs

The command-line input, assumption, and output contract is documented in [docs/cli-help.md](docs/cli-help.md). Coordinate and NIfTI rules are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md). Numerical acceptance thresholds are in [docs/validation-tolerances.md](docs/validation-tolerances.md).
