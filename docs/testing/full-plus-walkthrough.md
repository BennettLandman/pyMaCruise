# End-to-end PLUS testing: original MaCRUISE and PyMaCRUISE

Updated October 27, 2026. This is the current testing entry point. The selected Python LesionTOADS implementation is ready for exploratory testing. Original-versus-Python whole-run agreement remains unverified. The original recipe requires an installed MATLAB/MIPAV/JIST runtime; this machine cannot execute that runtime. The BASE v3.2.0 container has not been shown to provide PLUS.

## Test boundary and assumptions

Run both scientific pipelines independently from the same T1, FLAIR, MaSeg and supplied binary brain mask. Preserve all intermediates and failed logs. Atlas conversion is preparation of shared input data; segmentation, registration, memberships, level sets and surfaces from one pipeline must never initialize the other pipeline.

The following choices define this test, rather than claiming every historical product default:

| Choice | Assumption and consequence |
| --- | --- |
| Scientific profile | MaSeg-informed MaCRUISEplus V3 (`run_maCRUISEplus_v1.m` / SVN r9360 source sequence), with the user mask replacing BET/SPECTRE. Figures, FreeSurfer and MNI exports are excluded. |
| AtlasEM | Owner-selected `ibsr11objV3` replaces historical `ibsr8obj`. The kit preserves original source bytes, prepares a Java descriptor with redundant per-map metadata removed, and converts those same raw payloads to Python RAS NIfTI. This is a deliberate runtime adaptation. |
| LesionTOADS | The separate 2012 With-Lesion descriptor remains selected. Initial alignment uses 50 iterations and four pyramid levels; the testing script overrides the Python CLI convenience defaults of 10/two. Main limits are 99, tolerance 1e-4, thinning limit 3. |
| Geometry | Preserve the physical FLAIR affine through rigid registration. Use `--regrid-inputs`; never use the old human run's `--no-regrid-inputs` profile to assess oblique FLAIR anatomy. This kit requires an axis-aligned RAS T1 and matching mask/MaSeg so the comparator can use native T1 coordinates. For other T1 grids, agree on a common physical input grid before preparing a kit. |
| Registration | Python uses deterministic ANTsPy rigid registration; original uses NiftyReg Aladin. Backend equality is unproven. Record transforms and registration overlays separately. |
| AtlasEM covariance | Python uses a 1e-6 covariance floor to handle singular classes. Record this difference from source singular-matrix behavior. |
| Dura transfer | Python uses physical-affine transfer. Original JIST transfer semantics can differ. Retain both handoffs; use `jist-unchanged` only as a separate declared diagnostic run. |
| TGDM | The selected PLUS wrapper passes no optional brain-mask bound. The input mask still controls preprocessing. |
| Acceptance | Per-label Dice >=0.97; each surface symmetric mean distance <=0.01 mm and sampled symmetric Hausdorff <=1 mm. Thickness errors are descriptive until tolerances are agreed. Exploratory agreement does not establish complete port acceptance. |

`testing-kit.json` records these choices and the resolved Python command. Native descriptor spacing/orientation takes precedence over the inconsistent legacy XML headers. The Java descriptor compatibility rewrite removes only per-map Dimensions/Resolutions/Center lines; global geometry, raw voxels and intensity profiles remain intact. Check the installed Java loader against that descriptor and retain its console output.

## 1. Machine, checkout and environment

Use Linux x86-64 for the original runtime. Python also runs on macOS. Plan for 32 GB RAM, at least 60 GB free disk and several hours; these are planning allowances, not measured minima. No GPU is required by the Python workflow. Keep the checkout and input files fixed throughout the experiment.

Use the Python-only `pymacruise-python-testing-2026-10-27.tar.gz` tester bundle, or a checkout containing the reviewed preparation changes. The tester bundle includes the Python wheel, inputs, utilities and guides; it excludes original software, MATLAB reference drivers and Java harnesses. Its `HANDOFF.json` records file hashes. The earlier full development-source snapshot is not the tester distribution. From the extracted bundle or source directory:

```bash
bash
set -euo pipefail
cd /path/to/pymacruise
python3.12 -m venv .venv
source .venv/bin/activate
unset PYTHONPATH
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
python -m pip install --upgrade pip
if [ -d wheels ]; then
  python -m pip install './wheels/pymacruise-0.1.0-py3-none-any.whl[ants,dev]'
else
  python -m pip install '.[ants,dev]'
fi
python -m pip check
python -m pymacruise.asset_inventory
mkdir -p build/plus-e2e-logs
if [ -d .git ]; then
  git rev-parse HEAD > build/plus-e2e-logs/revision.txt
  git diff --binary > build/plus-e2e-logs/working-tree.patch
  git ls-files --others --exclude-standard > build/plus-e2e-logs/untracked-files.txt
else
  cp HANDOFF.json build/plus-e2e-logs/source-snapshot.json
fi
python -m pip freeze > build/plus-e2e-logs/packages.txt
```

Use Python 3.10–3.13; the example uses 3.12. Do not ignore a failed `pip check`; recreate the environment with the declared extras if older packages conflict. The fresh authoring-host installation was verified with Python 3.13, NumPy 2.2.6, SciPy 1.14.1 and ANTsPy 0.6.3. The generated script records the absolute interpreter path. Recreate the kit if moving the checkout or environment. The Python-only tester bundle includes the current installed-package build. Source preparation changes remain in the working tree: a fresh Git clone must include those reviewed changes before using this guide. A diff does not include untracked files; retain a complete checkout/archive when testing an unpublished revision.

Verify the bundled input packet (or use your own correctly paired files with exactly these filenames):

```bash
python - <<'PY'
import hashlib
from pathlib import Path
for line in Path('source-reference-run/SHA256SUMS').read_text().splitlines():
    if line.strip():
        expected, name = line.split(maxsplit=1)
        assert hashlib.sha256(Path(name.strip()).read_bytes()).hexdigest() == expected, name
print('Packet hashes verified')
PY
python scripts/prepare_plus_testing.py \
  --input-dir source-reference-run/inputs/100260 \
  --output-dir build/plus-e2e
export PLUS_KIT="$PWD/build/plus-e2e"
```

The directory must be new or empty. Preparation writes native Java atlas resources, converted Python atlas resources, `adaptation-manifest.json`, `source-config.json`, `testing-kit.json`, and `run-python.sh`. It needs no neighboring VandyCRUISE checkout and does not modify the installed MIPAV atlas. PLUS case 100260 differs from the BASE example; do not mix their files.

## 2. Original scientific pipeline, end to end

Obtain a working original V3+ MATLAB/MIPAV/JIST/VandyCRUISE installation (including Statistics and Machine Learning Toolbox `tabulate` and the source NIfTI `load_untouch_nii_gz`/`save_untouch_nii_gz` helpers), the matching full-multi-atlas helpers and NiftyReg binary, ANTs `N4BiasFieldCorrection`, and the three original LesionTOADS descriptors with all native payloads. Record revisions, plugin hashes, Java/MATLAB versions and exact executable versions. Python's converted `.nii.gz` atlas must not be passed to a Java raw-atlas loader.

Use paths without spaces for the legacy source runtime, inputs and output directory: upstream wrappers assemble shell commands without consistently quoting paths. The selected-mask MATLAB reference drivers are maintained separately from the Python-only bundle. Obtain them with the original testing environment from the project maintainer; `PLUS_MATLAB_HELPERS` must point to that separate driver directory. Set these paths to actual installed resources:

```bash
export MACRUISE_SOURCE=/path/to/MaCRUISE
export MIPAV_INSTALL=/path/to/MIPAV-install
export MULTIATLAS_SOURCE=/path/to/full-multi-atlas
export N4_BINARY=/path/to/ANTs/bin/N4BiasFieldCorrection
export PLUS_PYTHON="$PWD/.venv/bin/python"
export PLUS_MATLAB_HELPERS=/path/to/separately-obtained/reference-drivers
export PLUS_SOURCE_CONFIG="$PLUS_KIT/source-config.json"
test -x "$MIPAV_INSTALL/jre/bin/java"
test -x "$MULTIATLAS_SOURCE/niftyreg/bin/reg_aladin"
test -x "$N4_BINARY"
```

The MIPAV `plugins/Atlas/lesionToads-atlas-2012/` directory must contain `cruise-atlas-12obj-lesiontoads2012.txt`, `cruise-atlas-12obj-toads2012.txt` and `cruise-atlas-10obj-toads2010.txt`, with all referenced files. The kit's separate native AtlasEM path is supplied through `source-config.json`; no `ibsr8obj` installation is needed by the adapted driver.

```bash
set +e
matlab -batch "addpath(getenv('PLUS_MATLAB_HELPERS')); run_plus_selected_mask_reference(getenv('MACRUISE_SOURCE'),getenv('MIPAV_INSTALL'),getenv('MULTIATLAS_SOURCE'),getenv('N4_BINARY'),getenv('PLUS_PYTHON'),getenv('PLUS_SOURCE_CONFIG'));" \
  > "$PLUS_KIT/source-console.log" 2>&1
source_status=$?
set -e
printf '%s\n' "$source_status" > "$PLUS_KIT/source-exit-code.txt"
if [ "$source_status" -ne 0 ]; then
  echo "Original PLUS failed; retain source-console.log and partial source-run."
  exit "$source_status"
fi
```

The adapted driver calls original Aladin, N4, AtlasEM, initial Dura, LesionTOADS, cortical separation, transfer/crop/fusion, ACE/GVF/TGDM, final refinement, volume thickness and closest-point surface thickness. The AtlasEM command retains original scientific options and replaces only its selected atlas path. The thickness adapter calls the original helpers while retaining logs/meshes that the legacy convenience wrapper deletes. It writes `source-run/driver_outputs.json` and `.mat`, including final scientific product paths. Missing final products cause failure. The original stage wrappers can conceal errors; retain their nested logs and run the auditor below.

This is an adapted selected-mask scientific reference, not an unchanged execution of the historical top-level application. The MATLAB helpers have source/wiring checks here; full MATLAB parsing and execution still require the original host. If the installed release emits different filenames, return that failure and its logs rather than substituting old outputs.

## 3. Python scientific pipeline, end to end

Run independently even if collecting a failed original attempt; a source failure prevents the final independent comparison but need not prevent Python testing.

```bash
bash "$PLUS_KIT/run-python.sh" > "$PLUS_KIT/python-console.log" 2>&1
bash "$PLUS_KIT/run-python.sh" --resume > "$PLUS_KIT/python-resume.log" 2>&1
python scripts/audit_python_plus_run.py "$PLUS_KIT/python-run" --report "$PLUS_KIT/python-audit.json"
```

The script runs full settings, including AtlasEM/Dura and final refinement. Initial Numba compilation can take time. Completion requires all eleven stage records `complete`, top-level `status: complete` and `terminal: true`. The Python-only auditor verifies output integrity and lists any nonfinite intermediate arrays for review; final products must be finite. The primary record is `python-run/plus-full-workflow.json`. The final segmentation is indexed by `final_outputs.SEG_refine`; level sets, meshes and thickness are indexed under the child `plus_surface_workflow` manifest. Use those active paths rather than guessing attempt/revision folders.

The same script accepts `--stop-after-stage STAGE` and `--resume` for stage inspection; see [individual tools](individual-tools.md). A partial prefix is recorded as `checkpointed`, with `terminal: false`, and must never be reported as a finished run. Resume requires unchanged inputs/settings; source changes preserve the old generation and recompute into a new one.

## 4. Audit, record and compare

Audit the original run before recording it. The auditor hashes regular files, inspects NIfTI headers/finite values, and reports failures. Fill in `RUN_INFO.md` using the [template](../../source-reference-run/RUN_INFO.template.md), including the actual launch command and effective settings.

```bash
cp source-reference-run/RUN_INFO.template.md "$PLUS_KIT/source-run/RUN_INFO.md"
# Fill in RUN_INFO.md with actual runtime/source/atlas and stopping records.
python scripts/audit_source_reference_run.py "$PLUS_KIT/source-run" \
  --require-run-info --require-ras
python scripts/record_plus_reference.py \
  --kit "$PLUS_KIT/testing-kit.json" --source-root "$PLUS_KIT/source-run" \
  --lesion-descriptor "$MIPAV_INSTALL/plugins/Atlas/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt" \
  --log "$PLUS_KIT/source-console.log" \
  --exit-code "$(cat "$PLUS_KIT/source-exit-code.txt")" \
  --runtime 'REPLACE with MATLAB, Java, MIPAV, JIST and plugin identities' \
  --source-revision 'REPLACE with actual source revision and hashes' \
  --command 'REPLACE with exact MATLAB launch command' \
  --mesh-frame voxel
```

`record_plus_reference.py` reads original products, the source driver record and the shared input kit; it does not read Python results. It creates the nine-product map and reference-info JSON automatically and records all positive starting MaSeg labels as required labels. The three `REPLACE` strings must be replaced with actual information before use. Mesh frame `voxel` matches the reviewed `MaTGDMcrop` export (crop-offset translation of voxel coordinates). Confirm that against the installed source; if it exports physical RAS millimetres, use `ras-mm` here and below. The comparator applies the known T1 affine for voxel meshes; it never fits a transform between result meshes.

```bash
set +e
python scripts/compare_external_plus_run.py \
  --external-root "$PLUS_KIT/source-run" --python-run "$PLUS_KIT/python-run" \
  --t1 source-reference-run/inputs/100260/T1.nii.gz \
  --flair source-reference-run/inputs/100260/FLAIR.nii.gz \
  --brain-mask source-reference-run/inputs/100260/brain_mask.nii.gz \
  --maseg source-reference-run/inputs/100260/MaSeg.nii.gz \
  --reference-info "$PLUS_KIT/source-run/reference-info.json" \
  --reference-log "$PLUS_KIT/source-console.log" \
  --product-map "$PLUS_KIT/source-run/product-map.json" \
  --source-mesh-frame voxel --atlas-policy report-differences \
  --report "$PLUS_KIT/comparison.json"
comparison_status=$?
set -e
printf '%s\n' "$comparison_status" > "$PLUS_KIT/comparison-exit-code.txt"
```

For your own case, replace the four input paths with the same files used to prepare the kit. Exit 0 means segmentation/surface numeric gates passed, 1 means a numeric gate failed, and 2 means the tool could not perform a valid comparison. Every report leaves overall port acceptance unverified; inspect thickness, transforms, preprocessing, atlas representation and omitted comparisons separately.

`--atlas-policy strict` is the default and rejects different descriptor/payload bytes. Java raw versus Python NIfTI resources necessarily differ. The explicit exploratory `report-differences` policy records these discrepancies and computes product metrics without claiming atlas equivalence. It never relaxes input identity, output hashes, physical grids or numeric thresholds. Preserve native and adapted atlas manifests for a separate representation audit. Do not promote an exploratory numeric pass to source-equivalence acceptance.

Review registered FLAIR/T1 overlays and each final segmentation/surface in all three planes. Inspect the earliest differing intermediate before attributing a final difference to LesionTOADS or TGDM. Return the whole kit, complete source/Python trees, exit-code files, command/environment records, atlas provenance, comparison JSON and screenshots. A failed reference run is a valid test finding, not an oracle.

## Current evidence and limits

The October 18 Python human run completed/resumed eleven stages and naturally converged at lesion iteration 26, but used orientation-only FLAIR handling. It is historical execution evidence. The October 27 validation receipt is maintained in [the readiness report](../reports/plus-testing-readiness-2026-10-27.md). Original full-runtime execution, Aladin agreement, Java-versus-Python complete LesionTOADS agreement and whole-PLUS agreement remain unverified. No clinical-validity claim is made.
