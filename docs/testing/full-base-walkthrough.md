# Full BASE test: external MaCRUISE versus PyMaCRUISE

Run this from a fresh clone using a normal account with read-only GitHub access. No fork, commit, or push is required. This runs the complete public BASE container and the complete Python BASE workflow, then produces numerical comparisons. The small synthetic case is not used as the primary test.

**Reference selected by the project owner:** MaCRUISE **v3.2.0**, provided by [`masidocker/harmonizeapp:macruise`](https://hub.docker.com/r/masidocker/harmonizeapp), pinned to the registry digest verified October 16, 2026. The actual launch scripts were read from the immutable image's small script layers to determine the command below. See the [inspection record](../reports/external-v320-reference-inspection-2026-10-16.md). This supersedes the initial v3.1.0 recipe.

**Execution status:** the external tag/digest is available; the comparison tool has nine passing synthetic tests. Full external container execution and fresh Linux dependency installation have not been executed on the authoring host. Treat failures as test findings and retain logs. Python has full-human BASE execution evidence, but overall external-versus-Python equivalence is not established.

## 1. Computer and prerequisites

Use an x86-64 Linux workstation or server; Ubuntu 24.04 with Python 3.12 is a concrete starting environment. Budget 32 GB RAM and approximately 60 GB free disk space for the image, both output trees, and the result archive. These are planning resources, not measured minimum requirements. Allow several hours for the external workflow. The T1 and initial MaSeg segmentation are already bundled, so do not rerun SLANT as part of this test.

The Python package supports Python 3.10 through 3.13. Use 3.12 for these Ubuntu instructions. Do not use Python 3.14 or later. The invocation retains the upstream Singularity `--nv` flag; record any NVIDIA library/driver setup warnings. This recipe does not run SLANT. GPU-runtime behavior on a CPU-only host has not been verified here. Do not use Apple Silicon emulation for the initial reference test.

On Ubuntu, install the basic tools:

```bash
sudo apt-get update
sudo apt-get install -y git python3 python3-venv ca-certificates curl
python3 --version
uname -m
```

If Docker is already installed and `sudo docker info` succeeds, skip Docker installation. Otherwise use the current [official Docker Engine Ubuntu instructions](https://docs.docker.com/engine/install/ubuntu/). On a fresh supported Ubuntu host with no conflicting container packages, the installation is:

```bash
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

. /etc/os-release
printf 'Types: deb\nURIs: https://download.docker.com/linux/ubuntu\nSuites: %s\nComponents: stable\nArchitectures: %s\nSigned-By: /etc/apt/keyrings/docker.asc\n' \
  "$VERSION_CODENAME" "$(dpkg --print-architecture)" \
  | sudo tee /etc/apt/sources.list.d/docker.sources

sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo systemctl start docker
sudo docker run --rm hello-world
```

Have your system administrator handle Docker installation on a managed server. The virtual environment below needs no administrator permission. These commands use `sudo docker`; if your site provides Docker access without sudo, consistently replace it with `docker`.

## 2. Clone and create the virtual environment

Use a new parent folder and a Bash shell. Run all later commands from the repository root in this shell.

```bash
bash
set -euo pipefail

mkdir -p "$HOME/pymacruise-testing"
cd "$HOME/pymacruise-testing"
git clone https://github.com/BennettLandman/pymacruise.git
cd pymacruise
git checkout --detach HEAD

python3 -m venv .venv
source .venv/bin/activate
unset PYTHONPATH
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mkdir -p build/e2e/logs build/e2e/external-input build/e2e/external-output
python -m pip install --upgrade pip
python -m pip install '.[ants,dev]' 2>&1 | tee build/e2e/logs/install.log
python -m pip check
python -c 'import numpy, scipy, nibabel, numba, ants; print("All scientific dependencies imported")'
pymacruise --version

git rev-parse HEAD > build/e2e/logs/python-revision.txt
python --version > build/e2e/logs/python-version.txt
python -m pip freeze > build/e2e/logs/packages.txt
uname -a > build/e2e/logs/system.txt
python -m pytest -q tests/test_external_base_comparison.py \
  2>&1 | tee build/e2e/logs/comparator-tests.log
```

Expected: dependencies import, `pip check` succeeds, and nine comparator tests pass. Stop if installation or these checks fail; preserve the log. Do not update the checkout during the experiment. To continue in a new terminal later, return to this repository, run `source .venv/bin/activate`, and restore the environment exports above. Run `deactivate` when finished.

## 3. Verify and stage the same input scan

The BASE example is `source-reference-run/inputs/base-root/`. The `100260` packet is a different scan and must not be mixed into this test.

```bash
python - <<'PY'
from pathlib import Path
import hashlib
for line in Path('source-reference-run/SHA256SUMS').read_text().splitlines():
    if line.strip():
        expected, name = line.split(maxsplit=1)
        path = Path(name.strip())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, path
print('Packet hashes verified')
PY

cp source-reference-run/inputs/base-root/T1.nii.gz build/e2e/external-input/T1.nii.gz
cp source-reference-run/inputs/base-root/MaSeg.nii.gz build/e2e/external-input/orig_target_seg.nii.gz
```

Both workflows start from those same T1 and MaSeg arrays. The external application performs its normal skull stripping. Python additionally receives the bundled brain mask, as required by its supported interface. Preserve this difference in the report: final disagreement can include preprocessing differences. Do not feed external N4, Dura, segmentation, or surface intermediates into Python for this primary run.

## 4. Download and run the complete external v3.2.0 application

The image contains a BIDS convenience wrapper around a complete embedded Singularity application. Its wrapper stages `T1.nii.gz` and `orig_target_seg.nii.gz`, executes the application, copies outputs, and deletes its temporary input/output folders. It does not propagate the scientific application's failure status reliably.

For this single-scan comparison, call the same complete embedded application directly, with the same scientific invocation. This bypasses BIDS discovery/copy/cleanup only; it does not replace N4, skull stripping, TOADS, Dura, surfaces, refinement, or any other scientific stage. Both images start from raw T1/MaSeg. The application exit status is preserved and outputs remain on the host. This is not an unchanged invocation of the outer BIDS wrapper; its complete scientific application remains unchanged.

The nested Singularity runtime uses privileged Docker, as in the published BIDS-container usage. Use a host on which your administrator permits that runtime setup. Download the immutable image, retain its metadata, and record the launch script:

```bash
export REFERENCE_IMAGE='masidocker/harmonizeapp@sha256:96c6c7a033edb86f2be7b3a60c3e8cc658691de3f670b35a4147feaf287a28a6'
sudo docker pull --platform linux/amd64 "$REFERENCE_IMAGE" \
  2>&1 | tee build/e2e/logs/external-pull.log
sudo docker image inspect "$REFERENCE_IMAGE" > build/e2e/logs/reference-image.json
sudo docker version > build/e2e/logs/docker-version.txt

cat > build/e2e/run-external.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
sudo docker run --rm --privileged --platform linux/amd64 \
  --entrypoint singularity \
  -v "$PWD/build/e2e/external-input:/INPUTS/" \
  -v "$PWD/build/e2e/external-output:/OUTPUTS/" \
  "$REFERENCE_IMAGE" \
  exec --nv -B /INPUTS/:/INPUTS -B /OUTPUTS/:/OUTPUTS \
  /singularity/macruise_classern_v3.2.0.simg \
  xvfb-run -a --server-args="-screen 3 1920x1200x24 -ac +extension GLX" \
  /extra/MaCRUISE_v3_2_0_classern /extra/MaCRUISE_v3_2_0_classern
SH

set +e
bash build/e2e/run-external.sh 2>&1 | tee build/e2e/logs/external.log
external_status=${PIPESTATUS[0]}
set -e

python - "$external_status" "$REFERENCE_IMAGE" <<'PY'
import json
from pathlib import Path
import sys
Path('build/e2e/logs/reference-run.json').write_text(json.dumps({
    'workflow': 'BASE',
    'image': sys.argv[2],
    'exit_code': int(sys.argv[1]),
    'entrypoint': '/extra/MaCRUISE_v3_2_0_classern',
    'embedded_image': '/singularity/macruise_classern_v3.2.0.simg',
    'launch_mode': 'direct complete scientific application; outer BIDS discovery/copy/cleanup bypassed',
    'preprocessing': 'unchanged external application; internally generated skull stripping',
    'launch_script': 'build/e2e/run-external.sh',
    'console_log': 'build/e2e/logs/external.log'
}, indent=2) + '\n')
PY

sudo chown -R "$(id -u):$(id -g)" build/e2e/external-output
if [ "$external_status" -ne 0 ]; then
  echo "External workflow failed. Preserve build/e2e and return its logs; do not report a completed reference run."
  exit 1
fi

find build/e2e/external-output -name SEG_refine.nii.gz -print
```

The repeated application path at the end of the command is intentional: it matches the inspected v3.2.0 launcher. Do not replace it with the v3.1.0 executable.

Do not rely on a completion message alone. The final native output tree must include `SegRefine/SEG_refine.nii.gz`, `Raw/T1.nii.gz`, `Raw/Seg.nii.gz`, all three native `Surfaces/target_image_GMimg_*Vol.nii.gz` level sets and `*Surf.vtk` meshes, and the scientific thickness NIfTI/MAT products. The comparison command checks their presence. Preserve the original full output tree, including failed-stage logs. Rendering/PDF failures must be reported even if scientific products exist; a nonzero exit is not silently treated as success.

## 5. Run the complete Python BASE workflow

These explicit options use the selected V3 scientific settings. Do not omit them: several command defaults are intentionally small synthetic-test values.

```bash
cat > build/e2e/run-python.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
pymacruise base-run \
  --t1 source-reference-run/inputs/base-root/T1.nii.gz \
  --brain-mask source-reference-run/inputs/base-root/brain_mask.nii.gz \
  --maseg source-reference-run/inputs/base-root/MaSeg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/e2e/python-base \
  --no-regrid-inputs \
  --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 5 --dura-membership-level 0.3 \
  --crop-edge-width 5 --ace-diffusion-iterations 25 --gvf-iterations 200 \
  --tgdm-initial-iterations 6 --tgdm-inner-iterations 4 \
  --tgdm-central-iterations 7 --wm-iso-level 0.5 \
  --ace-use-updated-csf-in-ace 2 --source-faithful-refinement "$@"
SH

bash build/e2e/run-python.sh 2>&1 | tee build/e2e/logs/python-base.log
bash build/e2e/run-python.sh --resume 2>&1 | tee build/e2e/logs/python-resume.log
```

Expected: eight complete stages; then eight `Resume: verified ...; skipping` messages. Both TOADS stages should report whether they naturally converged or reached an iteration limit. Previous executions of this bundled case stopped naturally at iterations 20 and 18; record your actual values. If Python is interrupted, rerun its saved script with `--resume`, preserving the failure log.

`--no-regrid-inputs` is intentional here. Oblique directions are discarded while values and voxel ordering remain unchanged. Retain warnings and use the run's working T1 for visual overlays. The comparison restores the Python meshes to the original input's physical frame through the recorded voxel mapping.

## 6. Compare final scientific outputs

Run this before moving the output directories: the Python manifest contains absolute paths. The adapter follows the active attempts in that manifest and verifies all recorded Python stage-output hashes. It finds the native external output tree rather than requiring a machine-specific path.

```bash
set +e
python scripts/compare_external_base_run.py \
  --external-output build/e2e/external-output \
  --python-run build/e2e/python-base \
  --t1 source-reference-run/inputs/base-root/T1.nii.gz \
  --maseg source-reference-run/inputs/base-root/MaSeg.nii.gz \
  --reference-info build/e2e/logs/reference-run.json \
  --reference-log build/e2e/logs/external.log \
  --source-mesh-frame voxel --source-volume-grid raw-t1 \
  --report build/e2e/comparison.json \
  2>&1 | tee build/e2e/logs/comparison.log
comparison_status=${PIPESTATUS[0]}
set -e
printf 'Comparison exit code: %s\n' "$comparison_status"
```

The explicit native source export contract here is: source level-set/thickness arrays retain the source Raw/T1 voxel ordering, and the source native VTK meshes use voxel coordinates. This is the contract used by the existing downloaded-source comparison. The adapter verifies raw input identity, shapes, Python provenance and the axis-only mapping; it does not independently prove the external export's voxel-order convention. If the pinned container exports a different layout/frame, stop and return the products for adapter correction. Do not change frames or register outputs simply to improve metrics. Use `--source-volume-grid physical` only when both image headers describe the same actual physical grid.

Read `build/e2e/comparison.json`:

- Exit **0**: final label and three surface gates pass. This alone does not establish overall port acceptance.
- Exit **1**: comparisons completed, but at least one label/surface gate failed. This is a valid test result; retain and return it.
- Exit **2**: missing files, hash/input/grid problems, failed external execution, or another comparison error. Inspect `comparison_error`; do not interpret it as a numeric pass or fail.

The report includes every final nonbackground label's Dice, level-set differences, all three symmetric vertex-to-triangle surface distances in millimetres, volume-thickness errors/support Dice, and surface-thickness distribution summaries. Gates are Dice ≥0.97 for every label, surface mean ≤0.01 mm and sampled Hausdorff ≤1 mm for each surface. Scalar/thickness limits are not universally defined; vertexwise thickness correspondence between independently generated meshes remains unverified. Those limitations and preprocessing differences keep overall port acceptance explicitly unverified. The report is intended to expose end-to-end differences, not guarantee a pass.

## 7. Visual review and return package

Open the native outputs in your existing NIfTI/VTK-capable viewer. Review refined segmentation and inner/central/outer boundaries in axial, coronal, and sagittal planes. Capture discrepant regions with coordinates and screenshots. Use the Python working T1 under `python-base/resources/orientation-only-inputs/` when present; use the source Raw/T1 for external outputs. Do not overlay the axis-only Python outputs directly on an oblique source grid without the documented mapping.

Create `build/e2e/TESTER_NOTES.md` recording hardware/RAM, elapsed times, external/Python/resume completion, actual stop reasons, comparison exit code, viewer/frame choices, and observed problems. Put screenshots under `build/e2e/screenshots/`. Keep both complete output trees, logs, launch scripts, and the JSON report, including failures.

```bash
tar -czf build/pymacruise-full-base-results.tgz -C build e2e
deactivate
```

Return that archive through the project owner's agreed file-transfer channel. No GitHub write access is required. If a step fails before the comparison, archive the partial `build/e2e` tree and return it with the failure description.

## PLUS

This is a complete **BASE** test. The current Python PLUS route still lacks the full selected AtlasEM/pre-LesionTOADS Dura sequence and required atlas payload; substituting the existing nine-stage run would not constitute a complete PLUS reproduction. The external build, registration boundary, assets, and full Python PLUS route need to be resolved before a corresponding PLUS walkthrough can claim that scope. See [the source-stage audit](../reports/plus-source-stage-coverage-audit-2026-10-16.md).
