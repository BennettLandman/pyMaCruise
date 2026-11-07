Subject: PyMaCRUISE BASE 0.1.6 — updated testing package and instructions

Hi Adam and everyone,

Could you test the updated PyMaCRUISE T1-only BASE package?

Download: [INSERT DOWNLOAD LINK]
Package: pymacruise-testers-0.1.6-2026-10-29.tar.gz

The package uses full numerical defaults and exports final volumes on the original T1 grid with integer segmentation labels. It includes the Python wheel, atlas data, installation instructions, a run launcher, output checks, comparison tools, and automatic feedback collection. Original MaCRUISE programs and patient images are excluded.

[FINAL VALIDATION SUMMARY]

Please use a fresh Python 3.13 environment and a new output folder, keeping your previous environment and results for comparison.

1. Download and extract the archive, then open a terminal in the extracted folder. Create a new environment:

   python -m venv .venv-pymacruise-0.1.6

   Activate it on macOS/Linux:

   source .venv-pymacruise-0.1.6/bin/activate

   Or in Windows PowerShell:

   .\.venv-pymacruise-0.1.6\Scripts\Activate.ps1

2. Install the included wheel with the tested dependency constraints:

   python -m pip install --only-binary=:all: -c constraints-tested.txt "./wheels/pymacruise-0.1.6-py3-none-any.whl[ants]"
   python -m pip check
   pymacruise --version

   The version should be 0.1.6.

3. Prepare one input folder containing T1.nii.gz, brain_mask.nii.gz, and MaSeg.nii.gz. The mask and integer-valued MaSeg must match the original T1 grid. Run:

   python run-base.py --input-dir "/path/to/inputs" --output-dir "/path/to/new-output-0.1.6"

   Replace the example paths with your own. The launcher resolves the installed atlas automatically and runs with one numerical-library thread. Allow time for full processing and keep the terminal open. If interrupted, use the same command with --resume; use this only for the same 0.1.6 run and settings.

4. Check the final files listed under final_outputs in run.json, in stages/native-export/attempt-*/. Overlay these on the original T1 and display labels with nearest-neighbor interpolation. Primary surfaces use physical RAS coordinates in millimetres. Original MaCRUISE raw VTKs use voxel indices; the README explains how to create matching voxel-coordinate comparison copies.

5. If you have original MaCRUISE results for the same inputs, please also run:

   python compare-with-macruise.py --run "/path/to/new-output-0.1.6" --original-output "/path/to/MaCRUISEv3.2.0/Output" --surfaces
   python collect-feedback.py --run "/path/to/new-output-0.1.6"

   The surface comparison takes additional time. It reports measured agreement; a successful comparison command can still report failed agreement gates. If you have no original outputs, skip this step.

If installation fails before a run starts, please return the pip error, Python version, and OS/architecture.

Please reply with the generated *-feedback.zip, whether installation/full execution/visual review succeeded, and any alignment, label, or surface problems. The archive contains versions, manifests, checks and logs; it excludes images and meshes. Review the metadata before returning it because local paths remain in the logs. Screenshots with filenames and slice positions are helpful. Keep the run folder in place until checks and feedback collection are complete.

Thank you,
[Your name]
