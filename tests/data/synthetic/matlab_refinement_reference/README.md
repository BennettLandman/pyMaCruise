# Fresh original MATLAB refinement reference

Controlled 13×13×13 synthetic RAS fixture, 1 mm isotropic spacing, created 2026-10-14 for direct execution of the supplied `MAsegRefineBySurf.m` under MATLAB R2024b Update 5. Labels are synthetic WM/cortical/background values; all radius-2 candidate boxes are interior. The fixture exercises WM filling, outer clearing, and central filling, including a label gap (44 and 120 are not rows 1 and 2).

`matlab-refined.nii.gz` was written by the actual original MATLAB routine, not Python. `provenance.json` records source/input hashes and settings. `tabulate-probes.json` was written by MATLAB and covers positive integer uint16/double inputs, missing rows, and ties. Source code remains in the supplied checkout; generated logs are under ignored `build/matlab-refinement-reference-2026-10-14/`.

Original command: MATLAB `-batch` executes the source routine with these three input files and coefficients 0.2 mm, central processing enabled. Add the supplied MaCRUISE `pipeline` and `utils/matlab` paths first. The executable replay script is `run_reference.m` beside these files.

Scope: refinement on held-fixed interior synthetic inputs. This does not validate volume-border handling, historical compiled MATLAB versions, or a complete BASE/PLUS workflow. Data provenance is synthetic; project-generated fixture data contain no subject information.
