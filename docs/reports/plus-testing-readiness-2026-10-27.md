# PLUS exploratory testing readiness — October 27, 2026

The owner accepts the selected LesionTOADS implementation scope for exploratory testing and requests the best available independent end-to-end handoff. This receipt distinguishes testing readiness from proof of original-source agreement.

## Delivered testing path

- [Full original/Python walkthrough](../testing/full-plus-walkthrough.md): clean environment, shared raw inputs, explicit assumptions, both complete scientific runs, audits, product mapping, comparison and return records.
- [Individual tools](../testing/individual-tools.md) and [generated CLI reference](../tools-reference.md): all command usages/defaults, separate scientific tools and source helpers, and safe inspection at every full-workflow stage boundary.
- `scripts/prepare_plus_testing.py`: unpacks the checksum-verified native `ibsr11objV3` resources from the package, prepares a Java descriptor compatible with consecutive Structure/File loading, converts the same raw maps for Python, records assumptions and writes a fixed-settings run script. No neighboring atlas checkout is needed.
- Original driver: explicit atlas selection through `source-config.json`, complete final products, original volume/closest-point thickness helpers, preserved logs and source product JSON. Original algorithms and installed atlases remain untouched.
- `scripts/record_plus_reference.py`: creates the explicit nine-product map and source run record from original outputs plus shared input provenance; it never reads Python segmentation products.
- `scripts/audit_python_plus_run.py`: verifies complete/terminal stage order, hashes, finite RAS outputs, common final grids, integer segmentation and readable nonempty meshes without needing a source run.
- Comparator: strict atlas identity remains the default; exploratory `report-differences` records native/converted atlas discrepancies without relaxing image identity, output integrity or numeric gates. Explicit voxel versus RAS-mm source mesh coordinates avoid guessing or fitting a result transform.
- A concise [tester email](../testing/tester-email.md) is included; replace its shared-download placeholder before sending. No email has been sent.
- Active README, source packet, source contract, plan status and testing/I/O guides are refreshed. Dated reports/handoffs retain their original results and are marked historical where they contain superseded readiness statements.

## Assumptions

The scientific boundary uses the supplied binary brain mask in place of BET/SPECTRE and excludes visualization/MNI/FreeSurfer exports. AtlasEM uses owner-selected `ibsr11objV3`; LesionTOADS keeps its separate With-Lesion 2012 descriptor. The native atlas descriptor's geometry governs conversion, rather than conflicting legacy XML headers. Java preparation removes redundant per-map dimensions/resolutions/centers, preserving global geometry, intensity profiles and raw payload bytes. The installed original runtime must still verify its loader behavior.

The testing kit uses geometry-preserving `--regrid-inputs`, initial LesionTOADS alignment 50 iterations/four levels, 99-iteration limits, energy tolerance 1e-4, thinning limit 3 and physical-affine Dura transfer. The kit restricts T1 to an axis-aligned RAS physical grid shared with mask/MaSeg; oblique FLAIR retains its input affine during registration. Other T1 geometry needs an agreed common input grid. ANTsPy substitutes for original Aladin; AtlasEM's 1e-6 covariance floor and physical-affine transfer remain explicit source differences. PLUS TGDM receives no optional brain-mask bound. Tests do not establish clinical validity.

## Validation receipt

Validation uses the working tree based on `cbff9b09` (October 19 beta manual). Existing uncommitted PLUS work was preserved. The handoff source snapshot includes this complete working state; changes are not pushed or published.

| Check | Result and scope |
| --- | --- |
| Development full suite | **1,711 passed, four skipped, zero failed**, 30 orientation warnings, 683.33 s; 11.408 GiB peak sampled RSS under the 16 GiB watchdog. This older development environment has an ANTsPy/SciPy version mismatch and is supporting evidence. |
| Fresh installation and declared build dependencies | Isolated noneditable `.[ants,dev]` install succeeded; the clean distribution test exposed missing development build tools, so the `dev` extra now explicitly includes setuptools/wheel and the distribution test passes after that repair; `pip check` reports no broken requirements. Python 3.13.17, NumPy 2.2.6, SciPy 1.14.1, ANTsPy 0.6.3. The imported package path is the isolated environment's site-packages, rather than the source checkout. |
| Clean installed-package workflow/tool tests | **23 passed**, 58.51 s, ten warnings. Four warnings exercise deliberate orientation-only diagnostics; six floating-point warnings occurred during synthetic mesh mapping, whose finite output checks passed. Includes full synthetic PLUS execution, prefix checkpoint/resume and final refinement. |
| Final recorder/auditor/kit/comparator checks | **13 passed**, 5.76 s; source-driver static checks passed again after the final dependency preflight was added. A prior corruption-test assertion expected a different error message; it was corrected and the clean focused set rerun successfully. |
| Clean installed-package full suite and targeted repair | **1,712 passed, four skipped, one failed**, 105 warnings, 670.42 s; 11.434 GiB peak sampled RSS. The sole failure was the distribution build test because setuptools/wheel were missing from the declared development environment. Both are now in the `dev` extra; after installation the full distribution test **passed in 25.80 s**, including sdist/wheel creation, isolated wheel installation, atlas verification and workflow preflight. No scientific test failed. The full suite was not repeated after this metadata-only repair; the initial failure and targeted passing recheck remain separate receipts. |
| Documentation QA | Twelve active guides checked for local links; 27 Bash examples pass `bash -n`; 16 standalone invocation examples (including the separate LesionTOADS exporter) parse with the current CLI. |
| Fresh human preflight | Case 100260 completes geometry-preserving registration/N4, AtlasEM/initial Dura, crop and 50-iteration/four-level atlas alignment. Longer continuation was intentionally interrupted during registered recrop; cached prefixes were verified with `--resume --stop-after-stage initial_alignment`, exit zero. Parent is `checkpointed`, `terminal: false`; the interrupted later child is not claimed complete. Preprocessing peaked at 3.273 GiB, extended preflight at 4.527 GiB; final prefix verification at 2.112 GiB. |
| Human registration visual QA | Three-plane T1/registered-FLAIR checkerboards were rendered and inspected for gross alignment; this supports geometry inspection rather than source agreement or anatomical accuracy. |
| Existing complete human output integrity | The new auditor verifies the historical October 18 eleven-stage run: **101 hashed files**, **77 finite RAS NIfTIs**, three readable meshes, no nonfinite intermediate images. This retains the historical run's oblique-FLAIR limitation and is not a new full human run. |

Logs, registration QA and machine-readable receipts are under `build/plus-testing-validation-2026-10-27/`. The generated local kit/checkpoint is `build/plus-testing-2026-10-27-v2/`. The new complete human segmentation was deliberately left for the tester's end-to-end exercise; a local continuation can use the same `run-python.sh --resume` once no previous process is running. Atlas format/backend/physical-grid limits are explicit in the walkthrough. No background validation process is intended to remain after preparation finishes.

## Remaining external boundary

MATLAB, original installed MIPAV/JIST runtime and NiftyReg Aladin are unavailable on the authoring host. Source-driver checks are static/wiring evidence, not MATLAB execution or a successful original full run. The original stage wrappers may conceal nested failures; the handoff requires nonzero status handling, required products, finite/header audits and retained nested logs. A different installed runtime may need output-name adaptation, which must be reported with the original failure.

A fresh original PLUS run and independent Python result have not been compared here. Native/Python atlas representations, Aladin/ANTsPy registration, composed AtlasEM and complete LesionTOADS numerical agreement remain unverified. Exploratory final-product metrics cannot close those claims, and no universal thickness/scalar acceptance tolerance has been established.

## Transfer package

The current tester distribution is `build/pymacruise-python-testing-2026-10-27.tar.gz`: Python wheel (including required atlas input data), sample inputs, Python preparation/audit/comparison utilities, license notices and guides. It excludes original MaCRUISE software, MATLAB reference drivers, Java source/harnesses, JARs and containers. Original software and any selected-mask reference drivers must be supplied separately by the original testing environment. `HANDOFF.json` inventories every bundled file; the separate receipt records its SHA-256 and the zero-original-code audit.

The earlier `build/pymacruise-plus-testing-2026-10-27.tar.gz` is a full development-source snapshot containing validation/reference code. It is superseded for tester distribution and retained only as an internal development artifact. No original complete runtime was in that archive, but its contents were broader than the intended Python-only tester package. Changes remain local and unpushed.
