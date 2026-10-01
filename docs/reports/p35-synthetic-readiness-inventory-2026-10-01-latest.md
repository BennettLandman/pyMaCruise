# P35 Synthetic Readiness Inventory — 2026-10-01

## Purpose and status

This inventory was refreshed against implementation commit `006482f03a28809b0cf01b92b2f35f6f7b9f9111`. Its complete `pytest -q` run passed **1,028 tests**, with no failures or skips and 11 expected RAS warnings, in 1,290.93 seconds; see the [full-suite report](full-pytest-2026-10-01-lesion-output-io.md). The fixture corpus at this revision has 1,238 tracked files and no untracked files; the filesystem count is 1,240 because two generated cache files are ignored. No subject imaging data was used; the supplied human MaCRUISE run was deliberately held for the real-data gate. Bundled atlas/template files are fixed algorithm inputs, not synthetic phantoms or held-out subject cases. This readiness inventory does not establish whole-workflow equivalence or close P34/P37.

## Fixture and external-data inventory

| Area | Evidence | Readiness interpretation |
| --- | --- | --- |
| Synthetic corpus | `git ls-files tests/data/synthetic` lists **1,238 tracked files**; `git ls-files --others --exclude-standard tests/data/synthetic` finds **0 untracked files**; `find tests/data/synthetic -type f` finds **1,240 filesystem files**, including two ignored cache files; the tree has **297 directories** including the corpus root and contains **527 `.nii.gz` files** and **0 uncompressed `.nii` files**. Recent committed fixtures cover folded one-/two-channel TOADS continuations and the LesionTOADS terminal output-buffer bridge. | The synthetic test corpus is committed and uses `.nii.gz` for NIfTI data. The two ignored files are generated cache files, not fixtures. Counts are inventory evidence, not a measure of algorithm completeness. |
| Subject-data references | `rg -n --glob '*.py' --glob '*.sh' '(/Users/|/Volumes/|/mnt/|sub-[A-Z0-9])' tests` found one textual hit at `tests/test_lesion_initial_membership_small_distance.py:43`: `sub-1e-20f` is a small numerical value, not a subject identifier or path. No `/Users/`, `/Volumes/`, or `/mnt/` paths were found. | No subject-specific path or subject-case fixture was found by this static scan. |
| Test network references | The `requests`/`urllib`/`urlopen` scan found no matches. The `https?://` scan found one literal Git LFS pointer-format string in `tests/test_package_atlas_assets.py`; it is not a network fetch. | No test-source network-fetch call was found. This static scan is not proof from a separately isolated offline environment. |
| Bundled atlas resources | `test_atlas_assets.py` loads and checksum-checks the 35 image resources in `data/atlases/CONVERSION_AUDIT.json`; descriptor-resource tests load bundled TOADS and LesionTOADS metadata. | Atlas images, including the OASIS-named SPECTRE resources, are required fixed algorithm/template inputs and remain the explicit exception to synthetic-only test fixtures. |
| Reference regeneration | Some fixture provenance points to external source checkouts or a JIST jar. For example, `p18_distancefield_jar_reference/README.md` records the JIST jar used to generate one reference. | Normal tests use committed references and do not require regenerating them. Regeneration is separate from the readiness run and requires the recorded environment/source hashes. |
| Optional N4 | `tests/test_n4.py` uses `pytest.importorskip("ants")` for its synthetic N4 comparison. | The full run had zero skips, so this environment executed the N4 test with ANTsPy available. |

## Added component/composition coverage

- **TOADS (P18):** P18cyu composes one-channel membership through selected registration/control and an existing branch adapter. P18cyv then feeds that Python-produced grow/relation state to the next selected thinning adapter; a selected Java thinning reference consumes the same arrays and compares the queue/state outputs. Java does not recompute the grow/relation stage, and the wrapper-order harness is static metadata. These remain bounded slices, not full TOADS.
- **LesionTOADS → plus inputs (P39):** selected post-loop terminal dispatch and a bounded bridge crop the returned GM, WM-fill, and WM-mask channels with seven fixed companions. The Java crop/transfer comparison covers the identity RAS-grid path only; differing physical grids remain outside that comparison. The selected output arrays now also have a typed `.nii.gz` serializer that uses a RAS T1 reference affine and refuses existing destinations; MIPAV `ModelImage` lifecycle equivalence remains outside the port.
- **MaCRUISEplus (P40):** boundary, tissue update, ACE, GVF, and TGDM adapters are composed over supplied crop-ready inputs. GVF receives ACE enhanced GM; TGDM receives a supplied topology-corrected initial membership and the composed tissue, ACE, and GVF products. The test checks non-RAS warning/reorientation and RAS outputs. It does not implement or validate topology correction, crop setup, atlas registration, TOADS, or the complete plus workflow.
- **Registration literature (P59):** the paper guide and code annotations document the selected backend/substitution, transforms, and coordinate conventions without claiming Aladin numerical equivalence.

## Remaining gaps

P18 and P39 remain Started. P34 remains Blocked on assembling and validating the complete TOADS/Dura path. P41 whole-plus assembly, P37 whole-workflow tolerance comparisons, stage acceptance, biological review, and later real-data gates remain open. Surface extraction still uses scikit-image Lewiner marching cubes rather than JIST `IsoSurfaceOnGrid`; the difference is tracked as a risk. The committed synthetic components do not prove whole-workflow equivalence.

## Reproducible inventory commands

Run from the repository root at the implementation commit named above:

```sh
git ls-files tests/data/synthetic | wc -l
git ls-files --others --exclude-standard tests/data/synthetic | wc -l
find tests/data/synthetic -type f | wc -l
find tests/data/synthetic -type d | wc -l
find tests/data/synthetic -type f -name '*.nii.gz' | wc -l
git status --porcelain=v1 --untracked-files=all
rg -n --glob '*.py' --glob '*.sh' '(/Users/|/Volumes/|/mnt/|sub-[A-Z0-9])' tests
rg -n --glob '*.py' --glob '*.sh' '(requests\.|urllib\.|urlopen)' tests
rg -n --glob '*.py' --glob '*.sh' 'https?://' tests
```
