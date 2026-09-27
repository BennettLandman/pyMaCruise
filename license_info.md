# License and Attribution Working Notes

**Last reviewed:** 2026-09-27  
**Status:** Inventory for later review with the tech transfer office; no final PyMaCRUISE license has been selected.

## Project-owner direction

The project owner states that they are an author on all existing source code supplied for the MaCRUISE/PyMaCRUISE work. Preserve existing license files and meaningful author, revision, and attribution details as the Python port is written, including concise source-reference comments in corresponding Python modules. The tech transfer office will help place an appropriate license over the completed work. No compliance determination or final license choice is needed during implementation.

## Source and notice inventory

| Material | Observed license / attribution information | Preserved project copy and porting note |
|---|---|---|
| MaCRUISE MATLAB, SVN revision 9360 | No project-root LICENSE/COPYING/NOTICE file was found. Source behavior is anchored to `run_maCRUISE_v3` and the selected MaCRUISEplus wrapper. | Source revision and original entry points are recorded in `docs/source-baseline.md` and `docs/plus-baseline.md`. When a Python module follows a MATLAB routine, add a module comment naming the source routine and SVN revision. |
| VandyCRUISE | No repository-root license file was found in the supplied checkout. Individual class headers include author/version information; `ReorientVolume.java` identifies Blake Lucas. | Preserve applicable source author, version, and CVS revision facts in Python module comments and link to this overview. Do not erase relevant algorithm comments when translating behavior. |
| MaCRUISE-specific `toads-cruise` adapters | No repository-root LICENSE/COPYING/NOTICE was found. Many class headers identify Yuankai Huo as author and embed per-file CVS revisions. | Per-file revisions are catalogued in `docs/assets-and-dependencies.md`; retain class-level author/revision attribution in the corresponding Python module docstrings/comments. |
| MIPAV source distribution | `mipav-src/src/license/license.txt` is the MIPAV Source Code Use License. It states that some components have separate licenses, points to `CopyrightInfo.txt` and per-file notices, includes a research-purpose disclaimer, and asks that its license accompany copies or substantial portions of MIPAV software. | Original files copied unchanged to `licenses/upstream/mipav/`, including `license.txt`, `license.html`, and `CopyrightInfo.txt`. PyMaCRUISE does not require the MIPAV runtime; consult the copied notice if MIPAV-derived source is used in a port module. |
| JIST source distribution | `JIST/src/license/LICENSE.txt` is Apache License 2.0. `NOTICE.txt` attributes Apache Commons CLI. The tree also contains component-specific notices/licenses, including a JJ2000 notice and LGPLv3 text in `CopyrightInfo.txt` in the MIPAV tree. | Full supplied `JIST/src/license/` contents copied unchanged to `licenses/upstream/jist/` (including third-party component files). JIST is not a final PyMaCRUISE runtime dependency; retain relevant notices if a source-derived component is used. |
| Bundled TOADS/CRUISE, lesion-TOADS, and SPECTRE assets | Source paths and file hashes are tracked; the SPECTRE dataset includes OAS1-named image files. No standalone license file was found inside the selected atlas folders. The project-owner authorship statement applies to supplied source code; data provenance remains recorded separately from code authorship. | Bundled under `data/atlases/`; source paths, inventory, and checksums are recorded in `docs/assets-and-dependencies.md` and `data/atlases/SHA256SUMS`. |
| Python runtime dependencies | Final package dependencies and versions have not yet been pinned. | Record package name, pinned version, and the license metadata shipped with each selected distribution when the P12 environment is established. |

## Comment and notice practice during the port

For algorithm modules derived from a supplied source behavior, preserve a short attribution comment in the Python file. For example:

```python
"""Python port of the MaCRUISE boundary operation.

Behavior reference: `run_GetMABoundary.m`, MaCRUISE SVN r9360;
adapter reference: `MinFilter3D.java`, toads-cruise CVS r1.5.
Original source author listed in the adapter: Yuankai Huo.
"""
```

Keep comments that explain numerical choices, label definitions, coordinate transforms, and known historical behavior. Do not copy large source blocks merely to preserve attribution; use source references and comments. If an implementation directly incorporates source text, retain its original notice and license context alongside that module and record the incorporation here.

## Tech transfer handoff checklist

- Select the overall PyMaCRUISE license and any separate licensing for bundled source-derived components.
- Review the preserved MIPAV, JIST, and third-party notices against the code actually used in the final Python package.
- Review atlas/data provenance and package notices for all files that remain bundled.
- Inventory final pinned Python dependencies and their notices.
- Add the tech-transfer-approved project license and any required acknowledgments before the chosen release.

This file is a factual working inventory and handoff aid, not a legal analysis or a license grant.
