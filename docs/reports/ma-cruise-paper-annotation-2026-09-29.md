# P48: MaCRUISE and CCSE paper mapping

## Citation

Huo Y, Plassard AJ, Carass A, Resnick SM, Pham DL, Prince JL, Landman BA.
“Consistent cortical reconstruction and multi-atlas brain segmentation.”
*NeuroImage*. 2016;138:197–210. DOI
[10.1016/j.neuroimage.2016.05.030](https://doi.org/10.1016/j.neuroimage.2016.05.030);
PMID 27184203; PMCID PMC4927397. The paper is available from
[PubMed](https://pubmed.ncbi.nlm.nih.gov/27184203/) and
[PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC4927397/).

## Paper-to-code map

- Section 2.4 describes Cortical Consistent Segmentation Editing (CCSE), which
  reconciles volumetric labels with inner/outer CRUISE surfaces using their
  level sets and consistency distances. The selected MaCRUISE V3 helper is
  `src/pymacruise/segmentation_refine.py`; it edits labels outside the outer
  surface, on the WM side of the inner surface, and between surfaces.
- The function docstring and
  [`docs/segmentation-refinement.md`](../segmentation-refinement.md) cite that
  paper section and explain the mapping to the selected implementation.
- The paper's §2.2 includes NLSS multi-atlas segmentation and segmentation
  fusion with TOADS. The workflow currently accepts MaSeg as an input; neither
  NLSS nor the paper's fusion algorithm is implemented by this refinement
  function.
- MaACE is described separately in §2.3.1 and remains the distinct P51 paper
  mapping. The term is **MaACE**.

## Source-specific differences

The selected V3 wrapper and `MAsegRefineBySurf.m` define four 0.2 mm
coefficients, the label groups, a 5×5×5 neighborhood, the update order, and
the MATLAB `maxind` assignment behavior. The paper describes its typical
evaluation using 0.5 mm inner/outer consistency coefficients and discusses a
26-connected iterative cortical-label fill. Those paper settings are not
silently substituted for V3 source behavior. Python clips source neighborhood
windows at the image boundary; the original MATLAB helper may fail there.
Python's `source_faithful` setting preserves the historical `maxind` result by
default; corrected modal-label output remains opt-in pending owner approval.

## Validation evidence and limits

Synthetic refinement validation remains documented in
[`docs/segmentation-refinement.md`](../segmentation-refinement.md). The full
suite passed: `.venv/bin/python -m pytest -q` — **365 passed in 76.38 seconds**
on commit `0ecf26f5327253329f23bfc40c29ab62318f061a`; see
[`full-pytest-2026-09-29.md`](full-pytest-2026-09-29.md). The citation and
comments explain method lineage and source mapping; they do not establish
whole-workflow MaCRUISE equivalence.
