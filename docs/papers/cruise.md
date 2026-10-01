# CRUISE paper-to-code notes

## Published method

Han X, Pham DL, Tosun D, Rettmann ME, Xu C, Prince JL. “CRUISE: cortical
reconstruction using implicit surface evolution.” *NeuroImage*.
2004;23(3):997–1012. DOI
[10.1016/j.neuroimage.2004.06.043](https://doi.org/10.1016/j.neuroimage.2004.06.043);
[PubMed](https://pubmed.ncbi.nlm.nih.gov/15528100/);
[author-hosted paper](https://iacl.ece.jhu.edu/~prince/pubs/2004-NI-Han-CRUISE.pdf).

The paper's Methods section describes the coupled inner, central, and outer
cortical-surface evolution and their tissue-membership, force, and nesting
relationships. It is the algorithm-family context for the surface stages,
while the selected VandyCRUISE and MaCRUISE source revisions define the
specific port contract.

## Implemented-code crosswalk

| Python code | Paper relationship | Boundary |
| --- | --- | --- |
| [`tgdm.py`](../../src/pymacruise/tgdm.py): `tgdm_nested_arrays`, `tgdm_nested_images`, `tgdm_evolve`, force builders; [`plus_pipeline.py`](../../src/pymacruise/plus_pipeline.py): `run_plus_tgdm_stage`; [`cli.py`](../../src/pymacruise/cli.py): `tgdm` command | Coupled inner/central/outer CRUISE surfaces and tissue-driven nested level-set evolution, exposed through RAS-validating image and CLI adapters. | This path takes MaSeg/TOADS-derived fields and source-selected profile settings. The published CRUISE pipeline's FANTASM/AutoFill stages are not all implemented here. The MaCRUISE adapter's central-shrink pass and inner recomputation are source behavior, not attributed to the paper. |
| [`ace.py`](../../src/pymacruise/ace.py), [`gvf.py`](../../src/pymacruise/gvf.py), [`topology.py`](../../src/pymacruise/topology.py) | ACE, GVF forces, and topology correction are related CRUISE inputs/steps. | Each algorithm is mapped to its own paper and selected source in P51, P54, and P55. This citation does not make this module a complete reimplementation of those components. |
| [`surfaces.py`](../../src/pymacruise/surfaces.py), [`thickness.py`](../../src/pymacruise/thickness.py) | Han et al. describe triangulating the CRUISE level sets with Connectivity Consistent Marching Cubes (CCMC; Han et al. 2003), and give CRUISE thickness at central-surface points as the difference between the inner and outer signed-distance level sets (Eq. 13). | PyMaCRUISE substitutes scikit-image Lewiner marching cubes for JIST `IsoSurfaceOnGrid`/CCMC. Its volume level-set subtraction and selected Java closest-point surface measurement are source-mapped outputs, not asserted to reproduce the paper's central-surface Eq. 13. The JIST isosurface comparison remains open, and mesh substitution prevents claims of identical extracted surfaces or thickness results. |

## Validation scope

P29 and its synthetic reference reports provide source-comparison evidence for
selected TGDM operations and phantoms. Those reports should be described as
matching the Java harness reference for their tested cases; this paper citation
does not establish broad clinical or end-to-end equivalence.
