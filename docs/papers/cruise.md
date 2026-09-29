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
| [`tgdm.py`](../../src/pymacruise/tgdm.py): `tgdm_nested_arrays`, `tgdm_evolve`, force builders | Coupled inner/central/outer CRUISE surfaces and tissue-driven nested level-set evolution. | This path takes MaSeg/TOADS-derived fields and source-selected profile settings. The published CRUISE pipeline's FANTASM/AutoFill stages are not all implemented here. The MaCRUISE adapter's central-shrink pass and inner recomputation are source behavior, not attributed to the paper. |
| [`ace.py`](../../src/pymacruise/ace.py), [`gvf.py`](../../src/pymacruise/gvf.py), [`topology.py`](../../src/pymacruise/topology.py) | ACE, GVF forces, and topology correction are related CRUISE inputs/steps. | Each algorithm is mapped to its own paper and selected source in P51, P54, and P55. This citation does not make this module a complete reimplementation of those components. |
| [`surfaces.py`](../../src/pymacruise/surfaces.py), [`thickness.py`](../../src/pymacruise/thickness.py) | The paper describes topology-aware surface extraction and thickness calculation. | PyMaCRUISE uses scikit-image Lewiner marching cubes and its tested distance-based thickness operations. The JIST `IsoSurfaceOnGrid` comparison remains open; paper context is not evidence of identical extraction or thickness outputs. |

## Validation scope

P29 and its synthetic reference reports provide source-comparison evidence for
selected TGDM operations and phantoms. Those reports should be described as
matching the Java harness reference for their tested cases; this paper citation
does not establish broad clinical or end-to-end equivalence.
