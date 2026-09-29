# Published-paper code annotations — P49–P57

## Validation run

The nine paper-to-code rows below were completed in implementation/documentation
commit `d936ef52e4c2744ea1a9a0227b89aa53b8a58b96`. The full repository suite
passed with `.venv/bin/python -m pytest -q`: **365 passed in 89.50 seconds**.
The matching full-suite entry is in
[`full-pytest-2026-09-29.md`](full-pytest-2026-09-29.md). These annotation-only
changes do not alter runtime algorithms. Paper mapping explains method context
and source relationships; it does not replace the separate numerical/source
validation criteria.

## Per-tool paper and code map

| Plan row | Paper | Annotation location | Scope/limit recorded |
| --- | --- | --- | --- |
| P49 TOADS | Bazin & Pham, 2007 and 2008 | [`docs/papers/toads.md`](../papers/toads.md), TOADS module docstrings | Maps only current preprocessing, membership, atlas, registration, topology, heap, thinning, and crop components; full TOADS integration and convergence remain open. `toads_dura.py` is excluded as a separate MaCRUISE stage. |
| P50 LesionTOADS | Shiee et al., 2010 | [`docs/papers/lesion-toads.md`](../papers/lesion-toads.md), three implemented `lesion_*.py` module docstrings | Documents descriptor priors, normalization/centroids, and a capped distance primitive; the iterative lesion engine, atlas loading, and full path remain incomplete, and deployed plugin version is unresolved. |
| P51 ACE / MaACE | Tosun et al., 2004; Huo et al., 2016 §2.3.1 | [`docs/papers/ace-maace.md`](../papers/ace-maace.md), `ace.py`, `skeleton.py` | Separates the ACE field, MaACE multi-atlas boundary prior, and source-specific ThinACE cleanup. The selected Python path does not yet supply Huo's distinct boundary-label prior and does not claim complete MaACE. |
| P52 CRUISE | Han et al., 2004 | [`docs/papers/cruise.md`](../papers/cruise.md), `tgdm.py` | Maps nested surface evolution and identifies MaSeg/TOADS input and adapter differences from the published FANTASM/AutoFill workflow. Mesh/thickness substitutions and equivalence limits remain explicit. |
| P53 TGDM | Han, Xu & Prince, 2003 §3 | [`docs/papers/tgdm.md`](../papers/tgdm.md), `tgdm.py` | Maps narrow-band evolution, simple-point checks, and reinitialization while deferring exact `(18,6)` connectivity and profiles to source. It does not infer volume-wide topology/equivalence from local tests. |
| P54 Fast-marching topology correction | Bazin & Pham, 2007 | [`docs/papers/topology-fast-marching.md`](../papers/topology-fast-marching.md), `topology.py` | Maps scalar-field propagation and local component checks to selected VandyCRUISE methods. Explicitly distinguishes this from Han et al.'s separate graph-based GTCA, which is not implemented here. |
| P55 Generalized GVF | Xu & Prince, 1998 | [`docs/papers/gvf.md`](../papers/gvf.md), `gvf.py` | Explains the GGVF model and its mapping; documents SciPy conjugate gradients as a replacement for the selected source's multigrid solve and records source-specific stencils/weights. |
| P56 N4 | Tustison et al., 2010 | [`docs/papers/n4.md`](../papers/n4.md), `n4.py` | Attributes the bias-estimation implementation to ANTsPy/ANTs; distinguishes selected wrapper/CLI-equivalent settings from paper parameters and limits synthetic validation claims. |
| P57 Surface extraction | Lewiner et al., 2003; scikit-image library paper | [`docs/papers/marching-cubes.md`](../papers/marching-cubes.md), `surfaces.py` | Documents the external Lewiner triangulator, near-zero source perturbation and RAS affine map; JIST `IsoSurfaceOnGrid` 18-connectivity comparison remains an open substitution risk. |

The complete set of bibliography and mapping notes is indexed in
[`docs/papers/README.md`](../papers/README.md). Each note links to its
publication record or accessible manuscript and names the mapped code.
