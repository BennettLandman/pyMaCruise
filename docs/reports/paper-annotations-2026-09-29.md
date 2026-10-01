# Published-paper code annotations — P49–P57

## Validation run

The initial P49–P57 paper-to-code annotation pass was added in
implementation/documentation commit `d936ef52e4c2744ea1a9a0227b89aa53b8a58b96`.
The full repository suite passed at that revision with
`.venv/bin/python -m pytest -q`: **365 passed in 89.50 seconds**. The matching
full-suite entry is in
[`full-pytest-2026-09-29.md`](full-pytest-2026-09-29.md). These annotation-only
changes do not alter runtime algorithms. Paper mapping explains method context
and source relationships; it does not replace the separate numerical/source
validation criteria.

## Per-tool paper and code map

| Plan row | Paper | Annotation location | Scope/limit recorded |
| --- | --- | --- | --- |
| P49 TOADS | Bazin & Pham, 2007 and 2008 | [`docs/papers/toads.md`](../papers/toads.md), TOADS module docstrings | Maps only current preprocessing, membership, atlas, registration, topology, heap, thinning, and crop components; full TOADS integration and convergence remain open. `toads_dura.py` is excluded as a separate MaCRUISE stage. |
| P50 LesionTOADS | Shiee et al., 2010 | [`docs/papers/lesion-toads.md`](../papers/lesion-toads.md), selected implemented `lesion_*.py` module docstrings | Documents descriptor priors, selected wrapper normalization/centroids, bounded growth/topology and iteration/statistics slices, and a capped distance primitive; these are source-specific components rather than the paper's complete iterative algorithm. Atlas loading, the complete workflow, and deployed plugin version remain open. |
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

## P50 documentation follow-up (2026-09-30)

The P50 crosswalk and module documentation now also connect
`lesion_growth.py` to the paper's topology-preserving competition at a
conceptual level, `lesion_fixed_input_loop.py` to iterative lesion-aware
segmentation, and `lesion_statistics_refresh.py` to statistical modeling and
channel-weight updates. The comments identify queue ordering, local topology
decisions, refresh formulas, and wrapper stop conditions as selected
VandyCRUISE source behavior. They limit evidence to the named synthetic
Java-harness slices and do not claim whole-volume topology preservation,
complete convergence, or end-to-end LesionTOADS equivalence.

This was a documentation-only follow-up; the full-suite result above predates
later implementation slices and is not evidence for those newer slices.

## P52 / P57 documentation follow-up (2026-09-30)

P52 now links the Han et al. (2004) coupled-surface context to the RAS image
entry point and MaCRUISEplus TGDM wrapper, while retaining the selected source
as authority for force profiles and wrapper-only central-shrink/recomputation
steps. P57 now names Lewiner et al. (2003) beside the concrete
`marching_cubes(method="lewiner")` operation and its plus/CLI wrappers. The
notes distinguish Lewiner's external triangulator from JIST's
`IsoSurfaceOnGrid` settings, Java near-zero perturbation, and RAS affine
conversion. The 18-connectivity difference remains an open substitution risk;
neither mesh equivalence nor downstream thickness parity is claimed.

These are documentation-only edits. Existing analytic mesh and selected TGDM
Java-harness evidence remains limited to the tests and reports named in the
crosswalk; the full-suite result above predates newer implementation slices.

## P51 documentation follow-up (2026-09-30)

The ACE/MaACE guide now distinguishes the source-derived binary atlas-label
boundary helper (`boundary.py` / `run_plus_boundary_stage`) from Huo et al.'s
paper-level ``mu_boundary`` term. The code comments trace exact boundary
neighborhood and cleanup behavior to selected V3 source; the plus ACE wrapper
accepts a precomputed boundary-like image and applies its source mode-3 merge.
No equivalence to the paper's boundary definition, max-based speed equation,
or complete MaACE method is claimed. The ACE distance/speed field remains
connected to Tosun et al.'s ACE equations, and ThinACE remains a separate
VandyCRUISE cleanup implementation rather than a paper-specified thinning
algorithm.

This documentation-only update adds no new numerical validation evidence;
source-comparison claims remain those cited in `ace-port.md` and its test
reports.

## P49 documentation follow-up (2026-09-30)

The TOADS paper crosswalk and module docstrings now cover recently added
ordinary-OBJECT statistical updates, class-pair similarity, default/atlas-type
membership compositions, optional top-k initialization and registration, the
post-registration template refresh, and the bounded continuing-iteration
control slice. The comments use Bazin and Pham (2007, 2008) only for the broad
statistical-classification, atlas-registration, and topology-aware iterative
roles. Exact formulas, candidate/class ordering, pyramid schedules, wrapper
call order, and stop conditions are attributed to the selected VandyCRUISE
source. The new annotations do not claim a complete TOADS workflow or
paper-level equivalence; synthetic Java-harness evidence remains limited to
the named implementation slices.

This was a documentation-only follow-up; the full-suite result above predates
later implementation slices and is not evidence for those newer slices.
