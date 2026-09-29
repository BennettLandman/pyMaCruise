# Topology-preserving TGDM paper-to-code notes

## Published method

Han X, Xu C, Prince JL. “A topology preserving level set method for geometric
deformable models.” *IEEE Transactions on Pattern Analysis and Machine
Intelligence*. 2003;25(6):755–768. DOI
[10.1109/TPAMI.2003.1201824](https://doi.org/10.1109/TPAMI.2003.1201824);
[author-hosted paper](https://iacl.ece.jhu.edu/pubs/p150j-ieee.pdf).

Section 3 presents a narrow-band level-set method that tests candidate sign
changes for simple-point topology admissibility, updates accepted changes,
limits non-simple changes, and reinitializes the distance representation and
band. Section 3.5 notes that voxel visit order can affect results. The paper
supports the algorithm rationale; the selected Java source supplies the
behavior to reproduce.

## Implemented-code crosswalk

| Python code | Paper relationship | Boundary |
| --- | --- | --- |
| [`tgdm.py`](../../src/pymacruise/tgdm.py): `tgdm_geometry`, `tgdm_force` | Geometric deformable-model geometry and force/update terms used by the level-set evolution. | The exact finite-difference stencils, force profiles, and wrapper parameters follow VandyCRUISE `GenericTGDM` and MaCRUISE adapters. |
| [`tgdm.py`](../../src/pymacruise/tgdm.py): `tgdm_evolve` | Synchronous narrow-band evolution, candidate sign-change gate, level-set reinitialization, and outer-pass rebuilds. | The selected source uses `CONNECT_18_6`; paper-listed 3-D connectivity pairs are alternatives, not a guarantee that this selected pair is the paper's universal choice. Java-harness cases establish only their reported local/synthetic scope. |
| [`tgdm.py`](../../src/pymacruise/tgdm.py): `tgdm_simple_point_618` | The local simple-point admissibility idea used when a voxel sign would change. | Its geodesic neighborhood encoding and `(18,6)` behavior are matched to the selected Java wrapper and tested against source cases; do not infer complete-volume topology or paper equivalence from the predicate alone. |

The selected MaCRUISE adapter additionally defines profile-specific outer
passes and nested-surface behavior. Those source details are documented in
[`../cruise-tgdm-port.md`](../cruise-tgdm-port.md), not attributed to the
published method.
