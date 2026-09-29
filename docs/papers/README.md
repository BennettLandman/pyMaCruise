# Published-method code annotations

This folder contains the paper citations, method summaries, and bounded
paper-to-code maps requested in closeout steps P48–P57. Papers explain the
algorithmic context; the selected MaCRUISE/VandyCRUISE source remains the
authority for exact version-specific behavior. A citation or comment is not
evidence of numerical equivalence.

| Plan row | Tool or method | Paper-to-code guide |
| --- | --- | --- |
| P48 | MaCRUISE / CCSE | [`../reports/ma-cruise-paper-annotation-2026-09-29.md`](../reports/ma-cruise-paper-annotation-2026-09-29.md), [`segmentation-refinement.md`](../segmentation-refinement.md) |
| P49 | TOADS | [`toads.md`](toads.md) |
| P50 | LesionTOADS | [`lesion-toads.md`](lesion-toads.md) |
| P51 | ACE / MaACE | [`ace-maace.md`](ace-maace.md) |
| P52 | CRUISE | [`cruise.md`](cruise.md) |
| P53 | Topology-preserving TGDM | [`tgdm.md`](tgdm.md) |
| P54 | Fast-marching topology correction | [`topology-fast-marching.md`](topology-fast-marching.md) |
| P55 | Generalized gradient vector flow | [`gvf.md`](gvf.md) |
| P56 | N4 bias correction (ANTsPy dependency) | [`n4.md`](n4.md) |
| P57 | Lewiner marching cubes (scikit-image) | [`marching-cubes.md`](marching-cubes.md) |

Each guide identifies the relevant publication, code locations, what the
implementation actually does, and boundaries between paper concepts and
source/library-specific choices. The guides do not replace the source
comparison reports referenced by completed plan rows.
