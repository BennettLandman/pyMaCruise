# MaCRUISE second-pass dura removal: paper-to-code notes

## Primary publication for this stage

Landman BA, Bogovic JA, Carass A, Chen M, Roy S, Shiee N, Yang Z, Kishore B,
Pham D, Bazin P-L, Resnick SM, Prince JL. “System for Integrated
Neuroimaging Analysis and Processing of Structure.” *Neuroinformatics*.
2013;11:91–103. Published online 30 August 2012.
[doi:10.1007/s12021-012-9159-9](https://doi.org/10.1007/s12021-012-9159-9) ·
[publisher record](https://link.springer.com/article/10.1007/s12021-012-9159-9) ·
[open full text (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3511612/) ·
[author-hosted article PDF](https://iacl.ece.jhu.edu/proceedings/iacl/2012/LanxNI12c-SINAPS.pdf).

This paper is the identifiable primary publication for the specific
second-pass dura-removal procedure associated with the JIST module “Remove
Dura 2009” (the paper's software listing identifies version 1.6 RC). Its
introduction describes the dura stage as previously unpublished; the paper
then gives a method-level description in its “Dura Removal” section (printed
page 95). It explains that the stage follows initial tissue
segmentation and SPECTRE brain extraction, initializes a topology-preserving
geometric deformable model from the corrected white-matter membership, evolves
between corrected WM and the WM+GM union, adds cortical-thickness
regularization, uses fast-marching topology correction to prune thin or
irregular dura connections, and applies a resulting pial-level-set-derived
mask to the original image.

The article does not document every parameter, intermediate, arithmetic
choice, or code revision of `MedicAlgorithmRemoveDura09`. It is therefore the
primary published account of the *procedure*, while the selected VandyCRUISE
and JIST source checkout is the authority for the exact wrapper behavior being
ported. The paper's named software release does not by itself establish that
every detail in the selected source revision is identical to that release.
The component papers below explain algorithms used by the stage but should not
be mistaken for standalone publications of the complete dura procedure.

## Component method papers

- Han X, Xu C, Prince JL. “A topology preserving level set method for
  geometric deformable models.” *IEEE Transactions on Pattern Analysis and
  Machine Intelligence*. 2003;25(6):755–768.
  [doi:10.1109/TPAMI.2003.1201824](https://doi.org/10.1109/TPAMI.2003.1201824).
  This is the cited TGDM foundation; it does not specify the complete
  dura-removal stage or its source-specific thickness force and wrapper.
- Bazin P-L, Pham DL. “Topology correction of segmented medical images using
  a fast marching algorithm.” *Computer Methods and Programs in Biomedicine*.
  2007;88(2):182–190.
  [doi:10.1016/j.cmpb.2007.08.006](https://doi.org/10.1016/j.cmpb.2007.08.006).
  This is the fast-marching topology-correction method cited for pruning
  dura-like attachments; it is distinct from the complete Dura wrapper.

## Source-to-code crosswalk

| Operation described by Landman et al. | Selected source / PyMaCRUISE code | Scope boundary |
| --- | --- | --- |
| Start from post-segmentation white/gray membership images and a topologically correct WM initialization | `MedicAlgorithmRemoveDura09.execute()`; [`toads_dura.py`](../../src/pymacruise/toads_dura.py), `run_dura_removal`, `dura_correct_white_matter_topology`, `dura_correct_gray_matter_topology` | Python follows the selected wrapper order and its connectivity/settings. It does not provide BET, SPECTRE, or upstream TOADS as part of this operation; PyMaCRUISE takes the selected T1-grid tissue inputs and brain mask explicitly. |
| Evolve the pial level set between WM and WM+GM | `MedicAlgorithmRemoveDura09.execute()`, `DynamicPressureTGDM`; `toads_dura.py`, `dura_initial_level_sets`, `dura_tgdm_evolve` | The paper establishes the role of TGDM and the tissue bounds. Exact profiles, update order, narrow band, bounds, and stopping behavior are selected-source details. See [`tgdm.md`](tgdm.md) for TGDM's separate method/source mapping. |
| Apply thickness regularization and suppress dura attached by irregular/thin bridges | `MedicAlgorithmRemoveDura09`'s dynamic-pressure/profile and topology-correction calls; `toads_dura.py`, `dynamic_pressure_tgdm_point_update`, `dura_tgdm_evolve`, and the topology adapters | The published section describes this effect but not the full Java equations/parameterization. Attribute the specific force and iteration settings to the source unless a line-level paper correspondence is established. |
| Build final brain mask and apply it to original T1 | `MedicAlgorithmRemoveDura09.execute()`; `toads_dura.py`, `apply_dura_output_mask`, `run_dura_removal` | Output masking and threshold order follow the selected implementation. They are not a substitute for initial brain extraction. |
| Signed-distance initialization and bounded distances | JIST `DistanceField`; `toads_dura.py`, `dura_signed_distance_field` | The publication discusses the broader evolution, while six-neighbor traversal, float32 storage, voxel-unit distances/caps, and Java boundary conventions come from JIST/source comparison. |

## Validation boundaries

Existing validation is source-slice validation on committed synthetic
phantoms, not a new test of the publication's reported performance and not
clinical-data validation. The method-level distance, dynamic-pressure point
update, topology correction, TGDM evolution, and wrapper-composition evidence
is recorded in the P18 Dura reports, including:

- [`full-pytest-2026-09-28.md`](../reports/full-pytest-2026-09-28.md) for the
  early JIST distance and dynamic-pressure method slices;
- [`full-pytest-2026-09-29-dura-removal-composition.md`](../reports/full-pytest-2026-09-29-dura-removal-composition.md)
  for the selected composed wrapper on synthetic RAS inputs;
- [`full-pytest-2026-09-30-p18-dura-tgdm.md`](../reports/full-pytest-2026-09-30-p18-dura-tgdm.md)
  for the later folded-cortex synthetic TGDM case.

Those reports retain their individual measured harness scopes. A paper
annotation is not evidence of equivalence, and synthetic wrapper agreement
does not validate the Dura method on actual imaging data.
