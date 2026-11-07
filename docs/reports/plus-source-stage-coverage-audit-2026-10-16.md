# MaCRUISEplus source-stage coverage audit — 2026-10-16

> Historical dated record. Current testing readiness, native/Python atlas preparation and run instructions supersede stale open-work statements below; see [the October 27 readiness report](plus-testing-readiness-2026-10-27.md).

> Historical snapshot. The finding below that AtlasEM and pre-LesionTOADS
> Dura were absent from the Python runner was superseded on 2026-10-17 by the
> optional AtlasEM/Dura integration reported in
> [the integration report](plus-atlasem-dura-integration-2026-10-17.md).
> The current workflow has ten baseline stages and eleven with that optional
> stage; selected-atlas behavior and a complete independent source run remain
> unverified. See the [2026-10-18 continuation report](plus-continuation-2026-10-18.md).

## Finding

The current restartable Python PLUS workflow is a runnable nine-stage
LesionTOADS-based workflow, but it does not reproduce the full selected
`run_maCRUISEplus` source sequence. The key omitted source behavior is the
AtlasEM initialization and its Dura-removal pass before LesionTOADS. The
user-supplied mask replaces BET/SPECTRE skull stripping; it does not replace
AtlasEM or this Dura stage.

This is a source-coverage gap, not a numeric mismatch that can be solved by
loosening a tolerance. Synthetic PLUS convergence, Java comparisons of
LesionTOADS methods, and the nine-stage Python run remain valid evidence for
those narrower scopes, but they do not establish end-to-end PLUS equivalence.

## Source call chain

The selected MATLAB source call chain is:

1. `run_maCRUISEplus.m` rigidly registers raw FLAIR to T1 with
   `FlairRegistration`/NiftyReg Aladin (`rigOnly=1`).
2. `maCRUISEpluspipeline.m` prepares inputs and MaSeg-derived tissue maps.
3. `run_LesionToads_pipeline.m` N4-corrects T1 and registered FLAIR, prepares
   the SPECTRE/BET mask, strips both modalities, and clamps interpolation
   negatives to zero.
4. `run_EM.m` invokes `MedicAlgorithmAtlasEMSegmentation` on T1 only and
   returns CSF/classification plus GM, WM-fill, and WM-mask products.
5. `run_RemoveDura.m` evolves the Dura boundary from those AtlasEM products
   and returns the Dura-stripped T1 plus its mask.
6. The Dura mask is applied to registered FLAIR; the paired Dura-stripped
   T1/FLAIR then feed `run_LesionToads.m`.
7. The outer pipeline transfers/crops LesionTOADS maps, combines them with
   MaSeg tissues, runs the plus CRUISE surface workflow, and refines MaSeg.

The selected `run_EM_jist.m` overrides the plugin's built-in AtlasEM default
with `Atlas/ibsr8obj/cruise-atlas-8obj-noBG-EM.txt`. Its wrapper settings are
T1 MPRAGE, smoothing `10.0`, relative-energy threshold `0.01`, maximum 99
iterations, `dura removal inputs`, Atlas priors, rigid registration, demons
scale `2.0`, demons smoothing `1.0`, and inhomogeneity correction disabled.
The following Dura stage uses the selected MaCRUISE parameters already
recorded in [the TOADS/Dura plan](../toads-port-plan.md).

## Python coverage found

`src/pymacruise/plus_full_workflow.py` composes preprocessing, crop, atlas
alignment, registered recrop, the LesionTOADS loop, LesionTOADS Dura outputs,
native Dura transfer, MaSeg cortical separation, and the plus surface
workflow. No AtlasEM or pre-LesionTOADS Dura stage is present in this stage
list or its selected Python call chain. Current PLUS results therefore start
from a different algorithmic handoff than the original wrapper.

The Python selected-mask contract still intentionally omits BET/SPECTRE. Its
mask may be applied after source-compatible N4 and to registered FLAIR on the
T1 grid; this does not account for the separate AtlasEM/Dura initialization.
The current `run_plus_selected_mask_reference.m` retains those original
scientific stages and is intended to expose the gap in a fresh source run.

## Atlas availability

The available MIPAV and VandyCRUISE source trees contain the AtlasEM wrapper,
algorithm, preprocessing, segmentation, and `DemonToadDeformableAtlas` source.
The required selected 8-object AtlasEM descriptor and payloads are not present
under this repository's `data/atlases/`, nor in the inspected supplied source
checkouts. The external PLUS source runner checks for the descriptor in the
MIPAV plugin installation and records this resource as a required run input.
To make the Python port self-contained, the actual selected descriptor and all
referenced atlas assets must be staged under `data/atlases/` as regular Git
files, with hashes and attribution recorded. Do not substitute the bundled
12-object LesionTOADS atlas or a generic tissue classifier.

## Source provenance

The relevant supplied source files were read-only audited and SHA-256 hashed:

| Source file | SHA-256 |
|---|---|
| `MaCRUISE/pipeline/maCRUISEpluspipeline.m` | `df826582d539735274ebc76407545bd04423cb23247e337b500d53d37ccac65b` |
| `MaCRUISE/pipeline/run_LesionToads_pipeline.m` | `a2ec7e01b20df4a167147388df48a8fcf459d087358cadd6f769702902cf628f` |
| `MaCRUISE/pipeline/run_EM_jist.m` | `badc38dc6a867cfa403f459e3ce65ee07f6243820a741583a8af508b7cd2e446` |
| `MaCRUISE/pipeline/run_RemoveDura.m` | `b1df9a6852ffb42b470851a13ac60d54a743c367d562165759dd7861082345b7` |
| `MaCRUISE/pipeline/run_RemoveDura_jist.m` | `57473ab813e81b89ffed433de89b2fad8e111f6b68fa41e30c920639cd864e80` |
| `MaCRUISE/pipeline/run_LesionToads.m` | `bd26a93133c6063d4c18b62cc1abfe5fb14d2e712e6a8ce9f936881ed5452ee6` |
| `MaCRUISE/pipeline/run_LesionToads_jist.m` | `39f4318cbe8aa905e92bf5a6e3c543a312c377d0cdf685020bbe1694f3ae43fe` |
| `MedicAlgorithmAtlasEMSegmentation.java` | `96ad55161d772995135604792054c2f00d159b01159f82159a9fd06a33c29c10` |
| `AlgorithmAtlasEMSegmentation.java` | `923b3decdb4826351e3285eaa437648763d93af75e34e61fd1565eda744e54d2` |
| `AtlasEMPreprocess.java` | `acec23f916b3134e53381e83ff3de7ca569f785e6b866bfe135fb409c2926d2b` |
| `AtlasEMSegmentation.java` | `b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a` |
| `DemonToadDeformableAtlas.java` | `96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a` |

The hash table uses the source file paths supplied at the start of this report;
full absolute checkout paths are in the repository history and source-run
packet. This was a source-scope audit only. No MATLAB, Java plugin, or image
algorithm was executed for this report.
