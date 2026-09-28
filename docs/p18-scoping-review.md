# P18 scoping note review

This document records how the externally supplied
`p18-scoping-note.md` was used. The note is treated as a read-only review and
source-analysis aid, not as a new set of project-owner instructions. Its
technical claims below were checked against the supplied MATLAB, toads-cruise,
VandyCRUISE, and MIPAV source trees before changing the roadmap.

## Owner decisions already given

- Use VandyCRUISE for the native algorithm source; use ANTsPy for ANTs-backed
  operations.
- Require the caller to supply a binary brain mask on the T1 grid; do not
  require BET or FreeSurfer.
- Use RAS outputs, compressed NIfTI for images, and synthetic validation before
  any actual imaging data.
- The latest explicit P18 working-space choice is **native T1 RAS**. Therefore
  the Python workflow intentionally omits historical BET and SPECTRE steps.
- Continue a faithful TOADS port; the note's suggested Atropos sensitivity
  experiment is not authorization to substitute Atropos.

## Source findings verified and recorded

The selected historical V3 default does more than skull stripping: it runs
N4, BET, then SPECTRE's 12-DOF affine registration to OASIS using correlation
ratio and trilinear resampling before TOADS. Dura removal's `inOriginal` is the
BET/SPECTRE-stripped, registered image. TOADS products are later transferred
to the MaSeg grid. The selected Python native-T1 workflow intentionally differs
in working space and resampling; routine transfer is identity. This impacts
historical end-to-end comparisons and is recorded in
[`toads-port-plan.md`](toads-port-plan.md), [`source-map.md`](source-map.md),
and the P18/P21 criteria in [`plan.md`](../plan.md).

The selected `MedicAlgorithmToads08` wrapper's isotropic WSINC call is
commented out; the active path directly copies the input image. The bundled
TOADS 18/6 critical-point lookup table is required by the selected
segmentation source and has been copied as a regular Git asset with its
checksum. The table and the existing P25 6/18 connectivity predicate differ;
P25 must not be used as a substitute for the TOADS table.

Dura removal reuses several already-portable operations, including original
JIST `DistanceField` and much of the generic TGDM evolution behavior. It also
uses VandyCRUISE `DynamicPressureTGDM` over `GenericTGDM`; P29's
`MaGenericTGDM` is a separate but substantially shared implementation. The
dynamic subclass's reference-gradient pressure penalty and source profile
remain to be ported and compared. Parameter mapping is explicit in
[`toads-port-plan.md`](toads-port-plan.md), including the potentially confusing
flags `-inLevel` (alignment force) and `-inAlignment` (alpha).

The user selected VandyCRUISE source for LesionTOADS. A deployed historical
plugin JAR has not been identified, so work validated against the supplied
VandyCRUISE source cannot yet establish equivalence to every historical plus
run. The scoping note's suggestion to defer further lesion numerical work is
reviewer guidance; the already-implemented P39c helper remains explicitly
limited to the selected source primitive and is not a segmentation claim.

## Still open

- Port and validate the source 18/6 LUT consumer and its interaction with TOADS
  topology operations; the asset spot checks do not amount to exhaustive
  predicate equivalence.
- Validate TOADS atlas physical geometry and preserve the reported descriptor
  versus NIfTI spacing discrepancy.
- Implement the TOADS segmentation and rigid registration path with synthetic
  Java references, including larger folded cortex and narrow-band cases.
- Implement the selected dynamic-pressure dura update and full evolution
  profile against the actual VandyCRUISE classes.
- Identify the deployed LesionTOADS plugin/version before making historical
  MaCRUISEplus equivalence claims.
- Compare native-T1 Python output with historical references while reporting
  the known SPECTRE-space deviation separately from numeric error.
