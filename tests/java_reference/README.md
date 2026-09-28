# Java algorithm reference harnesses

## ThinACE reference

`ThinACEReferenceHarness.java` runs the copied VandyCRUISE `ThinACE.java`
implementation on committed analytic masks. Regenerate the line, cube,
nested GM/WM, and symmetric sphere `.nii.gz` source/output fixtures with
`python scripts/generate_java_thinace_reference.py`. The harness uses minimal
test-only JIST image and voxel stubs under `stubs/`; it does not run a MIPAV or
JIST application. Source attribution and provenance are documented in the
fixture README and `license_info.md`.

The selected source applies `DistanceField.solve(4.0)`; the reference harness
uses the supplied MaCRUISE `DistanceFieldyh.java`, whose source header states
that it is identical to `DistanceField`. Temporary binary interchange files
and compiled classes are written under ignored `build/` and are not committed.

## Generalized GVF reference

`GVFReferenceHarness.java` runs the selected VandyCRUISE generalized GVF
implementation using source copies under `upstream/.../algorithms/gvf/`.
Regenerate the analytic ramp, asymmetric step, and Gaussian blob reference
vectors with `python scripts/generate_java_gvf_reference.py`. The selected
GVF, gradient, Poisson, multigrid transfer, and red/black solver sources retain
their upstream comments and CVS revision 1.1 metadata. Minimal JIST/framework
stubs exercise only the required in-memory image operations; they avoid
requiring MIPAV or a JIST application. Vector arrays retain the source
component order x/y/z, used by VandyCRUISE TGDM's corresponding x/y/z gradient
terms. Comparison metrics and the voxel-axis/physical-space convention are
documented in [`docs/gvf-port.md`](../../docs/gvf-port.md).

This test-only harness runs the selected Java MaCRUISE TGDM adapter on the
committed, deterministic 25³ synthetic phantom and writes small compressed
NIfTI (`.nii.gz`) golden outputs. It requires a JDK but does not run MIPAV or
JIST as applications. Regenerate the fixtures from the repository root with
`scripts/generate_java_tgdm_reference.sh`.

## Source used

`upstream/` contains source copies from the supplied checkouts, with their
original code comments and author/version metadata retained:

- MaCRUISE `MaGenericTGDM.java`, CVS revision 1.5, from the supplied
  `toads-cruise` tree.
- `DistanceFieldyh.java`, CVS revision 1.1, from the same tree. Its header
  describes it as identical to the omitted `DistanceField` class.
- VandyCRUISE `ProfileTGDM`, inner/central/outer/smoothing profiles, and the
  6-, 18-, and 26-connectivity implementations from the supplied source tree.
- JIST `BinaryMinHeap` and `Indexable`, including their original attribution
  comments, from the supplied JIST source tree.

The test-only `stubs/` provide the numerical subset of JIST image storage,
progress reporting, mesh return types, and unused framework methods needed to
compile those sources without bringing in the GUI/application dependency tree.
The reference uses the original TGDM solver, profile force equations, topology
rules, and distance-field algorithm. The harness maps the MaCRUISE adapter's
documented call sequence to those methods; it does not call mesh extraction,
which is tracked separately in P30.

These Java files are for regression comparison only and are not imported by
the Python package. Attribution and notice review is tracked in
[`license_info.md`](../../license_info.md); the supplied MIPAV/JIST notices are
preserved under `licenses/upstream/`.

## Phantom and outputs

The phantom is generated analytically in `TGDMReferenceHarness.java`: a radius
8 membership sphere, radius 5 white-matter sphere, synthetic gray matter,
zero GVF, and a planar thin skeleton. It has a 1 mm identity RAS affine. There
is no subject data. `tests/data/synthetic/tgdm/java_reference/` contains the
initial signed seed, smoothed seed, tissue-sign field, initial inner surface,
central surface, central shrink, final inner surface, and outer surface.
