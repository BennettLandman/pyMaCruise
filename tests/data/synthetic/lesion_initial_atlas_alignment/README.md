# LesionTOADS initial atlas alignment fixture

This deterministic synthetic case supplies a padded, two-channel 33³ image,
three distinct 35×37×36 atlas shape maps, an explicit preprocessing-seeded
rigid transform, and crop bounds. All NIfTI files use an explicit RAS affine.
The image has asymmetric folded intensities, a zero-distance initial class,
and a first-channel centroid tie. The maps vary by class and coordinate so the
rigid pyramid has actual atlas support. `inputs.json` records the fixed input
contract; `java_reference.json` stores exact initial-membership float bits and
the wrapper call-order/configuration excerpt.

Regenerate image fixtures with `generate_fixture.py` and the reference excerpt
with `scripts/generate_lesion_initial_atlas_alignment_reference.sh`.

Source anchors are pinned in
`tests/java_reference/LesionInitialAtlasAlignmentReferenceHarness.java`:
selected VandyCRUISE `AlgorithmLesionToads.java` SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`,
`LesionToadAtlas.java` SHA-256
`718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e`, and
`ParametricTransform.java` SHA-256
`22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`.
The harness validates wrapper setup order and transform seeding. The numerical
coarse-to-fine optimizer remains the independently validated P39o component;
this fixture does not claim a second Java optimizer comparison or complete
LesionTOADS initialization.
