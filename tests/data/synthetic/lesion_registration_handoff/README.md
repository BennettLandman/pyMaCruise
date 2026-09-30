# P39zzzzt LesionTOADS membership-to-registration handoff

`registration_handoff_ras.nii.gz` is a synthetic 6×5×4 RAS fixture. Frames
0–3 are P39zzzzs Java-reference class memberships with fields disabled, frames
4–7 are the corresponding field-enabled membership output, and frames 8–11
are the atlas shape maps. The fixture is derived from the pinned P39zzzzs
`java_reference.json`; the test reruns the Python P39zzzzs bridge and requires
both membership arrays to match these frames bitwise before registration.
Its affine has positive RAS axes and anisotropic 0.9, 1.2, and 1.5 mm spacing.
Class 1 is folded/asymmetric, class 2 includes support on the image boundary,
and class 3 is empty to exercise the Java empty-bounds sentinels.

`LesionRegistrationHandoffReferenceHarness.java` applies the selected Java
registration coefficient, line-search, scale-1 `registerAllShapes` loop,
rotation refresh, and `computeTransformedShapeBoundingBox` excerpt to the
P39zzzzs Java membership outputs. The harness pins selected
`LesionToadAtlas.java` SHA-256
`718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e` and
`AlgorithmLesionToads.java` SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`.
Rebuild the RAS fixture and Java reference with
`scripts/generate_java_lesion_registration_handoff_reference.sh`.

The selected wrapper uses one rigid transform, a scale-1 update, `subsample=2`,
`itMax=initAlignIter` (default 10), lambda 1, factor 1.5, strict energy
improvement, `minEdiff=1e-6`, and `minLambda=0.001`. Tests compare every trial,
iteration state, final transform/rotation and class bound against Java and
verify that no membership or shape inputs are mutated. This row covers only
the fixed-input active registration handoff; initialization pyramid,
multi-transform mode, topology, and outer-loop convergence remain excluded.
