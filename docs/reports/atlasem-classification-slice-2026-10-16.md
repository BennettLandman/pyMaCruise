# AtlasEM classification output slice — 2026-10-16

`generate_atlas_em_classification` reproduces `AtlasEMSegmentation.generateClassification` from the pinned class (SHA-256 `b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a`). The selected output scans class memberships in order, replaces the selected class only on a strictly greater positive membership, maps through the Java signed-byte template label table, and writes label zero when all class memberships are nonpositive. This is distinct from the segmentation-index array used by `generateDuraRemovalOutputs`.

The Java harness and deterministic 3×2×2 phantom compare every output byte. The phantom includes a tied maximum and an all-zero voxel. Focused validation is included in the 14-test AtlasEM/morphology group recorded in P44b and the Gaussian membership report. This is an output-method slice, not evidence of atlas registration, source-configured AtlasEM, or full PLUS workflow equivalence.
