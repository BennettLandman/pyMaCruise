# AtlasEM transformed-prior normalization slice — 2026-10-16

`normalize_atlas_em_priors` reproduces the selected masked normalization loop in `AtlasEMSegmentation.computeAtlas` after transformed atlas maps have been generated. The source class SHA-256 is `b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a`.

On a deterministic 3×2×2, three-class fixture, all normalized prior float32 values match the Java excerpt bit-for-bit. The test also verifies that zero-sum priors remain unchanged and values outside the mask are left untouched. This validates only post-transform class normalization; it does not port or validate `DemonToadDeformableAtlas` resampling/registration or the missing atlas bundle.
