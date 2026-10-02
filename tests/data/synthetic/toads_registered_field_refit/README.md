# Registered TOADS field-fit integration fixture

A deterministic default 11-class atlas and asymmetric subject T1/brain-mask pair, all stored as RAS `.nii.gz` images. The atlas is compact; the subject has an anisotropic affine and synthetic folded intensity structure. This fixture exercises atlas registration, crop/constructor state, and the degree-three correction fit; it is not a clinical image or whole-TOADS reference. Regenerate with `python generate_fixture.py`.
