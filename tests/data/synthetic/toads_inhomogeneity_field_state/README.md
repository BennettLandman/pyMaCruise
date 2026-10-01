# TOADS inhomogeneity field attachment fixture

This fixture validates only the selected `MultipleToadSegmentation.addInhomogeneityCorrection(float[][][], int)` setter. It does not exercise field fitting, corrected-image multiplication, or segmentation iterations.

The 2×2×2 field is stored as `field.nii.gz` with a positive diagonal RAS affine. The selected setter assigns the supplied array object directly to `field[channel]`, then sets `useField=true`; it neither copies nor transforms field values. Therefore this check covers reference identity, mutation visibility after attachment, per-channel replacement, unchanged other channels, null acceptance, and the state flag.

Source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`, method `addInhomogeneityCorrection`, source SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`. The selected caller is `AlgorithmMultipleToads.java` line 677, SHA-256 `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`.

The standalone Java harness is `tests/java_reference/ToadsInhomogeneityFieldAttachmentReferenceHarness.java`. This is a method-state comparison only; it does not establish end-to-end workflow equivalence.
