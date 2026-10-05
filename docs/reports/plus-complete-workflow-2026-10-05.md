# Raw-input MaCRUISEplus synthetic workflow checkpoint

- Tested implementation commit: `3051b821ed93163f750b4e2b81fb98b0bad22a52`
- Focused command: `./.venv/bin/python -m pytest -q --tb=line tests/test_plus_full_workflow.py`
- Focused result: **1 passed, 0 failed, 0 skipped in 97.04s**
- Full-suite command and result: see [full suite report](full-pytest-2026-10-05-plus-complete-workflow.md)

The committed 25³ folded-cortex synthetic T1, shifted-grid FLAIR, user mask, MaSeg labels, and mini atlas run through the nine-stage `plus-full-run` CLI: preprocessing, crop, initial alignment, registered recrop, LesionTOADS, Dura outputs, native-grid Dura transfer, MaSeg cortical separation, and the plus surface workflow. The test verifies top-level stage order and child hashes, terminal classification/subtype, Dura channel values, native-grid transfer against the existing physical-affine adapter, RAS level-set/thickness images, VTK surface products, thickness JSON, resume, and corrupted-tail repair.

Surface-thickness statistics are valid on this phantom (mean **0.0639 mm** over **376** nonzero vertices). The volume-thickness calculation has **0 crossing voxels**, so its summary statistics are JSON `null` with `available: false`; it is not reported as a measured zero. The test establishes that these unavailable statistics remain serializable and explicit.

This is synthetic workflow-composition evidence. It does not establish whole-workflow Java equivalence, equivalence of ANTsPy registration or arbitrary-grid orientation transfer to NiftyReg/JIST, or real-MRI validation. The P41t connection maps **zero additional selected Java method lines** into the existing 13,651-line TOADS denominator. The 25³ phantom is suitable for integration behavior but does not replace brain-resolution accuracy, convergence, or source-comparison work.
