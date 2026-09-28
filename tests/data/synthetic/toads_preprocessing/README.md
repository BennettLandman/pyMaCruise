# TOADS robust-normalization phantom

This deterministic 13×11×9 RAS fixture has a rectangular foreground brain
mask, zero-valued exterior, a repeatable intensity ramp, and three bright
outliers. It validates the source `ImageFunctions.robustMinimum` /
`robustMaximum` histogram and `MultipleToadPreprocess.normalizeImages` span
division before the atlas/deformable-segmentation stages are ported.

`t1_n4.nii.gz` and `brain_mask.nii.gz` are the only committed image files.
`java_reference.json` stores expected scalar normalization results from
`tests/java_reference/ToadsIntensityReferenceHarness.java`. The generator
`make_phantom.py` recreates the input volumes; the harness source and command
are described in [docs/toads-port-plan.md](../../../../docs/toads-port-plan.md).
