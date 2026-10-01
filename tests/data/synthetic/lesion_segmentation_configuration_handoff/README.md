# LesionTOADS segmentation configuration handoff

This fixture stores two-channel float32 intensity maxima and scales as exact
bit patterns. They are per-channel scalar statistics, not spatial images, so
JSON is the appropriate fixture format.

The selected wrapper branches to `setLesionWeight` or `reSetLesionWeight`,
then calls `setIntensityMax`, `setIntensityScale`, and `setGMflag` in that
order. The Java harness checks both lesion-weight branches, enabled GM flag,
array reference retention, and exact float bits. It does not validate how
preprocessing computes the intensity arrays.
