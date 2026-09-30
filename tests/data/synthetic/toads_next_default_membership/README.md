# P18czj synthetic continuation case

This fixture advances one actual P18czh registered handoff through the
non-stopping P18czi branch, then runs source-ordered thinning and the next
default 11-class membership pass.

`subject_images.nii.gz` is a 25×25×25, two-channel, asymmetric folded phantom
in true RAS. It is built by periodic extension of the smaller source-faithful
folded synthetic volume so that the selected source's 100-voxel minimum
membership threshold is reached without changing that threshold. It is a
topology/queue stress phantom, not clinical data. The template atlas is a
separate 11×10×9×11 true-RAS grid in `template_atlas_shapes.nii.gz`.

The saved continuing membership, transformed prior, grown segmentation,
refreshed closest-label map, and relation adjacency are fixed handoff inputs
for the Java excerpt. The generator also executes the current Python
P18czh/P18czi adapters and then compiles the selected-source excerpt against
the pinned VandyCRUISE `BinaryHeap4D`, `CriticalPointLUT`, and `PowerTable`.

The Java and Python queue traces match through 8,962 pops, including 978
segmentation-label updates, skeleton/availability changes, and final class
volumes. The following 11-class maps, denominator, and objective match Java
float bits. Refreshed closest-label input is separately shown to change the
membership output relative to the stale pre-registration map.

Run `PYTHONPATH=tests ./.venv/bin/python
tests/data/synthetic/toads_next_default_membership/generate_reference.py` from
the repository root to regenerate the committed NIfTI snapshot and reference.
The existing P18czh/P18czi tests continue to validate their own registration
and wrapper branches; this case does not estimate class statistics or repeat
registration.
