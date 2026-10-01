# Default TOADS continuation to next thinning pass

This synthetic 9×9×9 fixture chains the selected Java growth harness output
through relation-list construction and the next initThinningLabels /
propagateThinning pass. It uses four classes, the deterministic folded
phantom and exact membership values already used by the continuation reference,
zero transformed-shape values, volumeRatio=0.5, and limit=3.

The selected short constructor uses exact memberships. Atlas registration is
held fixed and remains outside this seam. The first next-loop operation after
the non-stopping growth/relation branch is thinning. This fixture uses the
selected 100-voxel minimum-volume floor. One class reaches that floor during
processing, and the resulting segmentation differs from the post-growth
labels. The chained reference therefore exercises queue initialization, the
volume guard, topology checks, and the source stop limit alongside the state
handoff.

The generated reference records the continuation grown labels, refreshed
class adjacency, initialized volumes/boundaries, processed queue trace,
skeleton/availability, and final state. Regenerate it with
sh scripts/generate_java_toads_default_continue_thinning_reference.sh.
This validates one post-growth next-thinning transition only; registration,
shape transformation, subsequent memberships, and full TOADS are out of scope.
