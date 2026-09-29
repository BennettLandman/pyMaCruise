# TOADS growing-queue initialization reference

Three deterministic 7×7×7 synthetic skeletons exercise the selected source
initGrowingLabels method. The mixed case places four different skeleton labels
around a common empty voxel; the single-seed case has one class-0 point; and
the no-boundary case fills every interior voxel with class 0. Memberships use
the established integer formula divided by 30 in float32. Ordering is an
asymmetric negative formula on the same grid.

The Java reference records the source's interior availability reset,
availability clearing for seeded candidates, every priority/payload, and the
entire Java min-heap pop order. Border availability is initially false because
the selected Java scan visits border coordinates and assumes they cannot be
available. No images or subject data are involved.

Regenerate with
scripts/generate_java_toads_growing_queue_initialization_reference.sh. This
validates queue initialization only; topology-gated growth is a later step.
