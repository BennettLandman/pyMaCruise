# LesionTOADS BLACKHOLE smoothing sweep fixture

The synthetic `5×4×6` float32 fixture contains two OBJECT classes and one
BLACKHOLE class, plus separate lesion and black-hole membership/data-term
maps. It exercises source six-face sums (including BLACKHOLE self-neighbors),
mixed direct/prior relations, a zero-shape gate, in-place x/y/z traversal, and
denominator/normalization order. A selected-Java formula harness stores all
interior class and auxiliary membership values as float32 bit patterns.

Atlas-prior additions, distance modulation, other class sets, registration,
topology, and complete segmentation are excluded. Regenerate with
`scripts/generate_java_lesion_toads_blackhole_sweep_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
