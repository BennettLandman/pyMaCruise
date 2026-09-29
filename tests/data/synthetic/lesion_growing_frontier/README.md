# LesionTOADS growing-frontier synthetic reference

This 7³ synthetic fixture checks only the growing initialization selected
from VandyCRUISE `LesionToadSegmentation`: interior availability reset,
six-face frontier detection, lesion-aware best-class/ratio calculations, and
the exact source `BinaryHeap4D` pop order. The Java harness transcribes
`initGrowingLabels`, `bestMembership`, and `maxMembershipRatio`, while using
the committed selected-source `Numerics` and `BinaryHeap4D` classes.

Cases exercise additive and max lesion scoring, strict-tie/zero-score
behavior, and optimized lesion/black-hole object types. The ordinary-source
`bestMembership` running-maximum quirk is intentionally retained. This is a
formula/harness match, not a complete `LesionToadSegmentation` invocation;
no queue events or lesion segmentation are claimed.

Regenerate the JSON reference with
`scripts/generate_java_lesion_growing_frontier_reference.sh`.
