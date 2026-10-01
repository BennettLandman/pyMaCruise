# P39zzzzx two-pass loop fixture

This bounded loop reuses the synthetic asymmetric true-RAS multichannel input
at `../lesion_continuing_iteration/iteration_inputs_ras.nii.gz`. It begins
from a continuing P39w-style branch and runs exactly two P39v/P39w passes,
using fixed spatial and registration inputs. The selected Java wrapper
excerpt is pinned to VandyCRUISE `AlgorithmLesionToads.java` SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`.

The fixture uses positive `n_iter_max=2` and a zero convergence threshold:
pass one continues and refreshes relations; pass two stops at the signed
iteration ceiling and runs final growth/classification. The reference records
the selected wrapper's energy-ratio bits and per-pass call order. Membership,
registration, growth, and relation arithmetic retain their individual
selected-source Java harness evidence; the loop harness validates controls
and stage ordering, not a new Java execution of all numerical kernels.

`loop_states_ras.nii.gz` is generated from the two actual Python loop-pass
inputs and pins their starting segmentation, four class memberships, lesion
membership, and black-hole membership (seven channels per pass). The Java
wrapper harness reads this NIfTI directly and independently recomputes the
selected `initThinningLabels()` class volumes in x/y/z scan order. This checks
that class volume is reset from the current segmentation/memberships on each
pass; the previous growth volume is carried as state but is not reused by
thinning initialization. The selected source adds lesion/black-hole auxiliary
memberships only for the configured `lesionClass` index, even when another
class is named `BLACKHOLE`.

The Java loop excerpt also includes convergence-stop cases with a threshold
above the first pinned ratio for both positive and negative signed iteration
ceilings. Positive ceilings run final growth; negative ceilings skip it. The
corresponding class-volume reference bits are recorded in
`java_reference.json`.

Regenerate the pass-state NIfTI and Java wrapper/volume reference together
with:

    sh tests/java_reference/generate_lesion_repeated_iteration_reference.sh
