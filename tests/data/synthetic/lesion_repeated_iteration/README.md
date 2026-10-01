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
