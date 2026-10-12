# P41z ribbon-support loss localization — 2026-10-11

## Finding

The held-fixed P41z iteration-11 crop loses much of its WM support during the
selected TGDM initialization and inner-surface evolution. The earliest
measured drop is the six-pass initial smoothing; the larger loss occurs in
`initial_inner`. Both outputs match selected Java TGDM on the same serialized
inputs, including surface sign and distance gates. This localizes the present
failure to the source-sequence response to P41z's seed and tissue fields, not
to a Python-vs-Java TGDM discrepancy. It does not establish which upstream
input or anatomical property causes that response.

| P41z stage | WM label 44 inside (≤0) | GM label 100 inside (≤0) | Negative voxels | 26-connected components |
| --- | ---: | ---: | ---: | ---: |
| Initial level set | 99.379% | 1.374% | 54,657 | 1 |
| Smoothed initial | 86.636% | 4.102% | 51,819 | 5 |
| Initial inner | 59.747% | 0.102% | 22,681 | 12 |

The saved crop is 73×91×63 at approximately 2.01 mm spacing. The full
iteration-11 Java comparison reports zero sign changes for `smoothed_initial`
and `initial_inner`; smoothed-stage Dice is 1.0, symmetric mean distance is
`8.46e-8 mm`, and sampled Hausdorff is `1.75e-5 mm`. The full metrics, hashes,
label counts, Euler values, and connectivity counts are in the
[machine-readable localization](p41z-tgdm-support-loss-localization-2026-10-11.json)
and the [iteration-11 Java comparison](p41z-tgdm-java-iteration11-2026-10-11.md).

## Controlled-case context

On the 2 mm decimation of the 1 mm folded-cortex phantom, selected Python TGDM
keeps 100% of WM inside after smoothing and after `initial_inner`; GM support
is 87.64% inside after smoothing and 0% inside after `initial_inner`. That
case uses WM iso 0.6 and a simpler analytic shell. P41z uses the plus entry
point's WM iso 0.5 and atlas-derived tissue/seed fields. Thus the control
shows that 2 mm spacing is not by itself sufficient to force this support
loss, but it does not isolate P41z's anatomy, seed, or parameter interaction.
The level-set thickness sensitivity at 2 mm is documented separately in the
[resolution report](tgdm-resolution-sensitivity-2026-10-11.md).

## Next diagnostic step

Compare the topology-refined P41z seed and its MaSeg-derived WM field before
TGDM, including their physical boundaries and label overlap, against a
full-resolution source-selected case. Preserve the selected 0.5 WM iso and
P41z source inputs in that comparison. Do not change TGDM forces or weaken the
ribbon acceptance thresholds to make this checkpoint pass.

This is a synthetic held-fixed stage comparison. It is not an independently
generated Java MaCRUISEplus workflow or real-MRI validation.
