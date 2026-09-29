# Full pytest evidence — 2026-09-29 LesionTOADS growing event

| Field | Value |
|---|---|
| Command | .venv/bin/pytest from the repository root |
| Result | **425 passed in 54.31 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | b607fa74457e476da5337e4b7318a31dad2890f4 |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39r processes one LesionTOADS growing queue item. Five committed synthetic
Java reference cases cover lesion volume updates in additive and max modes,
non-lesion volume updates, a topology-rejected insertion, and a stale item.
Python matches every output label and availability value, float32 ordering
and volume bit, inserted heap tuple, and post-event heap order. The focused
P39r tests report **2 passed**.

The harness composes selected-method formulas with the previously validated
18/6 local topology method and uses the source CriticalPointLUT, Numerics.max,
and BinaryHeap4D implementations. It does not instantiate the whole
LesionTOADS class and validates one queue event only, not full propagation or
segmentation. Fixture details are in the committed lesion-growing-event
synthetic fixture README.
