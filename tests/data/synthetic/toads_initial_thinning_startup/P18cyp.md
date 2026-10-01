# P18cyp rigid registration and wrapper control slice

This slice consumes P18cyo's exact 25×25×25 RAS membership map and objective,
which are checked against `continuation_java_reference.json` before use. It
registers the 11-class map against the folded synthetic shape field with the
selected wrapper defaults: rigid transform, class 0 excluded, offset 0,
subsample 2, up to five updates, minimum energy change `1e-6`, and minimum
lambda `1e-3`. It then applies the selected wrapper's float32 energy-ratio and
stop/iteration control kernel.

No membership pass is recomputed. The source's post-thinning
`computeCorrectionField(deg)` refit is deliberately not ported in this slice;
P18cyo correction fields remain fixed.

## Java source and evidence boundary

The selected call order is pinned to the supplied VandyCRUISE sources:

* `AlgorithmMultipleToads.java` (`b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`), lines 766–810: exact membership/objective, energy ratio, optional `atlas.registerShapes()`, stop/iteration decisions, growth and classification/relation branch.
* `MultipleToadDeformableAtlas.java` (`e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`), lines 2718–2742: `registerShapes()` dispatch and `registerAllShapes()` reset/iteration/finalization.

The exact Java registration run could not be instantiated in this repository's
lightweight harness setup. I attempted to compile the selected
`MultipleToadDeformableAtlas.java` with the existing first-thinning Java
classes on the classpath. `javac` failed with 153 unresolved dependency errors,
including the absent `gov.nih.mipav.view` package, `ParametricTransform`,
`DemonsToadsAtlasWarping`, and the Medic `ImageFunctions`/utility classes.
Those dependencies are outside the committed minimal Java reference harness;
compiling the full MIPAV/JIST application is required to run the actual class.

Consequently, current Java evidence validates P18cyo's membership/objective on
this exact fixture, selected rigid coefficient/line-search behavior on its
separate small harness phantom, and wrapper control on separate Java control
cases. There is **no combined Java numeric registration digest** for this 25³
P18cyo state. The Python P18cyp result is therefore a composed implementation
using those validated kernels, not a demonstrated Java end-to-end match.

Run the focused Python composition check with:

```sh
.venv/bin/pytest -q tests/test_toads_initial_registration_control.py
```
