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
P18cyo correction fields remain fixed. The Java references currently validate
P18cyo's membership and objective on this exact fixture, the rigid optimizer's
coefficient/line-search behavior on a separate small Java harness phantom, and
wrapper control on its own selected-source cases. They do **not** provide a
combined Java numeric registration digest for this 25³ P18cyo state, so this
composition is not yet claimed to match a Java registration run end to end.

Run the focused Python composition check with:

```sh
.venv/bin/pytest -q tests/test_toads_initial_registration_control.py
```
