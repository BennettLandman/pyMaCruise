# P39ar post-registration NaN branch reference

This small Java harness pins the selected VandyCRUISE control-flow directly
after `registerShapes()`: energy convergence comparison, iteration ceiling,
optional `propagateGrowing`, then terminal classification or relation refresh.
It covers the source-faithful P39ak → P39am NaN objective/ratio and a finite
below-threshold control. The integration test rebuilds P39ak/P39am from the
folded synthetic RAS fixture, uses the P39aq registration fixture/reference,
and applies the existing Python growth, classification, and relation kernels.
It uses no actual imaging data and does not claim full workflow equivalence.

Regenerate the checked-in Java decision reference with:

```sh
sh scripts/generate_lesion_nan_post_registration_reference.sh
```

The Java predicates mirror the selected wrapper's order and original signed
`nIterMax` behavior. `Numerics.abs(float)` is expressed exactly, including its
negative NaN sign; Java still evaluates `NaN < threshold` as false.
