# P39ar post-registration NaN branch reference

The Java control harness pins selected VandyCRUISE flow directly after
`registerShapes()`: energy convergence comparison, iteration ceiling, optional
`propagateGrowing`, then terminal classification or relation refresh. It
covers the source-faithful P39ak → P39am NaN objective/ratio and a finite
below-threshold control. The output harness compares the final label/subtype
arrays and continuing relation map/table against the selected Java methods.
The integration test rebuilds P39ak/P39am from the folded synthetic RAS
fixture, uses the P39aq registration fixture/reference, and applies the
existing Python growth, classification, and relation kernels. For stop cases
without growth, the final segmentation is compared with the pinned P39ak Java
result. For the growth-on branch, a dynamic Java harness starts from the same
folded 17³, 13-class post-thinning state and compares segmentation, skeleton,
availability, ordering, class volumes, all 1,860 queue events, and all 991
neighbor insertions byte-for-byte. This pins selected `initGrowingLabels`,
`bestMembership`, `maxMembershipRatio`, `isHomeomorphicSegmentation`, and
`propagateGrowing` behavior on this state. The LesionTOADS branch calls
`computeRelations` only; it does not refresh base-TOADS boundary-distance or
closest-label maps. No actual imaging data is used and this is not full
workflow equivalence.

Regenerate the checked-in Java decision reference with:

```sh
sh scripts/generate_lesion_nan_post_registration_reference.sh
```

The Java predicates mirror the selected wrapper's order and original signed
`nIterMax` behavior. `Numerics.abs(float)` is expressed exactly, including its
negative NaN sign; Java still evaluates `NaN < threshold` as false.
