# LesionTOADS standard 13-class volume fixture

`folded_inputs_ras.nii.gz` is a synthetic 5×5×5 T1_MPRAGE+FLAIR image on an
identity RAS grid. The modular cross-axis intensity pattern creates repeated
fold-like transitions; it is not imaging data. The test initializes 13
positive transformed atlas shapes per voxel, fixed direct/prior relation
flags, and fixed ventricle, GM, and brainstem distance maps. Initial class,
lesion, black-hole, outlier, and lesion-relation memberships are zero, as at
the selected first membership pass.

The wrapper's UI distance values 3/25/2 resolve through the Java adapter to
maxGM=3, maxBstem=25, and maxVent=2. The initial call has `GMflag=false`, so
the GM map does not modulate the lesion. The harness holds transformed
shapes, relation flags, and distance maps fixed; it checks the source-ordered
x/y/z update, including six-face reads after earlier voxel mutations, all
membership/energy outputs, border preservation, running objective, and final
objective. The selected FLAIR PowerTable lookup values are pinned in the
harness and their scalar lookup has separate tests.

Regenerate the Java reference with:

```sh
tests/java_reference/generate_lesion_toads_standard13_volume_reference.sh
```

This validates one fixed-state membership pass. It does not generate the
transformed atlas, relation maps, or distance maps, and does not establish
full LesionTOADS equivalence.
