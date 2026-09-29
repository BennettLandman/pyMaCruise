# Dura two-pass topology reference fixture

This 9×9×9 deterministic synthetic case follows the selected
`MedicAlgorithmRemoveDura09.execute()` topology call sequence:

1. Run `AlgorithmTopologyCorrection` on the WM scalar field and original WM
   paint mask (`tcWM`).
2. Add `tcWM` to the original raw GM field voxelwise as float32.
3. Run a second `AlgorithmTopologyCorrection` on that union with the same
   original WM paint mask (`tcGM`).

The harness uses the selected VandyCRUISE Java topology implementation twice
and emits the input, paint mask, raw GM, first corrected WM, float32 union, and
second corrected GM as separate sections. It follows the selected wrapper's
call order and arithmetic but does not instantiate the full JIST
`MedicAlgorithmRemoveDura09` class. The selected wrapper source SHA-256 is
`1ed180daf4a4334c995ad310828ec91cfd5593413eec17ab15404a421d7805a4`.

Python compares each of `tcWM`, the intermediate union, and `tcGM` voxelwise
against the selected Java harness reference exactly. All products are 3D
float32 NIfTI images with the same RAS grid. A separate test rejects mismatched
input grids before union arithmetic. This is one bounded wrapper-composition
case, not full Dura removal or TOADS workflow equivalence.

Regenerate this fixture together with the topology cases by running
`scripts/generate_java_topology_correction_reference.sh` from the repository
root. The generator requires a JDK. Reference products are ordinary committed
synthetic `.nii.gz` test files.
