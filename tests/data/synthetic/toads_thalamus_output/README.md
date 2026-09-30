# Synthetic thalamus-output adapter fixture

This fixture validates only the selected VandyCRUISE
`MultipleToadSegmentation.generateThalamusOutput()` output adapter. It does not
validate thalamus evolution, MGDM, or the enclosing TOADS workflow.

The source file SHA-256 is
`63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884`; the
selected `AlgorithmMultipleToads.java` wrapper SHA-256 is
`b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`. A small,
standalone excerpt is pinned in
[`tests/java_reference/ToadsThalamusOutputReferenceHarness.java`](../../../java_reference/ToadsThalamusOutputReferenceHarness.java).
It can be re-run with a JDK by executing
`python tests/data/synthetic/toads_thalamus_output/generate_fixture.py` from
the repository root; the script compiles the harness in a temporary directory
inside this fixture and writes the Java float reference plus the NIfTI files.

The 4×3×5 membership and candidate-index images, and the two expected output
images, are compressed `.nii.gz` fixtures with 1 mm true-RAS affine geometry.
The asymmetric membership values contain both zeros and nonzero values. Atlas
names exercise exact case/space-sensitive matching, two duplicate exact
`Thalamus` entries, and near misses. Approximate candidate slots include a
voxel with both duplicate class indices and repeated matches for the later
duplicate, a repeated match for that duplicate alone, a no-match voxel, and a
voxel where only the first exact class is selected. `java_reference.json`
stores the excerpt output values, selected-source hashes, harness hash, and
fixture hashes. Each test verifies exact float32 bytes against both that
reference and the expected `.nii.gz` image. A separate check delegates
cropped-output restoration to the already validated P18cn uncrop helper.

No subject imaging data or external atlas is used by this fixture.
