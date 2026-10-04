# First-loop TOADS source diagnostic — 2026-10-03

## What the comparison covers

This is an exploratory Java-source/Python comparison on synthetic data. It is
separate from the downsampled Vandy-atlas runtime profile and does not establish
full TOADS equivalence or clinical validity. Java source ran the existing
`VandyCoreSyntheticRunExport` harness for one outer iteration. Python ran the
registered TOADS trace adapter on the same T1 values and the same 11-class raw
atlas fixture after converting the atlas to RAS NIfTI for the Python loader.

The Java harness requires MIPAV's x-fastest linear voxel order. Exporting a
NumPy NIfTI array with its default C order silently permuted the input axes and
made the earlier crop comparison invalid. Re-exporting the same T1 in Fortran
order corrected that harness error.

## Results

- Both implementations selected crop bounds `[4,33] [3,41] [1,28]` and a
  padded grid of `(50,59,48)`. This confirms crop geometry on this selected
  synthetic input after fixing the harness serialization.
- The transformed constructor topology differs at 17 of 141,600 voxels
  (99.988% exact agreement). This is a selected categorical-map comparison,
  not proof of transform equivalence.
- After the first thinning pass, 134,314 of 141,600 class-index labels agree
  exactly (94.85%).
- In the initial membership arrays, 264,798 of 1,557,600 values differ by
  more than `1e-6`; maximum absolute difference is `0.9996251`, and mean
  absolute difference is `0.0118932`.
- After the first exact membership update, 312,794 of 1,557,600 values differ
  by more than `1e-6`; maximum absolute difference is `0.9992604`, and mean
  absolute difference is `0.0120966`.
- Source initial membership energy is `6038.2886`; Python's corresponding
  objective is `5897.3169`, a 2.33% difference relative to the Java value.
  The arrays are on the same grid, but the uncertain atlas conversion means
  this is a diagnostic discrepancy, not an accepted source-equivalence score.

## Interpretation and next step

The comparison now has matching input voxel order, crop bounds, and membership
grid dimensions. The remaining discrepancy is concentrated in atlas shape
sampling/registration and the subsequent membership update, rather than the
crop extent. The 17-voxel topology disagreement alone does not explain the
membership error. The temporary NIfTI conversion of the Java atlas fixture was
diagnostic and its orientation mapping has not been independently certified;
shape-prior comparisons from that conversion are not accepted as equivalence
evidence.

Next, compare the selected Java and Python transformed shape priors at the same
registration state and validate the atlas NIfTI conversion's source axis order.
Then compare statistics and membership terms using identical constructor
labels and priors before testing any longer TOADS convergence run. No full
pytest suite was run for this diagnostic. Follow-up atlas-loader work corrected
the source-order shape check for true axis permutations and added an asymmetric
`(2,3,4)` round-trip case for both flip-only and permuted axes; the focused
`PYTHONPATH=src .venv/bin/pytest -q tests/test_atlas_assets.py` run passed **7
tests in 0.84s**. This validates the general conversion routine, not the
temporary diagnostic atlas's original axis semantics.
