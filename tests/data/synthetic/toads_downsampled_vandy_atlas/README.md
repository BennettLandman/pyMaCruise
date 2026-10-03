# Downsampled bundled VandyCRUISE atlas synthetic phantom

This deterministic, coarse whole-brain phantom is generated from the atlas maps already bundled in this repository. Its 36×45×34 grid preserves the source field of view and atlas labels; a synthetic T1 is synthesized from the descriptor's T1_MPRAGE class proportions with a small smooth bias. All stored images are compressed NIfTI in RAS orientation. Regenerate with `python generate_fixture.py`.

This is not clinical data and does not establish full-loop source equivalence or clinical performance. The coarse spacing is intended to make a folded whole-brain convergence/runtime experiment practical.
