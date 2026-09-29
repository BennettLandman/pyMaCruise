# Full test suite: MaCRUISEplus tissue update routing

- Date: 2026-09-29
- Tested implementation commit: `7fa402518c8490b2390f4ff06f685f7bf75d5d85`
- Command: `.venv/bin/pytest -q`
- Result: **478 passed in 58.40 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The selected plus MATLAB wrapper invokes `CombineMAandTOADS` with the six
cropped atlas/TOADS GM, WM, brain-mask, and segmentation products, `ErodeMethod`
3D, and `GMconfident=0.001`. `run_plus_tissue_update_stage` routes those six
images through the previously source-mapped image helper and returns its three
RAS products. The committed 7³ anisotropic RAS phantom covers confident
six-neighbor-eroded WM priority, the strict GM threshold (including equality),
TOADS fallback, outside-brain zeroing, output dtypes, and physical-grid
rejection. The adapter outputs match direct helper composition exactly.

Focused stage tests: 23 passed (plus tissue update, boundary, and plus option
tests). The six compressed NIfTI inputs, generator, provenance note, and
SHA-256 manifest are committed in `tests/data/synthetic/plus_tissue_update/`.
This validates plus-stage routing to the source-mapped combiner; it does not
validate the full plus orchestration or claim a direct Java harness comparison
for this adapter.
