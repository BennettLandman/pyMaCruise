# TOADS post-refresh exact registration (P18cyaj)

**Implementation commit:** `ee350877e61f97adb050d185bfb34f7bfd067ade`  
**Validation date:** 2026-10-02

The atlas-assets startup composition continues beyond the post-registration
crop and topology refresh. It reinitializes exact memberships on the refreshed
image grid, derives source-style transformed shape bounds, excludes the
descriptor background class from registration, and performs the selected
single-scale rigid `registerShapes()` update. Subject spacing is passed in
signed RAS form; template spacing remains the descriptor spacing. The wrapper
profile defaults to five registration iterations, while the compact synthetic
fixture explicitly overrides that value to keep validation fast.

Focused command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_alignment_startup.py tests/test_toads_post_registration_template_refresh.py tests/test_toads_registration_sampling.py tests/test_toads_registration.py
```

Result: **44 passed, 0 failed, 0 skipped in 6.37 seconds** on the implementation
commit above. The tests verify the synthetic atlas-to-registration call chain,
exact membership data handoff, background-class exclusion, final transform,
and refreshed shape bounds. This is a selected wrapper-composition check; it is
not an end-to-end comparison with Java, a TOADS segmentation iteration, or a
human-MRI run.

The full suite was also run on this exact commit: **1,049 passed, 0 failed,
0 skipped, 11 expected RAS warnings in 1,682.82 seconds**. See the [full-suite
report](full-pytest-2026-10-02-toads-atlas-registration.md).
