# P40n plus boundary-to-ACE branch evidence

Implementation commit under test: `869c90c93045d426dd24e0b2f8519ab4c21b0ec6`.

Command, run from the `pymacruise` repository root:

```text
./.venv/bin/pytest -q tests/test_plus_ace_branch.py
```

Result: **2 passed, 0 failed, 0 skipped, 0 warnings** in 3.62 seconds.

The tests cover two fixed synthetic branches on aligned RAS `.nii.gz` images:

- A segmentation containing lesion label 202, which is in the selected plus
  overlooked-label list.
- The same segmentation with that region retained as ordinary label 200.

Each branch routes the segmentation through the public plus boundary adapter,
compares raw and TOADS-WM-cleaned images against the existing source-mapped
boundary operations, routes six fixed atlas/TOADS tissue images through plus
tissue update and compares GM/WM/mask with the public combiner, then routes
the cleaned boundary and updated GM/WM through plus ACE and compares all ACE
outputs with the public ACE API using the selected mode-3 merge. It also checks
RAS output geometry, confirms that the lesion/no-lesion boundaries differ,
and verifies that all loaded input image arrays remain unchanged.

The second test calls ordinary `ace_enhance` with no plus options on the
established base ACE synthetic fixture. Enhanced GM, skeleton, and thinned
skeleton match the committed Java reference arrays exactly.

The checked-in fixture README and `SHA256SUMS` describe and verify the fixed
24³ inputs. This slice starts from synthetic upstream maps; it does not run
TOADS, LesionTOADS, N4, or full MaCRUISEplus.
