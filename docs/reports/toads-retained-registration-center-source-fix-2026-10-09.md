# TOADS retained registration center: source fix and held-fixed validation

## Finding

The selected `AlgorithmMultipleToads` sequence builds its
`ParametricTransform` during `initShapeRegistration`, using the first cropped
image dimensions `(166, 206, 160)`. After registration, the wrapper uncrops and
recrops the image to `(168, 207, 161)`, calls `MultipleToadDeformableAtlas.setImageInfo`,
and then computes the transformed template. That setter updates the atlas's
image metadata and centers but does not reconstruct `transformModel`; its
`imageToTemplate` calls therefore retain the original cropped-grid center.

The Python adapter had been using the refreshed crop center instead. The
initial Java/Python constructor maps differed at 183,510 of 5,598,936 voxels.
An exploratory one-voxel x shift raised class-1 Dice from 0.8753 to 0.9651,
which led to checking the source object's retained dimensions. The corrected
adapter now passes the initial crop dimensions into topology resampling after
the refresh. The initial-grid path continues to use its own crop dimensions.

The selected-source constructor map now matches exactly on the same held-fixed
input: 0 differences of 5,598,936 voxels; all 11 class-index Dice scores are
1.0. The Java/Python registration parameters and the recropped transform were
also checked at the corresponding boundaries; the latter is
`[-4.6936364, -1.3859992, -8.061735]` in both implementations.

## Held-fixed base TOADS comparisons

The selected Java V3 harness was run on the exact T1 volumes supplied to each
Python pass. The Python runs used the workflow's source settings (registration
levels/iterations `4/50`, post-refresh registration `5`, 18/6 connectivity,
and natural convergence). Both output volumes are in RAS.

| Pass | Natural stop | Java/Python output differences | Minimum per-label Dice | Python peak RSS |
|---|---:|---:|---:|---:|
| First TOADS, original N4 T1 and supplied brain mask | 11 / 11 | 0 / 6,092,928 | 1.0 | 8.495 GB |
| Second TOADS, held-fixed Dura-stripped T1 and matching mask | 24 / 24 | 0 / 6,092,928 | 1.0 | 7.380 GB |

The final energy ratios also match the Java trace: first pass
`9.77438758e-05` versus Java `9.7743876e-05`; second pass
`1.47487717e-05` versus Java `1.4748772e-05`. The first pass took 1,794.7 s
and the second 3,838.2 s. Neither run exceeded the 16 GB RSS watchdog.

| Pass | Java raw SHA-256 | Python RAS NIfTI SHA-256 |
|---|---|---|
| First | `83357ca2f01e01276d796c07c8bd984a2983e0c9f7df15390336169aecfbbba6` | `cc1da0ce5849050a1b3030da49b2acdd19e373fdd52e0229df68e1d6b391de95` |
| Second | `a220505ae8e073e1623cd4a1483b408415eccb83aeeec1592cd4cd067439bb7e` | `7f22478055b5fe3fe66d3e0856f8cf1dcebd5849ad7f8b4878a655302430d4be` |

The comparator commands were:

```sh
./.venv/bin/python scripts/compare_toads_java_classification.py \
  build/current-toads-selected-wrapper-2026-10-08/first-classification-java.raw \
  build/current-toads-selected-wrapper-2026-10-08/first-python-stale-center-correct-mask/classification.nii.gz
./.venv/bin/python scripts/compare_toads_java_classification.py \
  build/current-toads-selected-wrapper-2026-10-08/second-classification-java.raw \
  build/current-toads-selected-wrapper-2026-10-08/second-python-stale-center/classification.nii.gz
```

The focused reference command passed **5 tests**:
`./.venv/bin/pytest -q tests/test_toads_crop.py tests/test_toads_post_registration_template_refresh.py`.
The compact Java refresh harness was updated to retain the model's original
center across recropping. This is selected-source equivalence on the committed
brain-resolution synthetic workflow inputs. It does not establish complete
MaCRUISE composite equivalence, nor real-MRI validation. The corrected
seven-stage synthetic base workflow subsequently completed and resumed; see
[the integration report](base-workflow-retained-center-2026-10-09.md). The
full-suite result on the same implementation commit is recorded in
[the pytest report](full-pytest-2026-10-09-retained-center.md).

## Next integration milestone

The corrected restartable seven-stage synthetic base workflow and its
same-settings resume have now passed. Continue with the folded synthetic plus
workflow and the remaining selected-source and composite checks.
