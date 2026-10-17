# AtlasEM descriptor inventory command

## Scope

Commit `c8500ea35b1ff9a82af1c87a0ecb76892ccc7694` adds `pymacruise atlasem-atlas-inventory --descriptor PATH`. It prints JSON with descriptor class order, labels/topology names, registered-class flags, key dimensions/resolutions, every referenced payload path, file presence, and `.nii.gz` extension status. The command is read-only and does not load voxel data; a zero exit status means all listed resources exist with the expected storage extension, not that their NIfTI headers or contents have been validated.

This gives a low-memory first check for the selected `ibsr8obj` bundle when it becomes available. Header and voxel validation still use the separate loader API.

## Validation

The focused CLI/parser tests passed **13 tests, 3 expected RAS warnings, in 0.90 s**. The broader AtlasEM regression passed **38 tests, 3 expected RAS warnings, in 5.92 s** on tested commit `c8500ea35b1ff9a82af1c87a0ecb76892ccc7694`. Synthetic NIfTI files verify both complete and missing-resource reports. No full project suite was run.
