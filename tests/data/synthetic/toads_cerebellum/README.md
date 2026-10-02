# TOADS brainstem/cerebellum removal phantom

Small 9×8×7 RAS NIfTI fixture for the selected MaCRUISE second-pass tissue-map cleanup. `maseg.nii.gz` contains each source-excluded label (35, 38–41, 71–73) once, one retained label (4), and background. The GM, WM-fill, and WM-mask maps contain distinct constant float32 values. Voxel sizes are 0.8×1.1×2.0 mm. The expected cleanup zeros exactly the eight selected-label voxels in each map and preserves all other values.
