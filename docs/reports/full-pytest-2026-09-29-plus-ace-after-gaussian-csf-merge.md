# Full test suite: MaCRUISEplus ACE after-Gaussian CSF merge

- Date: 2026-09-29
- Tested implementation commit: `bd034aa282d600318d1cd2ec32723f3e8e08d7fa`
- Command: `.venv/bin/pytest -q`
- Result: **466 passed in 64.58 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The committed 24³ RAS synthetic phantom was run through the selected Java
MaACE wrapper with `placeToAddCSF=3` (`AfterGaussian`) and merge mode. Python
matches the after-Gaussian speed field with mean absolute error at most
`1e-6` and maximum absolute error at most `5e-6`. At positive updated-CSF
voxels, speed is exactly `0.1`; all other voxels retain the Gaussian speed
value. The resulting enhanced-GM output is within the same float tolerance,
and cleaned and thinned skeleton masks match exactly. The existing mode-2
after-diffusion Java comparison also passes unchanged.

The selected mode-3 outputs and intermediate stages are committed as
synthetic `.nii.gz` files in
`tests/data/synthetic/ace_after_gaussian/java_reference/`; all entries in
its `SHA256SUMS` verify. Focused ACE tests: 18 passed. This confirms the
selected plus wrapper's CSF merge placement on the synthetic phantom; it does
not validate full MaCRUISEplus workflow orchestration or the distinct MaACE
multi-atlas cortical-boundary prior from Huo et al.
