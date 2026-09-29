# Full pytest evidence — 2026-09-29 LesionTOADS centroid similarities

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest -q` from the repository root |
| Result | **427 passed in 63.35 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `91fc287a4b1dffc1f4acc701e0c44beca2816d61` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39t implements the opening similarity-coefficient setup of the selected
`LesionToadSegmentation.computeMemberships` method. Four committed Java
formula cases compare the class-to-class matrix and lesion-to-class vector
bit-for-bit, including channel-weighted terms, T1-specific Sulcal-CSF terms,
the no-T1 path, lesion-type exclusions, and canonical IEEE NaN results for
zero standard deviations. Black-hole similarity remains at its
constructor-initialized value because the update in the selected Java method
is commented out. The focused P39t test reports **1 passed**.

The harness is a direct excerpt of the selected arithmetic, not a complete
`LesionToadSegmentation` run. It does not calculate voxel memberships,
energies, atlas priors, or segmentation. Fixture and regeneration details are
in [`tests/data/synthetic/lesion_toads_similarity/README.md`](../../tests/data/synthetic/lesion_toads_similarity/README.md).
