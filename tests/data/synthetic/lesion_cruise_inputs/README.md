# LesionTOADS CRUISE-input adapter fixture

This fixture exercises the selected `LesionToadSegmentation.generateCruiseInputs()` Java method excerpt. The source file SHA-256 recorded by the harness is `77016c000b8bbbfe55bfbe8d665ddd8be914c8ecbae8563ccce59aaa834bcfe0`.

The synthetic grid is asymmetric (`5 × 4 × 3`) with RAS axis codes and voxel sizes `(0.8, 1.2, 2.5)` mm. It includes duplicate aliases for both CSF and GM to validate the Java last-matching-class overwrite behavior, plus cerebellar-prefix classes, Brainstem, Background, ordinary classes, lesion and black-hole membership maps. The label image cycles through all 13 atlas names. Membership values are distinct float32 fractions so channel routing and sequential accumulation can be checked by their exact IEEE-754 bits.

Files:

- `memberships.nii.gz`: `(x, y, z, class)` float32 memberships.
- `lesion_membership.nii.gz`, `black_hole_membership.nii.gz`: 3D float32 maps added to WM before class memberships.
- `segmentation.nii.gz`: 3D int16 atlas indices used to form the hard-WM channel.
- `java_reference.json`: exact Java float32 output bits in `[channel][x][y][z]` order.
- `generate_fixture.py`: deterministic fixture generator. Its output files are hash-pinned by the test.

The harness is a pinned method excerpt, not a build or execution of the full Java application. The pure Python adapter assumes callers already loaded all arrays in RAS voxel order; the fixture's NIfTI affines are checked independently by the test.
