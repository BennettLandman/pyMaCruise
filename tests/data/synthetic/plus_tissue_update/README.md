# MaCRUISEplus atlas/TOADS tissue update phantom

This synthetic 7³ RAS fixture uses six aligned inputs with anisotropic spacing.
It includes confident eroded atlas WM, a low-valued TOADS GM fill candidate,
a GM value exactly at the strict 0.255 threshold, a fallback voxel, and a
voxel outside the atlas brain mask. The adapter routes these products through
the existing source-checked `CombineMAandTOADS` implementation with the
selected 3D erosion and `GMconfident=0.001` settings.

Regenerate all compressed NIfTI inputs with `make_phantom.py`.
