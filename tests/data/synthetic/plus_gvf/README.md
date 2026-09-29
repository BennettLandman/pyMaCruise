# MaCRUISEplus GVF parameter-routing phantom

The 13³ scalar image is a synthetic RAS-oriented smooth field with two
ellipsoidal intensity components and anisotropic voxel spacing. It is used to
validate that the plus stage forwards `GVF_Lambda=0.2` and 200 iterations to
the existing Python GVF adapter, while leaving vectors unnormalized. It is
not an image-based comparison with the VandyCRUISE multigrid solver.

Regenerate `scalar.nii.gz` with `make_phantom.py`.
