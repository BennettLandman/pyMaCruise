# MaCRUISEplus boundary preparation phantom

This anisotropic 13³ RAS fixture contains two nested segmentation labels, one
configured overlooked label (120), and a partial TOADS WM mask. It validates
the plus wrapper composition: 3D multi-atlas boundary construction using the
configured overlooked-label set, followed by `toads_wm` cleanup.

Regenerate both compressed NIfTI inputs with `make_phantom.py`.
