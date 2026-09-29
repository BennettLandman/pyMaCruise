# MaCRUISEplus nested surface phantom

Three 17³ synthetic RAS level sets form concentric physical spheres with
radii 3, 4, and 5 mm on an anisotropic 1.0 × 1.5 × 2.0 mm grid. They validate
the plus adapter that routes uncropped level sets through the existing
scikit-image Lewiner surface extractor. Meshes are expected in RAS+ millimetres
and to be closed; this fixture does not compare against JIST geometry.

Regenerate the compressed NIfTI level sets with `make_phantom.py`.
