# Historical N4 automatic Otsu masks

Synthetic bimodal, discrete and zero-background images with independently
computed SimpleITK masks and thresholds. The reference fixes 200 histogram
bins, a returned bin midpoint, inside value zero and outside value one to
match the historical ANTs/ITK behavior. Exact settings, reference version and
thresholds are in `provenance.json`; arrays and expected masks are in
`reference.npz`. SimpleITK is used only to generate the reference, not at runtime.
