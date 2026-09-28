# Folded cortex synthetic fixture

This committed 41×41×41 RAS fixture is analytically generated and contains no
subject data. It uses a radially defined outer boundary with six azimuthal
lobes and four polar undulations, plus a 4.5 mm nested cortical shell. The
shape remains star-shaped and connected so tests can detect unintended
components or topology changes while exercising curved, narrow-band fronts.

`make_phantom.py` deterministically regenerates the NIfTI inputs and `spec.json`;
`SHA256SUMS` lists their payload checksums. Image files are compressed NIfTI.
The ACE and TGDM tests use this same phantom at separate algorithm stages.
