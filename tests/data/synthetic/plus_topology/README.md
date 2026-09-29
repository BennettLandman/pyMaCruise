# MaCRUISEplus MaWM topology-refinement phantom

This synthetic 9³ RAS case reuses the committed topology-correction scalar
volume and paint mask, and adds a binary MaSeg-derived WM cube as `MaWM`.
It validates only the plus wrapper's input routing to the already-portable
topology correction and its selected multi-atlas WM refinement. It is not a
whole MaCRUISEplus run.

Regenerate the three compressed NIfTI inputs from
`tests/data/synthetic/topology_correction/` with `make_phantom.py`.
