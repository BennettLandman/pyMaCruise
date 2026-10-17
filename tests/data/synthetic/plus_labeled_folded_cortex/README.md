# Labeled folded-cortex plus input case

This 65³ RAS synthetic case reuses the dedicated folded-cortex geometry and
provides same-grid T1, FLAIR, user brain mask, and explicit MaSeg-like labels:
0 background, 5 sulcal CSF, 44 white matter, and 100 cortical gray matter.
T1/FLAIR use deterministic tissue-specific values with a positive FLAIR
baseline to remain valid after cubic registration interpolation; these values
solely exercise the selected plus pipeline. They are not clinical intensity models. `provenance.json`
records image hashes, datatypes, geometry, and the assumptions.

Regenerate the four committed NIfTIs with
`./.venv/bin/python tests/data/synthetic/plus_labeled_folded_cortex/generate_fixture.py`.

After running the complete PLUS workflow on this case, audit the labeled
ribbon and surfaces with `scripts/audit_plus_ribbon_quality.py`. Add
`--require-quality` when the calling workflow should fail if any explicit
label-support, nesting, mesh-topology, or thickness diagnostic fails. The
JSON report is written either way. These checks measure synthetic ribbon
quality; they do not establish Java/MATLAB source equivalence.
