# Dura attached and no-dura synthetic reference

This fixture exercises the selected `MedicAlgorithmRemoveDura09.execute()` stage sequence on two generated 17 × 17 × 17 volumes. The `attached_dura` case adds a connected GM extension to the synthetic cortical shell; `no_dura` uses the same clean shell without that extension. Inputs and reference images are compressed NIfTI files on a 1 mm RAS grid. `make_inputs.py` regenerates the committed inputs.

Each case contains four inputs (`t1`, `wm`, `gm`, and `brain_mask`) and Java reference NIfTIs under `java_reference/`. The references include both topology fields, capped initial distance fields, each of the five TGDM pass level sets, the final level set, the brain mask, and the cleaned T1/WM/GM products. `trace.json` records Java per-pass update counts and hashes. `provenance.json` pins selected source and harness hashes, comparison tolerances, and scope; `fixture_hashes.json` checks the bytes of every committed NIfTI.

`tests/test_dura_attached_no_dura_java_reference.py` compares the Python adapter to these selected Java references. Topology and initial distance fields compare exactly; TGDM pass/final fields compare within 1e-6 mm; masks and cleaned products compare exactly. This is synthetic wrapper-sequence evidence only, not clinical validation or a full MaCRUISE equivalence claim.
