# Larger folded-cortex TOADS iteration phantom

This deterministic 27×25×23 RAS-oriented synthetic case has a closed, wavy
cortical ribbon, an adjacent tissue class, a one-sided notch, nonuniform
multiclass memberships, and spatially varying shape priors. The inputs are
compressed NIfTI files; `generate_fixture.py` recreates them.

`test_toads_folded_cortex_iteration.py` sends this volume through a bounded
continuing TOADS slice: queue-initialized thinning, ordinary membership
update, the selected non-stop decision, growing, and relation/boundary/
closest-label refresh. The committed digest/count snapshot is a Python
integration regression reference, not a Java harness output and not whole-run
equivalence evidence. The component kernels retain their existing selected
Java comparisons in the tests cited by `docs/reports/toads-folded-cortex-iteration-2026-10-02.md`.

The phantom is intentionally synthetic and does not establish clinical
performance, atlas registration, field correction, or convergence of a full
MaCRUISE/TOADS run.
