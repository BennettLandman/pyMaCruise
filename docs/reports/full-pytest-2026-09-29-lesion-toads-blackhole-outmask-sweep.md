# Full pytest evidence: P39zza BLACKHOLE+OUTMASK volume composition

- Date: 2026-09-29
- Implementation commit tested: `1c89071df616d47135f2ee5c1690e36259815309`
- Command: `.venv/bin/pytest -q`
- Result: **515 passed**, 2 expected RAS orientation warnings, in 77.35 seconds.
- Focused BLACKHOLE/OUTMASK volume suite: 11 passed.

Two committed synthetic Java formula references exercise BLACKHOLE before and
after OUTMASK in the same in-place sweep. Both compare every interior class,
lesion, black-hole, and outlier float32 membership bit, including class and
auxiliary atlas priors, distinct OUTMASK neighbor terms, ventricle/GM/brainstem
distance behavior, shape gating, and shared normalization. The fixture
preserves boundaries and caller inputs. This validates the selected combined
branch, not the complete iterative LesionTOADS workflow.
