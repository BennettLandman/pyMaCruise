# TOADS multi-label homeomorphic-change reference

`java_reference.json` contains 130 synthetic 3×3×3 label neighborhoods and
the selected 18/6 accept/reject result. Cases include the historical
same-label unconditional acceptance and 128 deterministic random patterns
with labels 0–3. Neighborhoods are embedded at the center of a 5×5×5 short
volume so every source-style neighbor access is in bounds.

The Java harness is a focused excerpt of
`MultipleToadSegmentation.isHomeomorphicSegmentation` (source lines
1554–1673, SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`), compiled
against the actual selected `CriticalPointLUT` source and the bundled
`critical186LUT.raw.gz`. It fixes `checkComposed=false`, matching the selected
18/6 constructor configuration. The harness does not instantiate the full
segmentation plugin, `ObjectProcessing`, membership evolution, update order,
or volume-boundary behavior. It is evidence for the local label-change
predicate only, not complete TOADS segmentation equivalence.

The labels are synthetic and the fixture contains no imaging data. The
candidate is compared to the current center label; relation checks consider
every distinct neighboring label except the old and proposed labels, then
each unordered pair of such labels.

Regenerate with a JDK installed:

```sh
scripts/generate_java_toads_homeomorphic_change_reference.sh
```
