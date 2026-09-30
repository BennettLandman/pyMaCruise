# TOADS similarity/classification and LesionTOADS classification maps — 2026-09-30

## Full-suite evidence

- Command: `.venv/bin/python -m pytest -q`
- Tested implementation commit: `e0db13bc4ff1661e47bde8411e958da2b3398a49`
- Result: **693 passed, 2 expected RAS reorientation warnings, 174.50 seconds**.
- The two warnings are from the existing topology RAS test, which deliberately supplies LPS-oriented topology and paint-mask inputs and verifies conversion to RAS.
- Focused command: `.venv/bin/python -m pytest -q tests/test_toads_pair_similarity.py tests/test_toads_generate_classification.py tests/test_lesion_membership_classification.py tests/test_lesion_classification_update.py` — **12 passed in 0.69 seconds**.

The harnesses below are selected-method Java excerpts, not integrated workflow references. Synthetic fixtures are committed NIfTI files with RAS axis codes.

## P18ct — class-pair similarity dispatch

The selected VandyCRUISE `MultipleToadSegmentation.computeExactMemberships`
prelude (source SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`)
sets similarity distance to zero for LESION/LESION pairs. For mixed
LESION/non-lesion pairs it calculates `lesionDist`, but the subsequent
similarity expression still consumes `centroidDist`. The port preserves this
source behavior and float32 operation order; it does not silently substitute
the apparently intended lesion distance. Java and Python match coefficient
bits and the zero-deviation IEEE behavior. Focused module: 3 passed. This
validates the similarity prelude, not the complete membership update.

## P18cu — TOADS classification output

The selected `MultipleToadSegmentation.generateClassification` helper emits
template labels. Exact mode scans classes in ascending order from a zero
maximum and uses strict greater-than; approximate mode uses the supplied
best-index map. Pinned Java and Python match signed byte outputs for ties,
zero/negative/NaN memberships, positive infinity, and template labels distinct
from class indices. The separate wrapper uncrop is covered by P18cn and was
not duplicated here. Focused module: 2 passed. This is output-map behavior,
not segmentation equivalence.

## P39zzzzg — membership-derived lesion classification output

The selected VandyCRUISE `updateLesionClassifcationformMemberships` method
(source SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`)
selects the last case-insensitive `Cerebellum-GM` atlas index, defaulting to
index 0, then marks voxels only when lesions are enabled, current class
matches the lesion template label, and lesion membership plus that GM
membership strictly exceeds ordinary membership. It returns the marked
count. Map values and count match a pinned excerpt for float32 and alias edge
cases. Focused module: 5 passed. No selected-plugin caller was found, so this
is an available output helper, not active workflow behavior.

## P39zzzzh — active membership-to-classification update

The selected VandyCRUISE `updateClassificationFromMemberships` method
(source SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`; active
wrapper SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`)
updates ordinary byte classification from the greatest ordinary membership,
then gives the lesion label a strict-greater-than comparison. The selected
wrapper enables this update by default. Java and Python match the byte map and
changed count, including first-class ties, lesion ties/wins, template labels,
and existing classifications. Focused module: 2 passed. This validates the
helper, not a full LesionTOADS run.
