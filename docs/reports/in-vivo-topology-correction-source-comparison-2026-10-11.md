# Human in-vivo topology-correction source comparison — 2026-10-11

## Result

The selected Java `AlgorithmTopologyCorrection` replay matches the saved
Python topology-corrected WM-fill field exactly on the 161×195×154 human crop:
**0 differing voxels of 4,834,830**, maximum and mean absolute differences both
zero, and identical value range `[0, 255.00006103515625]`.

The replay used the selected source settings: scalar image, regularization
`1e-4f`, 3D, object connectivity 6, background connectivity 18, paint mask
with strict threshold `>0.5`, levels 0–1, min/max normalization, and
object-to-background direction. At each measured superlevel threshold (0,
63.75, 127.5, 191.25, 254), Java and Python masks were identical, with 209
foreground 6-connected components and one background 18-connected component.
The paint mask already has that topology; the multi-component foreground is
present in the source input, not introduced by a Python/Java discrepancy.

The earlier read-only stage audit found the first connectivity change before
crop and topology correction: zeroing selected brainstem/cerebellum MaSeg
labels removes 159,330 voxels and changes the second-TOADS WM-fill from one to
147 6-connected components (12 at 18-connectivity). This localizes where the
connectivity transition occurs but does not establish that the selected label
removal is erroneous. Topology correction increases the count to 209 and
matches Java exactly. JIST and Python also extract the same final 209 surface
components. The remaining 97% GM/WM support gap is therefore not explained by
a human-sized topology-correction mismatch or by surface extraction; the next
source comparison should focus on TGDM's treatment of these exact inputs.

## Execution and evidence scope

The temporary Java harness, binary inputs, and full comparison are under the
ignored directory:

`build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/topology-java-human/`

The selected Java replay completed under an 8 GiB RSS watchdog with `-Xmx6g`
and **0.918 GB peak sampled direct-process RSS**. This is a held-fixed human
component source comparison. It does not establish whole-workflow equivalence,
ribbon anatomical validity, or clinical validity. No human data were committed.
