# Human TGDM island lineage, attempt 02

## Result

An isolated replay of nested TGDM on the saved attempt-02 human crop localizes the small negative components to the **initial membership level set**, before surface evolution. Its `<= 0` mask is bitwise identical to the saved `wmfill_topology_corrected > 126.9` mask: **408,915 voxels, XOR 0**. Both masks contain 131 six-connected components and 82 singleton components. The surface extractor is not the source of these components.

The six saved/replayed evolution fields retain 131 six-connected components, while singleton counts change through the steps:

| TGDM field | 6-connected components | Singleton components | Negative voxels |
|---|---:|---:|---:|
| Initial membership level set | 131 | 82 | 408,915 |
| Smoothed initial level set | 131 | 130 | 418,475 |
| Initial inner level set | 131 | 126 | 394,982 |
| Central level set | 131 | 119 | 625,251 |
| Central shrink level set | 131 | 126 | 562,359 |
| Delivered inner level set | 131 | 127 | 398,034 |

The delivered inner field's 127 singleton coordinates and the central field's 119 singleton coordinates exactly match the saved attempt-02 outputs. The diagnostic also records 18- and 26-connectivity counts, crop/input hashes, RAS headers, and all six NIfTI output hashes in the ignored build artifact `build/diagnostics/human-tgdm-lineage-attempt-02-2026-10-11/lineage.json`.

## Saved-output comparison

All three replayed final fields (`inner`, `central`, `outer`) match the saved full-grid results bitwise on the region written by the selected JIST uncrop, indices `[:140, :172, :130]`: **3,130,400 voxels, zero differences**. The remaining 65,083 voxels are the union of the final x/y/z crop planes. The source uncrop uses exclusive high bounds and leaves those planes unwritten at the configured background value 20; comparing them as TGDM-produced voxels would create a false discrepancy.

This is a replay of the Python TGDM stage against the saved Python attempt, not an independent Java comparison or a new full workflow run. It identifies the earliest saved input that carries the topology, but does not establish whether the upstream WM map's small components are anatomically valid. Do not remove them without further source-based evidence.

## Reproduction and resource use

From `pymacruise/`:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_human_tgdm_lineage_watchdog.py \
  --core-attempt-dir build/invivo-base-source-v3-human-rasgrid-2026-10-11/stages/cruise-core/attempt-02 \
  --output-dir build/diagnostics/human-tgdm-lineage-attempt-02-2026-10-11 \
  --limit-gb 8 --timeout-seconds 600
```

The replay took **327.6 s**, peaked at **2.970 GiB** sampled RSS under the 8 GiB ceiling, and was neither interrupted nor terminated for memory. The focused synthetic diagnostic tests passed **2/2**. No full pytest suite was run; the selected TOADS source-line estimate is unchanged.
