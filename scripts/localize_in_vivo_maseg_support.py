#!/usr/bin/env python3
"""Localize source-faithful MaSeg support before/after refinement.

This is a read-only analysis of already-saved RAS NIfTIs. It does not invoke
the MaCRUISE pipeline or alter any image. The output directory receives CSV
and JSON tables suitable for a reproducible report.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import nibabel as nib
import numpy as np

from pymacruise.cortical_separation import GM_LABELS, WM_LABELS


DISTANCE_EDGES = (0.0, 0.2, 0.5, 1.0, 2.0, 5.0, float("inf"))
DISTANCE_NAMES = ("0_to_0.2", "0.2_to_0.5", "0.5_to_1", "1_to_2", "2_to_5", "over_5")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_ras(path: Path, ref: nib.spatialimages.SpatialImage | None = None) -> nib.Nifti1Image:
    image = nib.load(str(path))
    if len(image.shape) != 3:
        raise ValueError(f"Expected 3D NIfTI: {path} ({image.shape})")
    if nib.aff2axcodes(image.affine) != ("R", "A", "S"):
        raise ValueError(f"Expected RAS-oriented saved image: {path}")
    if ref is not None and (
        image.shape != ref.shape
        or not np.allclose(image.affine, ref.affine, rtol=0, atol=1e-4)
    ):
        raise ValueError(f"Image does not share the reference RAS grid: {path}")
    return image


def values(image: nib.spatialimages.SpatialImage, *, dtype: np.dtype) -> np.ndarray:
    array = np.asarray(image.dataobj, dtype=dtype)
    if not np.isfinite(array).all():
        raise ValueError("Input image contains non-finite values")
    return array


def class_masks(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.isin(labels, tuple(GM_LABELS)), np.isin(labels, tuple(WM_LABELS))


def support_mask(kind: str, inner: np.ndarray, outer: np.ndarray) -> np.ndarray:
    # These reproduce audit_in_vivo_base_products.py's threshold convention:
    # signed level sets <= 0 are inside; GM is between inner/outer and WM inside inner.
    if kind == "GM":
        return (inner > 0.0) & (outer <= 0.0)
    return inner <= 0.0


def summarize_classes(labels: np.ndarray, inner: np.ndarray, outer: np.ndarray) -> dict[str, dict[str, int | float]]:
    gm, wm = class_masks(labels)
    result: dict[str, dict[str, int | float]] = {}
    for kind, mask in (("GM", gm), ("WM", wm)):
        region = support_mask(kind, inner, outer)
        total = int(np.count_nonzero(mask))
        supported = int(np.count_nonzero(mask & region))
        result[kind] = {
            "voxels": total,
            "supported_voxels": supported,
            "unsupported_voxels": total - supported,
            "support_fraction": supported / total if total else 0.0,
        }
    return result


def per_label_rows(before: np.ndarray, after: np.ndarray, inner: np.ndarray, outer: np.ndarray) -> list[dict[str, int | float | str]]:
    rows: list[dict[str, int | float | str]] = []
    after_gm, after_wm = class_masks(after)
    for kind, allowed in (("GM", GM_LABELS), ("WM", WM_LABELS)):
        region = support_mask(kind, inner, outer)
        after_class = after_gm if kind == "GM" else after_wm
        for label in sorted(allowed):
            original_cohort = before == label
            after_label = after == label
            before_n = int(np.count_nonzero(original_cohort))
            after_n = int(np.count_nonzero(after_label))
            if before_n == 0 and after_n == 0:
                continue
            before_supported = int(np.count_nonzero(original_cohort & region))
            after_supported = int(np.count_nonzero(after_label & region))
            retained_same = int(np.count_nonzero(original_cohort & after_label))
            after_as_any_expected_tissue = int(np.count_nonzero(original_cohort & after_class))
            after_cohort_supported = int(np.count_nonzero(original_cohort & region & after_class))
            rows.append({
                "tissue_class": kind,
                "label": label,
                "before_voxels": before_n,
                "before_supported_voxels": before_supported,
                "before_support_fraction": before_supported / before_n if before_n else "",
                "after_voxels": after_n,
                "after_supported_voxels": after_supported,
                "after_support_fraction": after_supported / after_n if after_n else "",
                "original_cohort_retained_same_label": retained_same,
                "original_cohort_retained_expected_tissue_class": after_as_any_expected_tissue,
                "original_cohort_supported_after_refinement": after_cohort_supported,
                "original_cohort_support_fraction_after_refinement": after_cohort_supported / before_n if before_n else "",
            })
    return rows


def unsupported_bins(
    labels: np.ndarray,
    inner: np.ndarray,
    outer: np.ndarray,
    affine: np.ndarray,
    brain_mask: np.ndarray,
) -> dict[str, dict[str, object]]:
    """Count unsupported voxels by surface-side/field magnitude and RAS octant.

    The field magnitudes are deliberately called level-set values, not certified
    physical distances. TGDM reinitializes signed fields, but this report does
    not assert exact distance calibration for each saved field voxel.
    """
    gm, wm = class_masks(labels)
    gm_support = support_mask("GM", inner, outer)
    wm_support = support_mask("WM", inner, outer)
    # GM outside the ribbon is split into the inner side, outer side, and the
    # logically inconsistent overlap if the two evolved fields cross.
    gm_unsupported = gm & ~gm_support
    gm_inner = gm_unsupported & (inner <= 0) & (outer <= 0)
    gm_outer = gm_unsupported & (inner > 0) & (outer > 0)
    gm_crossed = gm_unsupported & (inner <= 0) & (outer > 0)
    wm_unsupported = wm & ~wm_support
    groups = {
        "GM_inner_side": (gm_inner, np.abs(inner)),
        "GM_outer_side": (gm_outer, outer),
        "GM_crossed_fields": (gm_crossed, np.maximum(np.abs(inner), outer)),
        "WM_outside_inner": (wm_unsupported, inner),
    }

    brain_coords = np.argwhere(brain_mask)
    if brain_coords.size:
        brain_ras = nib.affines.apply_affine(affine, brain_coords)
        y_mid = float(np.median(brain_ras[:, 1]))
        z_mid = float(np.median(brain_ras[:, 2]))
    else:
        y_mid = z_mid = 0.0

    result: dict[str, dict[str, object]] = {}
    for name, (mask, magnitude_field) in groups.items():
        n = int(np.count_nonzero(mask))
        magnitudes = magnitude_field[mask]
        hist, _ = np.histogram(magnitudes, bins=np.asarray(DISTANCE_EDGES))
        coords = np.argwhere(mask)
        if coords.size:
            ras = nib.affines.apply_affine(affine, coords)
            octants: dict[str, int] = {}
            for x_side, xcond in (("L", ras[:, 0] < 0), ("R", ras[:, 0] >= 0)):
                for ap_side, ycond in (("P", ras[:, 1] < y_mid), ("A", ras[:, 1] >= y_mid)):
                    for si_side, zcond in (("I", ras[:, 2] < z_mid), ("S", ras[:, 2] >= z_mid)):
                        octants[x_side + ap_side + si_side] = int(np.count_nonzero(xcond & ycond & zcond))
            bounds = {
                axis: [float(np.min(ras[:, i])), float(np.median(ras[:, i])), float(np.max(ras[:, i]))]
                for i, axis in enumerate(("R", "A", "S"))
            }
        else:
            octants = {}
            bounds = {}
        result[name] = {
            "unsupported_voxels": n,
            "field_magnitude_bins": {
                bin_name: int(count) for bin_name, count in zip(DISTANCE_NAMES, hist, strict=True)
            },
            "field_magnitude_note": "Bins are magnitudes of saved TGDM scalar fields in nominal mm scale; they are not an independently validated nearest-surface distance.",
            "ras_octants": octants,
            "ras_bounds_median": bounds,
        }
    result["GM_unsupported_total"] = {
        "unsupported_voxels": int(np.count_nonzero(gm_unsupported)),
        "inside_inner": int(np.count_nonzero(gm_inner)),
        "outside_outer": int(np.count_nonzero(gm_outer)),
        "crossed_fields": int(np.count_nonzero(gm_crossed)),
        "unsupported_but_neither_side": int(np.count_nonzero(gm_unsupported & ~(gm_inner | gm_outer | gm_crossed))),
    }
    return result


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def label_transition_rows(before: np.ndarray, after: np.ndarray) -> list[dict[str, int]]:
    changed = before != after
    pairs = np.stack((before[changed], after[changed]), axis=1)
    if pairs.size == 0:
        return []
    unique, counts = np.unique(pairs, axis=0, return_counts=True)
    return [
        {"before_label": int(pair[0]), "after_label": int(pair[1]), "voxels": int(count)}
        for pair, count in zip(unique, counts, strict=True)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maseg", type=Path, required=True)
    parser.add_argument("--refined", type=Path, required=True)
    parser.add_argument("--inner", type=Path, required=True)
    parser.add_argument("--outer", type=Path, required=True)
    parser.add_argument("--brain-mask", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    maseg_im = load_ras(args.maseg)
    refined_im = load_ras(args.refined, maseg_im)
    inner_im = load_ras(args.inner, maseg_im)
    outer_im = load_ras(args.outer, maseg_im)
    mask_im = load_ras(args.brain_mask, maseg_im)
    before = values(maseg_im, dtype=np.int32)
    after = values(refined_im, dtype=np.int32)
    inner = values(inner_im, dtype=np.float32)
    outer = values(outer_im, dtype=np.float32)
    mask = values(mask_im, dtype=np.float32) > 0
    if not np.equal(before, np.rint(before)).all() or not np.equal(after, np.rint(after)).all():
        raise ValueError("MaSeg maps must contain integer labels")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_label = per_label_rows(before, after, inner, outer)
    write_csv(args.output_dir / "per-label-support.csv", per_label)
    transitions = label_transition_rows(before, after)
    write_csv(args.output_dir / "label-transitions.csv", transitions)
    before_summary = summarize_classes(before, inner, outer)
    after_summary = summarize_classes(after, inner, outer)
    # Label outcomes that are generated by the historical MATLAB maxind bug
    # are made visible instead of silently being treated as missing tissue.
    before_values, before_counts = np.unique(before, return_counts=True)
    after_values, after_counts = np.unique(after, return_counts=True)
    before_count_map = {int(v): int(n) for v, n in zip(before_values, before_counts, strict=True)}
    known = GM_LABELS | WM_LABELS
    non_tissue_counts = {
        str(int(v)): {"before_voxels": before_count_map.get(int(v), 0), "after_voxels": int(n),
                      "change": int(n) - before_count_map.get(int(v), 0)}
        for v, n in zip(after_values, after_counts, strict=True)
        if int(v) not in known and int(v) != 0
    }
    introduced_label_counts = {
        str(int(v)): int(n) for v, n in zip(after_values, after_counts, strict=True)
        if int(v) not in before_count_map and int(v) != 0
    }
    increased_label_counts = {
        str(int(v)): {"before_voxels": before_count_map.get(int(v), 0), "after_voxels": int(n),
                      "increase": int(n) - before_count_map.get(int(v), 0)}
        for v, n in zip(after_values, after_counts, strict=True)
        if int(v) != 0 and int(n) > before_count_map.get(int(v), 0)
    }
    bins = {
        "before": unsupported_bins(before, inner, outer, maseg_im.affine, mask),
        "after_source_faithful_refinement": unsupported_bins(after, inner, outer, maseg_im.affine, mask),
    }
    report = {
        "schema_version": 1,
        "evidence_category": "read-only localization of saved human RAS products; not a pipeline rerun or source-equivalence comparison",
        "definition": {
            "gm_labels": sorted(GM_LABELS),
            "wm_labels": sorted(WM_LABELS),
            "gm_supported": "inner > 0 and outer <= 0",
            "wm_supported": "inner <= 0",
            "source_faithful_refinement": True,
            "field_bins": "magnitudes of saved TGDM field values; not independently validated Euclidean distance",
            "spatial_bins": "RAS left/right at x=0; anterior/posterior and superior/inferior split at respective median coordinate of supplied brain mask",
        },
        "inputs": {
            key: {"path": str(path.resolve()), "sha256": sha256(path.resolve())}
            for key, path in {
                "maseg": args.maseg,
                "refined": args.refined,
                "inner": args.inner,
                "outer": args.outer,
                "brain_mask": args.brain_mask,
            }.items()
        },
        "grid": {"shape": list(maseg_im.shape), "affine_ras": maseg_im.affine.tolist(), "zooms_mm": [float(v) for v in maseg_im.header.get_zooms()[:3]]},
        "support_before": before_summary,
        "support_after_source_faithful_refinement": after_summary,
        "per_label_csv": "per-label-support.csv",
        "label_transition_csv": "label-transitions.csv",
        "unsupported_spatial_and_field_bins": bins,
        "nonzero_refined_labels_outside_selected_gm_wm_lists": non_tissue_counts,
        "new_nonzero_label_values_not_present_before_refinement": introduced_label_counts,
        "nonzero_label_values_with_increased_counts": increased_label_counts,
        "limitations": [
            "Refinement changes labels and also changes the denominator; both current-map support and original-label cohort retention are reported in the CSV.",
            "A source-faithful historical maxind behavior can assign row-index values rather than MaSeg label IDs. These values are explicitly counted and can lower tissue-class support; this is a faithful bug, not evidence that tissue disappeared anatomically.",
            "Level-set field magnitudes are useful for localization but this analysis does not establish that every voxel's scalar value equals its exact Euclidean distance to the surface.",
            "RAS octants are coarse image-coordinate bins, not atlas-based anatomical regions.",
            "The report uses saved inputs and outputs only. It does not rerun refinement, modify any algorithm, compare against historical whole-pipeline output, or establish clinical validity.",
        ],
    }
    (args.output_dir / "support-localization.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir.resolve()), "support_before": before_summary, "support_after_source_faithful_refinement": after_summary}, indent=2))


if __name__ == "__main__":
    main()
