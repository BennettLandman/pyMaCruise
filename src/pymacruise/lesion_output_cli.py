"""Command-line materialization of already-computed LesionTOADS products.

Run with ``python -m pymacruise.lesion_output_cli --help``. This command does
not execute segmentation; it writes the selected terminal arrays produced by
the Python LesionTOADS APIs to source-named, typed RAS NIfTI files.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import nibabel as nib
import numpy as np

from .lesion_output_io import save_lesion_toads_outputs
from .lesion_post_loop_outputs import (
    LesionToadsOutputProduct,
    LesionToadsPostLoopOutputs,
)


def _load_array(path: Path) -> np.ndarray:
    image = nib.load(str(path))
    if len(image.shape) not in (3, 4):
        raise ValueError(f"expected a 3-D or 4-D NIfTI image: {path}")
    return np.asarray(image.dataobj)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pymacruise.lesion_output_cli",
        description=(
            "Write already-computed LesionTOADS terminal arrays as typed RAS "
            "NIfTI outputs. This is an output adapter, not a segmentation run."
        ),
    )
    parser.add_argument("--reference", required=True, type=Path, help="original T1 reference .nii.gz")
    parser.add_argument("--classification", required=True, type=Path, help="terminal class labels .nii.gz")
    parser.add_argument("--lesion-subtype", required=True, type=Path, help="lesion subtype labels .nii.gz")
    parser.add_argument(
        "--mode",
        required=True,
        choices=("classification-only", "class+member", "cruise inputs", "dura removal inputs"),
        help="selected terminal output branch",
    )
    parser.add_argument(
        "--membership-classification",
        type=Path,
        help="optional membership-derived class labels .nii.gz",
    )
    parser.add_argument(
        "--mode-data",
        type=Path,
        help=(
            "mode product .nii.gz: [x,y,z,channels] for class+member, or "
            "[4,x,y,z] in CSF/GM/WM-fill/WM-mask order for CRUISE and Dura"
        ),
    )
    parser.add_argument(
        "--inhomogeneity-field",
        type=Path,
        help="optional already-restored 3-D or 4-D field .nii.gz",
    )
    parser.add_argument("--input-name", required=True, help="input T1 filename used for output stem")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def _make_outputs(args: argparse.Namespace) -> tuple[LesionToadsPostLoopOutputs, np.ndarray | None]:
    classification = _load_array(args.classification)
    subtype = _load_array(args.lesion_subtype)
    if classification.ndim != 3 or subtype.ndim != 3:
        raise ValueError("classification and lesion-subtype inputs must be 3-D")
    membership_classification = (
        _load_array(args.membership_classification)
        if args.membership_classification is not None
        else None
    )
    if membership_classification is not None and membership_classification.ndim != 3:
        raise ValueError("membership-classification input must be 3-D")

    products = [LesionToadsOutputProduct("classification", classification)]
    if membership_classification is not None:
        products.append(
            LesionToadsOutputProduct(
                "classification_from_membership", membership_classification
            )
        )
    products.append(LesionToadsOutputProduct("lesion_subtype", subtype))

    mode_product = None
    if args.mode == "classification-only":
        if args.mode_data is not None:
            raise ValueError("--mode-data is not used by classification-only mode")
    else:
        if args.mode_data is None:
            raise ValueError(f"--mode-data is required for {args.mode!r} mode")
        mode_product = _load_array(args.mode_data)
        product_name = {
            "class+member": "memberships_plus_lesion",
            "cruise inputs": "cruise_inputs",
            "dura removal inputs": "dura_removal_inputs",
        }[args.mode]
        products.append(LesionToadsOutputProduct(product_name, mode_product))

    outputs = LesionToadsPostLoopOutputs(
        classification=np.asarray(classification, dtype=np.int8),
        membership_classification=(
            None
            if membership_classification is None
            else np.asarray(membership_classification, dtype=np.int8)
        ),
        lesion_classification=np.asarray(subtype, dtype=np.int16),
        result_mode=args.mode,
        products=tuple(products),
        source_order=tuple(product.name for product in products),
    )
    field = _load_array(args.inhomogeneity_field) if args.inhomogeneity_field else None
    return outputs, field


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        reference = nib.load(str(args.reference))
        outputs, field = _make_outputs(args)
        paths = save_lesion_toads_outputs(
            outputs,
            reference,
            args.output_dir,
            input_name=args.input_name,
            inhomogeneity_field=field,
        )
    except (OSError, ValueError, TypeError) as exc:
        parser.error(str(exc))
    for path in paths.values():
        print(path)
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised by subprocess test
    sys.exit(main())
