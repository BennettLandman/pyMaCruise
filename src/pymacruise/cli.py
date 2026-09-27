"""Small command line entry point for implemented PyMaCRUISE stages."""

from __future__ import annotations

import argparse
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import nibabel as nib

from pymacruise.cortical_separation import cortical_separation
from pymacruise.raw_preparation import prepare_raw
from pymacruise.volume_ops import calculate_nifti_files


def _package_version() -> str:
    try:
        return version("pymacruise")
    except PackageNotFoundError:
        return "0.1.0+source"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pymacruise")
    parser.add_argument("--version", action="version", version=f"%(prog)s {_package_version()}")
    commands = parser.add_subparsers(dest="command", required=True)

    prepare = commands.add_parser("prepare-raw", help="validate and prepare T1 plus MaSeg inputs")
    prepare.add_argument("--t1", required=True, type=Path)
    prepare.add_argument("--maseg", required=True, type=Path)
    prepare.add_argument("--output-dir", required=True, type=Path)
    prepare.add_argument("--overwrite", action="store_true")

    masks = commands.add_parser(
        "cortical-separation", help="derive the four MaSeg-based cortical masks"
    )
    masks.add_argument("--maseg", required=True, type=Path)
    masks.add_argument("--output-dir", required=True, type=Path)

    calculate = commands.add_parser("calculate", help="apply a generic aligned volume operation")
    calculate.add_argument("--input-a", required=True, type=Path)
    calculate.add_argument("--input-b", required=True, type=Path)
    calculate.add_argument("--operation", required=True)
    calculate.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "prepare-raw":
        result = prepare_raw(args.t1, args.maseg, args.output_dir, overwrite=args.overwrite)
        print(f"Prepared T1: {result.t1_path}")
        print(f"Prepared MaSeg: {result.maseg_path}")
    elif args.command == "cortical-separation":
        image = nib.load(str(args.maseg))
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise ValueError(f"MaSeg input is not a NIfTI image: {args.maseg}")
        outputs = cortical_separation(image, args.output_dir)
        for name in outputs:
            print(f"Wrote {args.output_dir / f'{name}.nii.gz'}")
    elif args.command == "calculate":
        calculate_nifti_files(
            str(args.input_a), str(args.input_b), str(args.output), args.operation
        )
        print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
