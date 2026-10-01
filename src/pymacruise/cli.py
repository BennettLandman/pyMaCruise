"""Small command line entry point for implemented PyMaCRUISE stages."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path

import nibabel as nib

from pymacruise.ace import ace_enhance
from pymacruise.cortical_separation import cortical_separation
from pymacruise.gvf import generalized_gvf_image
from pymacruise.image_io import save_nifti_ras
from pymacruise.raw_preparation import prepare_raw
from pymacruise.surfaces import (
    extract_levelset_surface,
    read_legacy_vtk,
    write_legacy_vtk,
)
from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_images
from pymacruise.thickness import closest_point_thickness, levelset_subtraction_thickness
from pymacruise.toads_dura import run_dura_removal
from pymacruise.volume_ops import calculate_volume


def _package_version() -> str:
    try:
        return version("pymacruise")
    except PackageNotFoundError:
        return "0.1.0+source"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pymacruise")
    parser.add_argument("--version", action="version", version=f"%(prog)s {_package_version()}")
    commands = parser.add_subparsers(dest="command", required=True)

    prepare = commands.add_parser(
        "prepare-raw", help="validate and prepare T1, MaSeg, and brain mask inputs"
    )
    prepare.add_argument("--t1", required=True, type=Path)
    prepare.add_argument("--maseg", required=True, type=Path)
    prepare.add_argument("--brain-mask", required=True, type=Path)
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
    _add_provenance_option(calculate)

    ace = commands.add_parser("ace", help="run the supported MaACE/ACE enhancement stage")
    ace.add_argument("--gray", required=True, type=Path)
    ace.add_argument("--white", required=True, type=Path)
    ace.add_argument("--atlas-csf", type=Path)
    ace.add_argument("--enhanced", required=True, type=Path, help="enhanced GM .nii.gz")
    ace.add_argument("--skeleton", required=True, type=Path, help="ACE skeleton .nii.gz")
    ace.add_argument("--thinned-skeleton", required=True, type=Path, help="thinned skeleton .nii.gz")
    ace.add_argument("--diffusion-iterations", type=int, default=25)
    _add_provenance_option(ace)

    gvf = commands.add_parser("gvf", help="compute generalized GVF from a scalar image")
    gvf.add_argument("--input", required=True, type=Path)
    gvf.add_argument("--output", required=True, type=Path, help="component-last vector .nii.gz")
    gvf.add_argument("--lambda", dest="lambda_value", type=float, default=1.0)
    gvf.add_argument("--max-iterations", type=int, default=200)
    gvf.add_argument("--solver-tolerance", type=float, default=1e-6)
    _add_provenance_option(gvf)

    tgdm = commands.add_parser("tgdm", help="run supported nested CRUISE/TGDM stages")
    tgdm.add_argument("--initial", required=True, type=Path)
    tgdm.add_argument("--white", required=True, type=Path)
    tgdm.add_argument("--gray", required=True, type=Path)
    tgdm.add_argument("--ace-gray", required=True, type=Path)
    tgdm.add_argument("--brain-mask", type=Path)
    tgdm.add_argument("--ace-skeleton", type=Path)
    tgdm.add_argument("--gvf", type=Path)
    tgdm.add_argument("--output-dir", required=True, type=Path)
    _add_provenance_option(tgdm)

    dura = commands.add_parser("dura", help="run supported Dura removal stages")
    dura.add_argument("--t1", required=True, type=Path)
    dura.add_argument("--white", required=True, type=Path)
    dura.add_argument("--gray", required=True, type=Path)
    dura.add_argument("--white-mask", required=True, type=Path)
    dura.add_argument("--output-dir", required=True, type=Path)
    dura.add_argument("--max-iterations", type=int, default=5)
    _add_provenance_option(dura)

    surface = commands.add_parser("surface-extract", help="extract a RAS-mm level-set surface as VTK")
    surface.add_argument("--level-set", required=True, type=Path)
    surface.add_argument("--output", required=True, type=Path, help="legacy VTK .vtk output")
    surface.add_argument("--level", type=float, default=0.0)
    _add_provenance_option(surface)

    level_thickness = commands.add_parser(
        "thickness-levelset", help="compute voxelwise level-set subtraction thickness"
    )
    level_thickness.add_argument("--inner", required=True, type=Path)
    level_thickness.add_argument("--outer", required=True, type=Path)
    level_thickness.add_argument("--output", required=True, type=Path, help="thickness .nii.gz")
    level_thickness.add_argument("--iso-level", type=float, default=0.0)
    _add_provenance_option(level_thickness)

    surface_thickness = commands.add_parser(
        "thickness-surface", help="compute closest-point thickness on RAS-mm VTK meshes"
    )
    surface_thickness.add_argument("--inner", required=True, type=Path)
    surface_thickness.add_argument("--outer", required=True, type=Path)
    surface_thickness.add_argument("--central", type=Path)
    surface_thickness.add_argument("--inner-output", required=True, type=Path)
    surface_thickness.add_argument("--central-output", type=Path)
    _add_provenance_option(surface_thickness)
    return parser


def _add_provenance_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--provenance",
        type=Path,
        help="optional JSON sidecar recording inputs, effective parameters, and outputs",
    )


def _load_nifti_gz(path: Path) -> nib.spatialimages.SpatialImage:
    if not path.name.lower().endswith(".nii.gz"):
        raise ValueError(f"Image inputs must use .nii.gz: {path}")
    image = nib.load(str(path))
    if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
        raise ValueError(f"Input is not a NIfTI image: {path}")
    return image


def _require_nifti_gz_output(path: Path) -> None:
    if not path.name.lower().endswith(".nii.gz"):
        raise ValueError(f"Image outputs must use .nii.gz: {path}")


def _save_nifti_gz(image: nib.spatialimages.SpatialImage, path: Path) -> None:
    # This helper enforces the compressed image contract and delegates RAS
    # normalization/warnings to the shared image I/O layer.
    save_nifti_ras(image, path)


def _write_provenance(
    sidecar: Path | None,
    *,
    command: str,
    inputs: dict[str, Path | None],
    parameters: dict[str, object],
    outputs: dict[str, Path | list[Path]],
) -> None:
    if sidecar is None:
        return
    if sidecar.suffix.lower() != ".json":
        raise ValueError("Provenance sidecar must use the .json extension")

    def resolved(value: Path) -> str:
        return str(value.expanduser().resolve())

    document = {
        "schema_version": 1,
        "command": command,
        "software": {"name": "pymacruise", "version": _package_version()},
        "spatial_convention": "RAS",
        "inputs": {name: resolved(path) for name, path in inputs.items() if path is not None},
        "parameters": parameters,
        "outputs": {
            name: [resolved(item) for item in value]
            if isinstance(value, list)
            else resolved(value)
            for name, value in outputs.items()
        },
    }
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    sidecar.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote provenance: {sidecar}")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "provenance", None) is not None and args.provenance.suffix.lower() != ".json":
        raise ValueError("Provenance sidecar must use the .json extension")
    if args.command == "prepare-raw":
        result = prepare_raw(
            args.t1,
            args.maseg,
            args.brain_mask,
            args.output_dir,
            overwrite=args.overwrite,
        )
        print(f"Prepared T1: {result.t1_path}")
        print(f"Prepared MaSeg: {result.maseg_path}")
        print(f"Prepared brain mask: {result.brain_mask_path}")
    elif args.command == "cortical-separation":
        image = nib.load(str(args.maseg))
        if not isinstance(image, (nib.Nifti1Image, nib.Nifti2Image)):
            raise ValueError(f"MaSeg input is not a NIfTI image: {args.maseg}")
        outputs = cortical_separation(image, args.output_dir)
        for name in outputs:
            print(f"Wrote {args.output_dir / f'{name}.nii.gz'}")
    elif args.command == "calculate":
        _require_nifti_gz_output(args.output)
        effective_operation = args.operation.strip().upper()
        result = calculate_volume(
            _load_nifti_gz(args.input_a),
            _load_nifti_gz(args.input_b),
            effective_operation,
        )
        _save_nifti_gz(result, args.output)
        print(f"Wrote {args.output}")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"image_a": args.input_a, "image_b": args.input_b},
            parameters={"operation": effective_operation},
            outputs={"result": args.output},
        )
    elif args.command == "ace":
        for path in (args.enhanced, args.skeleton, args.thinned_skeleton):
            _require_nifti_gz_output(path)
        ace_parameters = {
            "intensity_threshold": 126.9,
            "tensor_alpha": 1.6,
            "gaussian_alpha": 1.8,
            "diffusion_iterations": args.diffusion_iterations,
            "use_updated_csf_in_ace": 2,
            "merge_method": "Merge",
        }
        result = ace_enhance(
            _load_nifti_gz(args.gray),
            _load_nifti_gz(args.white),
            multi_atlas_csf=None if args.atlas_csf is None else _load_nifti_gz(args.atlas_csf),
            **ace_parameters,
        )
        for image, path in (
            (result.enhanced_gray_matter, args.enhanced),
            (result.skeleton, args.skeleton),
            (result.thinned_skeleton, args.thinned_skeleton),
        ):
            _save_nifti_gz(image, path)
            print(f"Wrote {path}")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"gray": args.gray, "white": args.white, "atlas_csf": args.atlas_csf},
            parameters=ace_parameters,
            outputs={
                "enhanced_gray_matter": args.enhanced,
                "skeleton": args.skeleton,
                "thinned_skeleton": args.thinned_skeleton,
            },
        )
    elif args.command == "gvf":
        _require_nifti_gz_output(args.output)
        gvf_parameters = {
            "lambda_value": args.lambda_value,
            "max_iterations": args.max_iterations,
            "solver_tolerance": args.solver_tolerance,
        }
        result = generalized_gvf_image(
            _load_nifti_gz(args.input),
            **gvf_parameters,
        )
        _save_nifti_gz(result.vector_field, args.output)
        print(f"Wrote {args.output} (RAS component-last vector field)")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"scalar_image": args.input},
            parameters=gvf_parameters,
            outputs={"vector_field": args.output},
        )
    elif args.command == "tgdm":
        tgdm_parameters = TGDMBaselineParameters()
        result = tgdm_nested_images(
            _load_nifti_gz(args.initial),
            _load_nifti_gz(args.white),
            _load_nifti_gz(args.gray),
            _load_nifti_gz(args.ace_gray),
            brain_mask=None if args.brain_mask is None else _load_nifti_gz(args.brain_mask),
            ace_skeleton=None if args.ace_skeleton is None else _load_nifti_gz(args.ace_skeleton),
            gvf=None if args.gvf is None else _load_nifti_gz(args.gvf),
            parameters=tgdm_parameters,
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        tgdm_outputs = {}
        for name in ("initial_level_set", "tissue_sign", "initial_inner", "inner", "central", "outer"):
            path = args.output_dir / f"{name}.nii.gz"
            _save_nifti_gz(getattr(result, name), path)
            print(f"Wrote {path}")
            tgdm_outputs[name] = path
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={
                "initial_membership": args.initial,
                "white_matter": args.white,
                "gray_matter": args.gray,
                "ace_gray_matter": args.ace_gray,
                "brain_mask": args.brain_mask,
                "ace_skeleton": args.ace_skeleton,
                "gvf": args.gvf,
            },
            parameters=asdict(tgdm_parameters),
            outputs=tgdm_outputs,
        )
    elif args.command == "dura":
        dura_parameters = {
            "membership_level": 0.4,
            "curvature_force": 0.5,
            "alignment_force": 1.0,
            "pressure_force": 0.5,
            "alignment_alpha": 0.5,
            "max_iterations": args.max_iterations,
            "mask_translation": 2.0,
            "level_set_translation": 0.0,
            "record_passes": False,
        }
        result = run_dura_removal(
            _load_nifti_gz(args.t1),
            _load_nifti_gz(args.white),
            _load_nifti_gz(args.gray),
            _load_nifti_gz(args.white_mask),
            **dura_parameters,
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        dura_outputs = {}
        for name in (
            "stripped_t1", "corrected_white_matter", "corrected_gray_matter",
            "brain_mask_level_set", "brain_mask", "wm_topology", "gm_topology",
        ):
            path = args.output_dir / f"{name}.nii.gz"
            _save_nifti_gz(getattr(result, name), path)
            print(f"Wrote {path}")
            dura_outputs[name] = path
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={
                "t1": args.t1,
                "white_matter": args.white,
                "gray_matter": args.gray,
                "white_matter_mask": args.white_mask,
            },
            parameters=dura_parameters,
            outputs=dura_outputs,
        )
    elif args.command == "surface-extract":
        if args.output.suffix.lower() != ".vtk":
            raise ValueError("Surface outputs must use .vtk")
        mesh = extract_levelset_surface(_load_nifti_gz(args.level_set), level=args.level)
        write_legacy_vtk(mesh, args.output)
        print(f"Wrote {args.output} (RAS mm)")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"level_set": args.level_set},
            parameters={"level": args.level, "coordinate_frame": "RAS", "units": "mm"},
            outputs={"surface": args.output},
        )
    elif args.command == "thickness-levelset":
        _require_nifti_gz_output(args.output)
        result = levelset_subtraction_thickness(
            _load_nifti_gz(args.inner), _load_nifti_gz(args.outer), iso_level=args.iso_level
        )
        _save_nifti_gz(result.image, args.output)
        print(
            f"Wrote {args.output}; mean={result.mean_mm:.8g} mm, "
            f"SD={result.standard_deviation_mm:.8g} mm, "
            f"crossing_voxels={result.crossing_voxel_count}"
        )
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"inner_level_set": args.inner, "outer_level_set": args.outer},
            parameters={"iso_level": args.iso_level},
            outputs={"thickness_volume": args.output},
        )
    elif args.command == "thickness-surface":
        if args.inner_output.suffix.lower() != ".vtk":
            raise ValueError("Surface outputs must use .vtk")
        if args.central_output is not None and args.central_output.suffix.lower() != ".vtk":
            raise ValueError("Surface outputs must use .vtk")
        inner = read_legacy_vtk(args.inner)
        outer = read_legacy_vtk(args.outer)
        central = None if args.central is None else read_legacy_vtk(args.central)
        if (central is None) != (args.central_output is None):
            raise ValueError("--central and --central-output must be supplied together")
        result = closest_point_thickness(inner, outer, central=central)
        write_legacy_vtk(result.inner_surface, args.inner_output)
        if result.central_surface is not None and args.central_output is not None:
            write_legacy_vtk(result.central_surface, args.central_output)
        print(
            f"Mean surface thickness={result.mean_mm:.8g} mm, "
            f"SD={result.standard_deviation_mm:.8g} mm, "
            f"nonzero_vertices={result.nonzero_vertex_count}"
        )
        surface_outputs: dict[str, Path] = {"inner_surface_with_measurements": args.inner_output}
        if args.central_output is not None:
            surface_outputs["central_surface_with_measurements"] = args.central_output
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"inner_surface": args.inner, "outer_surface": args.outer, "central_surface": args.central},
            parameters={"central_surface_enabled": central is not None},
            outputs=surface_outputs,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
