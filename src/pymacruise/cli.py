"""Small command line entry point for implemented PyMaCRUISE stages."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import time

import nibabel as nib
import numpy as np

from pymacruise.atlas_resources import get_bundled_atlas_descriptor
from pymacruise.ace import ace_enhance
from pymacruise.cortical_separation import cortical_separation
from pymacruise.gvf import generalized_gvf_image
from pymacruise.image_io import nifti_ras_on_reference, normalize_image_to_ras, save_nifti_ras
from pymacruise.n4 import N4Config, N4DependencyError, n4_bias_correct
from pymacruise.raw_preparation import prepare_raw
from pymacruise.segmentation_refine import refine_segmentation_by_surfaces
from pymacruise.surfaces import (
    extract_levelset_surface,
    read_legacy_vtk,
    write_legacy_vtk,
)
from pymacruise.tgdm import TGDMBaselineParameters, tgdm_nested_images
from pymacruise.thickness import closest_point_thickness, levelset_subtraction_thickness
from pymacruise.topology import topology_correct
from pymacruise.toads_dura import run_dura_removal
from pymacruise.toads_output_io import save_toads_dura_wrapper_outputs
from pymacruise.toads_cerebellum import remove_brainstem_and_cerebellum
from pymacruise.toads_tissue_crop import crop_toads_tissue_maps
from pymacruise.toads_cruise_preparation import prepare_toads_cruise_inputs
from pymacruise.toads_cruise_core import run_toads_cruise_core
from pymacruise.toads_alignment_startup import register_toads_from_atlas_assets
from pymacruise.toads_atlas import (
    load_toads_atlas,
    load_toads_intensity_variance_factors,
    load_toads_optimized_factors,
)
from pymacruise.toads_registered_loop import ToadsRegisteredLoopConfig, run_registered_toads_loop
from pymacruise.volume_ops import calculate_volume
from pymacruise.workflow import load_workflow_config_file, preflight_workflow


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

    topology = commands.add_parser(
        "topology-correct", help="apply source-selected (6,18) topology correction"
    )
    topology.add_argument("--volume", required=True, type=Path, help="scalar volume .nii.gz")
    topology.add_argument("--paint-mask", required=True, type=Path, help="object paint mask .nii.gz")
    topology.add_argument("--multi-atlas-wm", type=Path, help="optional multi-atlas WM mask .nii.gz")
    topology.add_argument(
        "--output", required=True, type=Path, help="topology-corrected volume .nii.gz"
    )
    topology.add_argument(
        "--atlas-refined-output",
        type=Path,
        help="optional V3 multi-atlas-refined volume .nii.gz; requires --multi-atlas-wm",
    )
    topology.add_argument("--paint-threshold", type=float, default=0.5)
    topology.add_argument("--lowest-level", type=float, default=0.0)
    topology.add_argument("--highest-level", type=float, default=1.0)
    topology.add_argument("--regularize", type=float, default=1e-4)
    _add_provenance_option(topology)

    refine = commands.add_parser(
        "segmentation-refine",
        help="refine MaSeg labels using inner and outer level sets",
    )
    refine.add_argument("--segmentation", required=True, type=Path, help="MaSeg labels .nii.gz")
    refine.add_argument("--inner-levelset", required=True, type=Path, help="inner level set .nii.gz")
    refine.add_argument("--outer-levelset", required=True, type=Path, help="outer level set .nii.gz")
    refine.add_argument("--output", required=True, type=Path, help="refined labels .nii.gz")
    refine.add_argument("--summary", required=True, type=Path, help="JSON per-rule result summary")
    refine.add_argument("--inner-coefficient-mm", type=float, default=0.2)
    refine.add_argument("--outer-coefficient-mm", type=float, default=0.2)
    refine.add_argument("--inner-central-coefficient-mm", type=float, default=0.2)
    refine.add_argument("--outer-central-coefficient-mm", type=float, default=0.2)
    refine.add_argument(
        "--deal-with-central",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="enable central-label handling (default: enabled)",
    )
    refine.add_argument("--neighborhood-radius", type=int, default=2)
    refine.add_argument("--max-sweeps", type=int, default=50)
    refine.add_argument(
        "--source-faithful",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="preserve the historical MATLAB maxind behavior (default: enabled)",
    )
    _add_provenance_option(refine)

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

    n4 = commands.add_parser("n4", help="correct smooth intensity bias with optional ANTsPy")
    n4.add_argument("--input", required=True, type=Path)
    n4.add_argument("--output", required=True, type=Path, help="corrected image .nii.gz")
    n4.add_argument("--shrink-factor", type=int, default=4)
    n4.add_argument("--iterations", type=int, nargs="+", default=[50, 50, 50, 50])
    n4.add_argument("--convergence-tolerance", type=float, default=0.0)
    n4.add_argument("--spline-mesh", type=int, nargs=3, default=[1, 1, 1])
    n4.add_argument("--rescale-intensities", action="store_true")
    _add_provenance_option(n4)

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

    cerebellum = commands.add_parser(
        "remove-cerebellum",
        help="zero brainstem/cerebellum labels from second-pass TOADS tissue maps",
    )
    cerebellum.add_argument("--segmentation", required=True, type=Path, help="MaSeg label map .nii.gz")
    cerebellum.add_argument("--gm", required=True, type=Path, help="second-pass TOADS GM map .nii.gz")
    cerebellum.add_argument("--wm-fill", required=True, type=Path, help="second-pass TOADS WM-fill map .nii.gz")
    cerebellum.add_argument("--wm-mask", required=True, type=Path, help="second-pass TOADS WM-mask map .nii.gz")
    cerebellum.add_argument("--output-dir", required=True, type=Path)
    _add_provenance_option(cerebellum)

    tissue_crop = commands.add_parser(
        "crop-tissue-maps",
        help="crop second-pass GM, WM-fill, and WM-mask maps using GM support",
    )
    tissue_crop.add_argument("--gm", required=True, type=Path)
    tissue_crop.add_argument("--wm-fill", required=True, type=Path)
    tissue_crop.add_argument("--wm-mask", required=True, type=Path)
    tissue_crop.add_argument("--output-dir", required=True, type=Path)
    tissue_crop.add_argument("--edge-width", type=int, default=5)
    _add_provenance_option(tissue_crop)

    cruise_inputs = commands.add_parser(
        "prepare-cruise-inputs",
        help="optionally clean, crop, and topology-correct second-pass TOADS maps before ACE",
    )
    cruise_inputs.add_argument("--gm", required=True, type=Path)
    cruise_inputs.add_argument("--wm-fill", required=True, type=Path)
    cruise_inputs.add_argument("--wm-mask", required=True, type=Path)
    cruise_inputs.add_argument("--segmentation", type=Path, help="required with --remove-cerebellum")
    cruise_inputs.add_argument("--remove-cerebellum", action="store_true")
    cruise_inputs.add_argument("--edge-width", type=int, default=5)
    cruise_inputs.add_argument("--output-dir", required=True, type=Path)
    _add_provenance_option(cruise_inputs)

    cruise_core = commands.add_parser(
        "toads-cruise-core",
        help="run second-pass preparation, ACE, GVF, and nested TGDM",
    )
    cruise_core.add_argument("--gm", required=True, type=Path)
    cruise_core.add_argument("--wm-fill", required=True, type=Path)
    cruise_core.add_argument("--wm-mask", required=True, type=Path)
    cruise_core.add_argument("--segmentation", type=Path)
    cruise_core.add_argument("--remove-cerebellum", action="store_true")
    cruise_core.add_argument("--edge-width", type=int, default=5)
    cruise_core.add_argument("--ace-diffusion-iterations", type=int, default=25)
    cruise_core.add_argument("--gvf-lambda", type=float, default=1.0)
    cruise_core.add_argument("--gvf-iterations", type=int, default=200)
    cruise_core.add_argument("--gvf-solver-tolerance", type=float, default=1e-6)
    cruise_core.add_argument("--tgdm-initial-iterations", type=int, default=6)
    cruise_core.add_argument("--tgdm-inner-iterations", type=int, default=4)
    cruise_core.add_argument("--tgdm-central-iterations", type=int, default=7)
    cruise_core.add_argument("--output-dir", required=True, type=Path)
    _add_provenance_option(cruise_core)

    toads = commands.add_parser(
        "toads-loop",
        help="run the selected registered base-TOADS loop from an N4-corrected RAS T1 and brain mask",
    )
    toads.add_argument("--t1-n4", required=True, type=Path)
    toads.add_argument("--brain-mask", required=True, type=Path)
    toads.add_argument("--output", required=True, type=Path, help="native-grid classification .nii.gz")
    toads.add_argument("--dura-output-dir", type=Path, help="optional directory for the four source-scaled Dura-removal inputs")
    toads.add_argument("--atlas-descriptor", type=Path)
    toads.add_argument("--modality", default="T1_MPRAGE")
    toads.add_argument("--registration-levels", type=int, default=4)
    toads.add_argument("--registration-iterations", type=int, default=50)
    toads.add_argument(
        "--post-refresh-registration-iterations",
        type=int,
        default=1,
        help="rigid atlas updates after template refresh (selected Java default: 1)",
    )
    toads.add_argument("--registration-max-iterations", type=int, default=5)
    toads.add_argument("--padding", type=int, default=10)
    toads.add_argument("--maximum-iterations", type=int, default=50)
    toads.add_argument("--convergence-threshold", type=float, default=1e-4)
    toads.add_argument("--relation-limit", type=int, default=3)
    toads.add_argument("--maximum-distance", type=float, default=1e-4)
    toads.add_argument(
        "--category-exponents",
        type=float,
        nargs=3,
        default=None,
        help="override atlas CSF/GM/WM optimized factors (default: use atlas values)",
    )
    _add_provenance_option(toads)

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

    workflow_preflight = commands.add_parser(
        "workflow-preflight",
        help="validate a workflow config and write a dry-run manifest; no workflow stages execute",
    )
    workflow_preflight.add_argument("--config", required=True, type=Path)
    workflow_preflight.add_argument(
        "--overwrite",
        action="store_true",
        help="replace an existing run.json manifest only; do not remove other files",
    )
    base_run = commands.add_parser(
        "base-run",
        help="execute the restartable synthetic/native-grid base MaCRUISE stage chain",
    )
    base_run.add_argument("--t1", required=True, type=Path, help="input T1 .nii.gz")
    base_run.add_argument("--brain-mask", required=True, type=Path, help="same-grid brain mask .nii.gz")
    base_run.add_argument("--atlas-descriptor", required=True, type=Path)
    base_run.add_argument("--maseg", required=True, type=Path, help="same-grid MaSeg labels for cleanup and surface refinement")
    base_run.add_argument("--output-dir", required=True, type=Path)
    base_run.add_argument("--resume", action="store_true", help="resume from verified stage outputs in run.json")
    base_run.add_argument("--registration-levels", type=int, default=1)
    base_run.add_argument("--registration-iterations", type=int, default=1)
    base_run.add_argument("--post-refresh-registration-iterations", type=int, default=1)
    base_run.add_argument("--registration-max-iterations", type=int, default=1)
    base_run.add_argument("--toads-iterations", type=int, default=2)
    base_run.add_argument("--dura-iterations", type=int, default=1)
    base_run.add_argument("--crop-edge-width", type=int, default=2)
    base_run.add_argument("--ace-diffusion-iterations", type=int, default=1)
    base_run.add_argument("--gvf-iterations", type=int, default=200)
    base_run.add_argument("--tgdm-iterations", type=int, default=1)

    plus_register = commands.add_parser(
        "plus-register-flair",
        help="register plus FLAIR to the T1 grid and persist a resumable ANTsPy checkpoint",
    )
    plus_register.add_argument("--t1", required=True, type=Path, help="fixed T1 reference .nii.gz")
    plus_register.add_argument("--flair", required=True, type=Path, help="moving FLAIR .nii.gz")
    plus_register.add_argument("--output-dir", required=True, type=Path)
    plus_register.add_argument("--resume", action="store_true")
    plus_register.add_argument("--metric", choices=("mattes", "meansquares", "GC"), default="mattes")
    plus_register.add_argument("--metric-bins", type=int, default=32)
    plus_register.add_argument("--sampling-rate", type=float, default=0.2)
    plus_register.add_argument("--iterations", type=int, nargs="+", default=[1000, 500, 250, 0])
    plus_register.add_argument("--shrink-factors", type=int, nargs="+", default=[6, 4, 2, 1])
    plus_register.add_argument("--smoothing-sigmas", type=float, nargs="+", default=[3.0, 2.0, 1.0, 0.0])

    plus_prepare_lesion = commands.add_parser(
        "plus-prepare-lesion",
        help="register FLAIR, N4-correct both modalities, and apply the supplied T1-grid brain mask",
    )
    plus_prepare_lesion.add_argument("--t1", required=True, type=Path, help="input T1 .nii.gz")
    plus_prepare_lesion.add_argument("--flair", required=True, type=Path, help="moving FLAIR .nii.gz")
    plus_prepare_lesion.add_argument("--brain-mask", required=True, type=Path, help="same-grid T1 brain mask .nii.gz")
    plus_prepare_lesion.add_argument("--output-dir", required=True, type=Path)
    plus_prepare_lesion.add_argument("--resume", action="store_true")
    plus_prepare_lesion.add_argument("--metric", choices=("mattes", "meansquares", "GC"), default="mattes")
    plus_prepare_lesion.add_argument("--metric-bins", type=int, default=32)
    plus_prepare_lesion.add_argument("--sampling-rate", type=float, default=0.2)
    plus_prepare_lesion.add_argument("--registration-iterations", type=int, nargs="+", default=[1000, 500, 250, 0])
    plus_prepare_lesion.add_argument("--shrink-factors", type=int, nargs="+", default=[6, 4, 2, 1])
    plus_prepare_lesion.add_argument("--smoothing-sigmas", type=float, nargs="+", default=[3.0, 2.0, 1.0, 0.0])
    plus_prepare_lesion.add_argument("--n4-shrink-factor", type=int, default=4)
    plus_prepare_lesion.add_argument("--n4-iterations", type=int, nargs="+", default=[50, 50, 50, 50])
    plus_prepare_lesion.add_argument("--n4-convergence-tolerance", type=float, default=0.0)
    plus_prepare_lesion.add_argument("--n4-spline-mesh", type=int, nargs=3, default=[1, 1, 1])

    plus_crop_lesion = commands.add_parser(
        "plus-crop-lesion-inputs",
        help="normalize source-ordered T1/FLAIR modalities and prepare the padded LesionTOADS crop",
    )
    plus_crop_lesion.add_argument("--modalities", required=True, type=Path, help="4-D T1-then-FLAIR .nii.gz")
    plus_crop_lesion.add_argument("--topology", required=True, type=Path, help="LesionTOADS topology atlas .nii.gz")
    plus_crop_lesion.add_argument("--output-dir", required=True, type=Path)
    plus_crop_lesion.add_argument("--padding", type=int, default=10)
    plus_crop_lesion.add_argument("--resume", action="store_true")

    plus_lesion_align = commands.add_parser(
        "plus-lesion-initial-align",
        help="initialize source-ordered 13-class LesionTOADS memberships and rigid atlas alignment",
    )
    plus_lesion_align.add_argument("--cropped-modalities", required=True, type=Path)
    plus_lesion_align.add_argument("--crop-manifest", required=True, type=Path)
    plus_lesion_align.add_argument("--atlas-descriptor", required=True, type=Path)
    plus_lesion_align.add_argument("--output-dir", required=True, type=Path)
    plus_lesion_align.add_argument("--init-align-iterations", type=int, default=10)
    plus_lesion_align.add_argument("--pyramid-levels", type=int, default=2)
    plus_lesion_align.add_argument("--resume", action="store_true")

    plus_lesion_recrop = commands.add_parser(
        "plus-lesion-registered-recrop",
        help="apply the post-registration uncrop, recrop, and atlas-grid refresh checkpoint",
    )
    plus_lesion_recrop.add_argument("--normalized-modalities", required=True, type=Path)
    plus_lesion_recrop.add_argument("--crop-manifest", required=True, type=Path)
    plus_lesion_recrop.add_argument("--initial-alignment-manifest", required=True, type=Path)
    plus_lesion_recrop.add_argument("--atlas-descriptor", required=True, type=Path)
    plus_lesion_recrop.add_argument("--output-dir", required=True, type=Path)
    plus_lesion_recrop.add_argument("--resume", action="store_true")

    plus_lesion_first = commands.add_parser(
        "plus-lesion-first-iteration",
        help="compute registered LesionTOADS initial memberships and first thinning pass",
    )
    plus_lesion_first.add_argument("--registered-recrop-manifest", required=True, type=Path)
    plus_lesion_first.add_argument("--atlas-descriptor", required=True, type=Path)
    plus_lesion_first.add_argument("--output-dir", required=True, type=Path)
    plus_lesion_first.add_argument("--iterations", type=int, default=99)
    plus_lesion_first.add_argument("--first-thinning-limit", type=float, default=3.0)
    plus_lesion_first.add_argument("--first-thinning-volume-ratio", type=float, default=0.5)
    plus_lesion_first.add_argument("--finite-zero-variance", action="store_true")
    plus_lesion_first.add_argument("--resume", action="store_true")

    plus_lesion_branch = commands.add_parser(
        "plus-lesion-registered-branch",
        help="continue registered plus LesionTOADS through one post-thinning branch",
    )
    plus_lesion_branch.add_argument("--registered-recrop-manifest", required=True, type=Path)
    plus_lesion_branch.add_argument("--atlas-descriptor", required=True, type=Path)
    plus_lesion_branch.add_argument("--output-dir", required=True, type=Path)
    plus_lesion_branch.add_argument("--iterations", type=int, default=99)
    plus_lesion_branch.add_argument("--first-thinning-limit", type=float, default=3.0)
    plus_lesion_branch.add_argument("--first-thinning-volume-ratio", type=float, default=0.5)
    plus_lesion_branch.add_argument("--max-distance", type=float, default=0.0)
    plus_lesion_branch.add_argument("--finite-zero-variance", action="store_true")
    plus_lesion_branch.add_argument("--resume", action="store_true")

    plus_lesion_loop = commands.add_parser(
        "plus-lesion-registered-loop",
        help="continue and checkpoint a bounded registered plus LesionTOADS loop segment",
    )
    plus_lesion_loop.add_argument("--registered-checkpoint-manifest", required=True, type=Path)
    plus_lesion_loop.add_argument("--atlas-descriptor", required=True, type=Path)
    plus_lesion_loop.add_argument("--output-dir", required=True, type=Path)
    plus_lesion_loop.add_argument("--iterations", type=int, default=99)
    plus_lesion_loop.add_argument("--thinning-limit", type=float, default=3.0)
    plus_lesion_loop.add_argument("--max-passes", type=int, default=1)
    plus_lesion_loop.add_argument("--finite-zero-variance", action="store_true")
    plus_lesion_loop.add_argument("--resume", action="store_true")

    plus_run = commands.add_parser(
        "plus-run",
        help="run the restartable plus continuation from persisted P39 Dura channels and aligned maps",
    )
    plus_run.add_argument("--t1-reference", required=True, type=Path)
    plus_run.add_argument("--segmentation", required=True, type=Path)
    plus_run.add_argument("--dura-csf", required=True, type=Path)
    plus_run.add_argument("--dura-gm", required=True, type=Path)
    plus_run.add_argument("--dura-wm-fill", required=True, type=Path)
    plus_run.add_argument("--dura-wm-mask", required=True, type=Path)
    plus_run.add_argument("--multi-atlas-gm", required=True, type=Path)
    plus_run.add_argument("--multi-atlas-wm", required=True, type=Path)
    plus_run.add_argument("--multi-atlas-mask", required=True, type=Path)
    plus_run.add_argument("--brain-mask", required=True, type=Path)
    plus_run.add_argument("--output-dir", required=True, type=Path)
    plus_run.add_argument("--resume", action="store_true")
    plus_run.add_argument("--crop-edge-width", type=int, default=5)
    plus_run.add_argument("--gvf-lambda", type=float, default=0.2)
    plus_run.add_argument("--wm-iso-level", type=float, default=0.5)
    plus_run.add_argument(
        "--refine-topology-by-maseg",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
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
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "provenance", None) is not None and args.provenance.suffix.lower() != ".json":
        raise ValueError("Provenance sidecar must use the .json extension")
    if args.command == "workflow-preflight":
        config = load_workflow_config_file(args.config)
        manifest_path, manifest = preflight_workflow(
            config, overwrite=args.overwrite
        )
        print(f"Preflight passed: {manifest_path}")
        print("Dry run only; no N4 correction or segmentation workflow was executed.")
        print(f"Workflow profile: {manifest['workflow_profile']}")
    elif args.command == "base-run":
        from pymacruise.base_workflow import BaseRunOptions, run_base_workflow

        manifest_path = run_base_workflow(
            t1_path=args.t1,
            brain_mask_path=args.brain_mask,
            atlas_descriptor=args.atlas_descriptor,
            maseg_path=args.maseg,
            output_dir=args.output_dir,
            resume=args.resume,
            options=BaseRunOptions(
                registration_levels=args.registration_levels,
                registration_iterations=args.registration_iterations,
                post_refresh_registration_iterations=args.post_refresh_registration_iterations,
                registration_max_iterations=args.registration_max_iterations,
                toads_iterations=args.toads_iterations,
                dura_iterations=args.dura_iterations,
                crop_edge_width=args.crop_edge_width,
                ace_diffusion_iterations=args.ace_diffusion_iterations,
                gvf_iterations=args.gvf_iterations,
                tgdm_iterations=args.tgdm_iterations,
            ),
        )
        print(f"Base run complete: {manifest_path} (synthetic/native-grid route; inspect stage checks and TOADS stop metadata)")
    elif args.command == "plus-register-flair":
        from pymacruise.plus_workflow import run_plus_flair_registration_checkpoint
        from pymacruise.registration import RegistrationConfig

        result = run_plus_flair_registration_checkpoint(
            args.t1,
            args.flair,
            args.output_dir,
            config=RegistrationConfig(
                transform_type="Rigid",
                metric=args.metric,
                metric_bins=args.metric_bins,
                sampling_rate=args.sampling_rate,
                iterations=tuple(args.iterations),
                shrink_factors=tuple(args.shrink_factors),
                smoothing_sigmas=tuple(args.smoothing_sigmas),
            ),
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Registered"
        print(f"{verb} FLAIR on the T1 RAS grid: {result.warped_flair_path}")
        print(f"Registration manifest: {result.manifest_path}")
    elif args.command == "plus-prepare-lesion":
        from pymacruise.plus_lesion_preprocessing import run_plus_lesion_preprocessing_checkpoint
        from pymacruise.registration import RegistrationConfig

        result = run_plus_lesion_preprocessing_checkpoint(
            args.t1,
            args.flair,
            args.brain_mask,
            args.output_dir,
            registration_config=RegistrationConfig(
                transform_type="Rigid",
                metric=args.metric,
                metric_bins=args.metric_bins,
                sampling_rate=args.sampling_rate,
                iterations=tuple(args.registration_iterations),
                shrink_factors=tuple(args.shrink_factors),
                smoothing_sigmas=tuple(args.smoothing_sigmas),
            ),
            n4_config=N4Config(
                shrink_factor=args.n4_shrink_factor,
                iterations=tuple(args.n4_iterations),
                convergence_tolerance=args.n4_convergence_tolerance,
                spline_mesh=tuple(args.n4_spline_mesh),
            ),
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Prepared"
        print(f"{verb} source-ordered T1/FLAIR inputs: {result.modalities_path}")
        print(f"Lesion-preprocessing manifest: {result.manifest_path}")
    elif args.command == "plus-crop-lesion-inputs":
        from pymacruise.plus_lesion_crop import run_plus_lesion_crop_checkpoint

        result = run_plus_lesion_crop_checkpoint(
            args.modalities,
            args.topology,
            args.output_dir,
            padding=args.padding,
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Normalized and cropped"
        print(f"{verb} LesionTOADS modalities: {result.cropped_path}")
        print(f"Inclusive source bounds: {result.crop_bounds}")
        print(f"Lesion-crop manifest: {result.manifest_path}")
    elif args.command == "plus-lesion-initial-align":
        from pymacruise.plus_lesion_atlas_init import run_plus_lesion_initial_atlas_alignment_checkpoint

        result = run_plus_lesion_initial_atlas_alignment_checkpoint(
            args.cropped_modalities,
            args.crop_manifest,
            args.atlas_descriptor,
            args.output_dir,
            init_align_iterations=args.init_align_iterations,
            pyramid_levels=args.pyramid_levels,
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Initialized"
        print(f"{verb} LesionTOADS initial atlas alignment: {result.memberships_path}")
        print(f"Initial rigid transform: {result.transform}")
        print(f"Initial-alignment manifest: {result.manifest_path}")
    elif args.command == "plus-lesion-registered-recrop":
        from pymacruise.plus_lesion_registered_recrop import run_plus_lesion_registered_recrop_checkpoint

        result = run_plus_lesion_registered_recrop_checkpoint(
            args.normalized_modalities,
            args.crop_manifest,
            args.initial_alignment_manifest,
            args.atlas_descriptor,
            args.output_dir,
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Completed"
        print(f"{verb} LesionTOADS registered recrop: {result.modalities_path}")
        print(f"Transformed topology: {result.topology_path}")
        print(f"Registered shape maps: {result.shapes_path}")
        print(f"Registered-recrop manifest: {result.manifest_path}")
    elif args.command == "plus-lesion-first-iteration":
        from pymacruise.plus_lesion_first_iteration import run_plus_lesion_first_iteration_checkpoint

        result = run_plus_lesion_first_iteration_checkpoint(
            args.registered_recrop_manifest,
            args.atlas_descriptor,
            args.output_dir,
            iterations=args.iterations,
            first_thinning_limit=args.first_thinning_limit,
            first_thinning_volume_ratio=args.first_thinning_volume_ratio,
            source_zero_variance_nan=not args.finite_zero_variance,
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Completed"
        print(f"{verb} LesionTOADS initial memberships: {result.memberships_path}")
        print(f"First thinning segmentation: {result.segmentation_path}")
        print(f"First-iteration manifest: {result.manifest_path}")
    elif args.command == "plus-lesion-registered-branch":
        from pymacruise.plus_lesion_registered_branch import run_plus_lesion_registered_branch_checkpoint

        result = run_plus_lesion_registered_branch_checkpoint(
            args.registered_recrop_manifest,
            args.atlas_descriptor,
            args.output_dir,
            iterations=args.iterations,
            first_thinning_limit=args.first_thinning_limit,
            first_thinning_volume_ratio=args.first_thinning_volume_ratio,
            source_zero_variance_nan=not args.finite_zero_variance,
            max_distance=args.max_distance,
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Completed"
        print(f"{verb} registered LesionTOADS branch: {result.manifest_path}")
        print(f"Updated segmentation: {result.segmentation_path}")
    elif args.command == "plus-lesion-registered-loop":
        from pymacruise.plus_lesion_registered_loop import persist_plus_lesion_registered_loop_checkpoint

        result = persist_plus_lesion_registered_loop_checkpoint(
            args.registered_checkpoint_manifest,
            args.atlas_descriptor,
            args.output_dir,
            n_iter_max=args.iterations,
            thinning_limit=args.thinning_limit,
            source_zero_variance_nan=not args.finite_zero_variance,
            max_passes=args.max_passes,
            resume=args.resume,
        )
        verb = "Verified and resumed" if result.resumed else "Completed"
        print(f"{verb} registered LesionTOADS loop segment: {result.manifest_path}")
        print(f"Updated segmentation: {result.segmentation_path}")
    elif args.command == "plus-run":
        from pymacruise.plus_file_workflow import run_plus_file_workflow
        from pymacruise.plus_pipeline import PlusPipelineOptions

        result = run_plus_file_workflow(
            t1_reference_path=args.t1_reference,
            segmentation_path=args.segmentation,
            dura_csf_path=args.dura_csf,
            dura_gm_path=args.dura_gm,
            dura_wm_fill_path=args.dura_wm_fill,
            dura_wm_mask_path=args.dura_wm_mask,
            multi_atlas_gm_path=args.multi_atlas_gm,
            multi_atlas_wm_path=args.multi_atlas_wm,
            multi_atlas_mask_path=args.multi_atlas_mask,
            brain_mask_path=args.brain_mask,
            output_dir=args.output_dir,
            crop_edge_width=args.crop_edge_width,
            options=PlusPipelineOptions(
                gvf_lambda=args.gvf_lambda,
                wm_iso_level=args.wm_iso_level,
                refine_topology_by_maseg=args.refine_topology_by_maseg,
            ),
            resume=args.resume,
        )
        print(f"Plus continuation complete: {result.manifest_path}")
        print(f"Verified stages skipped: {', '.join(result.resumed_stages) or 'none'}")
        for name, path in result.outputs.items():
            print(f"{name}: {path}")
    elif args.command == "toads-loop":
        _require_nifti_gz_output(args.output)
        if args.output.exists():
            raise FileExistsError(f"Refusing to overwrite TOADS output: {args.output}")
        progress_started = time.perf_counter()

        def report_toads_progress(message: str) -> None:
            elapsed = time.perf_counter() - progress_started
            print(f"TOADS progress +{elapsed:.1f}s: {message}", flush=True)

        report_toads_progress("loading descriptor and atlas registration inputs")
        descriptor = args.atlas_descriptor or get_bundled_atlas_descriptor(
            "cruise-atlas-10obj-toads2008.txt"
        )
        atlas = load_toads_atlas(descriptor)
        t1_n4 = _load_nifti_gz(args.t1_n4)
        brain_mask = _load_nifti_gz(args.brain_mask)
        report_toads_progress("starting atlas alignment")
        alignment = register_toads_from_atlas_assets(
            t1_n4,
            brain_mask,
            atlas,
            levels=args.registration_levels,
            max_iterations=args.registration_iterations,
            refinement_iterations=args.post_refresh_registration_iterations,
            modality=args.modality,
            padding=args.padding,
            relation_limit=args.relation_limit,
            energy_tolerance=1e-6,
            minimum_lambda=1e-3,
            progress_callback=report_toads_progress,
        )
        class_indices = {"optimized_csf": 0, "optimized_gm": 1, "optimized_wm": 2}
        if args.category_exponents is None:
            category_exponents = load_toads_optimized_factors(
                atlas.descriptor, [args.modality]
            )
        else:
            category_exponents = np.asarray([args.category_exponents], dtype=np.float32)
        exponent_by_class = np.empty((1, len(atlas.structures)), dtype=np.float32)
        for index, structure in enumerate(atlas.structures):
            topology = structure.topology.lower()
            if topology == "mask":
                exponent_by_class[0, index] = np.float32(1.0)
            elif topology in class_indices:
                exponent_by_class[0, index] = category_exponents[0, class_indices[topology]]
            else:
                raise ValueError(f"Unsupported base-TOADS structure topology: {structure.topology!r}")
        config = ToadsRegisteredLoopConfig(
            category_exponents=category_exponents,
            exponent_by_channel_class=exponent_by_class,
            modality_weights=np.ones(1, dtype=np.float32),
            class_variance_weights=load_toads_intensity_variance_factors(
                atlas.descriptor, [args.modality]
            ),
            smoothing_priors=np.ones(len(atlas.structures), dtype=np.float32),
            atlas_shape_priors=np.ones(len(atlas.structures), dtype=np.float32),
            smoothing=0.2 / 6.0 / len(atlas.structures),
            # AlgorithmMultipleToads hard-codes atlasCoupling=0.5; this is
            # forwarded as priorFactor to MultipleToadSegmentation.
            prior_factor=0.5,
            prior_coefficient=2.0 / len(atlas.structures),
            maximum_iterations=args.maximum_iterations,
            convergence_threshold=args.convergence_threshold,
            relation_limit=args.relation_limit,
            do_growing=True,
            maximum_polynomial_degree=3,
            maximum_distance=args.maximum_distance,
            max_registration_iterations=args.registration_max_iterations,
        )
        report_toads_progress("starting registered TOADS segmentation loop")
        result = run_registered_toads_loop(
            alignment, config, progress_callback=report_toads_progress
        )
        report_toads_progress("registered TOADS segmentation loop complete")
        if result.native_classification is None:
            raise RuntimeError("selected TOADS loop returned without terminal classification")
        ras_t1 = normalize_image_to_ras(t1_n4, source_name="N4-corrected T1 output reference")
        output_image = nifti_ras_on_reference(
            ras_t1,
            result.native_classification,
            dtype=np.uint8,
        )
        _save_nifti_gz(output_image, args.output)
        outputs: dict[str, Path | list[Path]] = {"classification": args.output}
        if args.dura_output_dir is not None:
            if result.wrapper_output_buffers is None:
                raise RuntimeError("selected TOADS loop did not return Dura-removal wrapper outputs")
            dura_outputs = save_toads_dura_wrapper_outputs(
                result.wrapper_output_buffers,
                ras_t1,
                args.dura_output_dir,
                input_name=args.t1_n4,
            )
            outputs["dura_removal_inputs"] = list(dura_outputs.values())
            for name, path in dura_outputs.items():
                print(f"Wrote Dura-removal input {name}: {path}")
        print(
            f"Wrote {args.output}; iterations={result.iteration_count}, "
            f"changed_classification_voxels={result.changed_classification_voxels}, "
            f"terminal_energy_ratio={float(result.terminal_energy_ratio):.9g}, "
            f"energy_threshold_reached={str(result.energy_threshold_reached).lower()}, "
            f"iteration_limit_reached={str(result.iteration_limit_reached).lower()} (RAS)"
        )
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"t1_n4": args.t1_n4, "brain_mask": args.brain_mask},
            parameters={
                "atlas_descriptor": str(descriptor),
                "modality": args.modality,
                "registration_levels": args.registration_levels,
                "registration_iterations": args.registration_iterations,
                "maximum_iterations": args.maximum_iterations,
                "maximum_distance": args.maximum_distance,
                "category_exponents": category_exponents[0].tolist(),
            },
            outputs=outputs,
        )
    elif args.command == "prepare-raw":
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
    elif args.command == "topology-correct":
        _require_nifti_gz_output(args.output)
        if args.atlas_refined_output is not None:
            _require_nifti_gz_output(args.atlas_refined_output)
        if (args.multi_atlas_wm is None) != (args.atlas_refined_output is None):
            parser.error("--multi-atlas-wm and --atlas-refined-output must be supplied together")
        parameters = {
            "paint_threshold": args.paint_threshold,
            "lowest_level": args.lowest_level,
            "highest_level": args.highest_level,
            "regularize": args.regularize,
        }
        result = topology_correct(
            _load_nifti_gz(args.volume),
            _load_nifti_gz(args.paint_mask),
            multi_atlas_wm=(
                None if args.multi_atlas_wm is None else _load_nifti_gz(args.multi_atlas_wm)
            ),
            **parameters,
        )
        _save_nifti_gz(result.corrected, args.output)
        outputs: dict[str, Path] = {"corrected_volume": args.output}
        print(f"Wrote {args.output} (RAS)")
        if result.atlas_refined is not None and args.atlas_refined_output is not None:
            _save_nifti_gz(result.atlas_refined, args.atlas_refined_output)
            outputs["atlas_refined_volume"] = args.atlas_refined_output
            print(f"Wrote {args.atlas_refined_output} (RAS)")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={
                "volume": args.volume,
                "paint_mask": args.paint_mask,
                "multi_atlas_wm": args.multi_atlas_wm,
            },
            parameters={
                **parameters,
                "atlas_refinement_enabled": args.multi_atlas_wm is not None,
            },
            outputs=outputs,
        )
    elif args.command == "segmentation-refine":
        _require_nifti_gz_output(args.output)
        if args.summary.suffix.lower() != ".json":
            parser.error("Segmentation refinement summary must use the .json extension")
        refinement_parameters = {
            "inner_coefficient_mm": args.inner_coefficient_mm,
            "outer_coefficient_mm": args.outer_coefficient_mm,
            "inner_central_coefficient_mm": args.inner_central_coefficient_mm,
            "outer_central_coefficient_mm": args.outer_central_coefficient_mm,
            "deal_with_central": args.deal_with_central,
            "neighborhood_radius": args.neighborhood_radius,
            "max_sweeps": args.max_sweeps,
            "source_faithful": args.source_faithful,
        }
        result = refine_segmentation_by_surfaces(
            _load_nifti_gz(args.segmentation),
            _load_nifti_gz(args.inner_levelset),
            _load_nifti_gz(args.outer_levelset),
            **refinement_parameters,
        )
        _save_nifti_gz(result.image, args.output)
        rule_names = sorted(
            set(result.changed_voxels)
            | set(result.unresolved_voxels)
            | set(result.sweeps)
            | set(result.threshold_voxels)
            | set(result.candidate_masks)
        )
        summary = {
            "schema_version": 1,
            "command": args.command,
            "spatial_convention": "RAS",
            "source_faithful": args.source_faithful,
            "output": str(args.output.expanduser().resolve()),
            "changed_voxels": dict(result.changed_voxels),
            "unresolved_voxels": dict(result.unresolved_voxels),
            "sweeps": dict(result.sweeps),
            "threshold_voxels": dict(result.threshold_voxels),
            "rules": {
                name: {
                    "candidate_voxels": (
                        int(result.candidate_masks[name].sum())
                        if name in result.candidate_masks
                        else 0
                    ),
                    "changed_voxels": result.changed_voxels.get(name, 0),
                    "unresolved_voxels": result.unresolved_voxels.get(name, 0),
                    "sweeps": result.sweeps.get(name, 0),
                    "threshold_voxels": result.threshold_voxels.get(name),
                }
                for name in rule_names
            },
        }
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"Wrote {args.output} (RAS)")
        print(f"Wrote rule summary: {args.summary}")
        print(f"source_faithful={str(args.source_faithful).lower()}")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={
                "segmentation": args.segmentation,
                "inner_levelset": args.inner_levelset,
                "outer_levelset": args.outer_levelset,
            },
            parameters=refinement_parameters,
            outputs={"refined_segmentation": args.output, "rule_summary": args.summary},
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
    elif args.command == "n4":
        _require_nifti_gz_output(args.output)
        config = N4Config(
            shrink_factor=args.shrink_factor,
            iterations=tuple(args.iterations),
            convergence_tolerance=args.convergence_tolerance,
            spline_mesh=tuple(args.spline_mesh),
            rescale_intensities=args.rescale_intensities,
        )
        try:
            result = n4_bias_correct(_load_nifti_gz(args.input), config=config)
        except N4DependencyError as exc:
            parser.error(str(exc))
        _save_nifti_gz(result, args.output)
        print(f"Wrote {args.output} (RAS)")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"intensity_image": args.input},
            parameters=asdict(config),
            outputs={"corrected_image": args.output},
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
    elif args.command == "remove-cerebellum":
        maps = {
            "gm": args.gm,
            "wm_fill": args.wm_fill,
            "wm_mask": args.wm_mask,
        }
        segmentation = _load_nifti_gz(args.segmentation)
        output_names = {
            "gm": "gm_cerebellum_removed.nii.gz",
            "wm_fill": "wmfill_cerebellum_removed.nii.gz",
            "wm_mask": "wmmask_cerebellum_removed.nii.gz",
        }
        output_paths = {key: args.output_dir / name for key, name in output_names.items()}
        collisions = [path for path in output_paths.values() if path.exists()]
        if collisions:
            raise FileExistsError(
                "Refusing to overwrite cerebellum-removal output(s): "
                + ", ".join(str(path) for path in collisions)
            )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for key, input_path in maps.items():
            result = remove_brainstem_and_cerebellum(
                _load_nifti_gz(input_path), segmentation
            )
            _save_nifti_gz(result.image, output_paths[key])
            print(
                f"Wrote {output_paths[key]}; transferred brainstem/cerebellum "
                f"voxels={result.removed_voxels} (RAS)"
            )
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"segmentation": args.segmentation, **maps},
            parameters={"removed_labels": sorted((35, 38, 39, 40, 41, 71, 72, 73)), "interpolation": "nearest_neighbor"},
            outputs=output_paths,
        )
    elif args.command == "crop-tissue-maps":
        cropped = crop_toads_tissue_maps(
            _load_nifti_gz(args.gm),
            _load_nifti_gz(args.wm_fill),
            _load_nifti_gz(args.wm_mask),
            edge_width=args.edge_width,
        )
        output_paths = {
            "gm": args.output_dir / "gm_cropped.nii.gz",
            "wm_fill": args.output_dir / "wmfill_cropped.nii.gz",
            "wm_mask": args.output_dir / "wmmask_cropped.nii.gz",
        }
        collisions = [path for path in output_paths.values() if path.exists()]
        if collisions:
            raise FileExistsError(
                "Refusing to overwrite crop output(s): "
                + ", ".join(str(path) for path in collisions)
            )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        images = {
            "gm": cropped.gray_matter,
            "wm_fill": cropped.white_matter_fill,
            "wm_mask": cropped.white_matter_mask,
        }
        for name, image in images.items():
            _save_nifti_gz(image, output_paths[name])
            print(f"Wrote {output_paths[name]} (RAS; bounds={cropped.crop_parameters.bounds})")
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={"gm": args.gm, "wm_fill": args.wm_fill, "wm_mask": args.wm_mask},
            parameters={
                "threshold": 0.0,
                "edge_width": args.edge_width,
                "bounds": [list(axis) for axis in cropped.crop_parameters.bounds],
            },
            outputs=output_paths,
        )
    elif args.command == "prepare-cruise-inputs":
        prepared = prepare_toads_cruise_inputs(
            _load_nifti_gz(args.gm),
            _load_nifti_gz(args.wm_fill),
            _load_nifti_gz(args.wm_mask),
            segmentation=(
                None if args.segmentation is None else _load_nifti_gz(args.segmentation)
            ),
            remove_cerebellum=args.remove_cerebellum,
            crop_edge_width=args.edge_width,
        )
        output_paths = {
            "gm_crop": args.output_dir / "gm_cropped.nii.gz",
            "wm_fill_crop": args.output_dir / "wmfill_cropped.nii.gz",
            "wm_mask_crop": args.output_dir / "wmmask_cropped.nii.gz",
            "wm_fill_topology_corrected": args.output_dir / "wmfill_topology_corrected.nii.gz",
        }
        collisions = [path for path in output_paths.values() if path.exists()]
        if collisions:
            raise FileExistsError(
                "Refusing to overwrite prepared CRUISE input(s): "
                + ", ".join(str(path) for path in collisions)
            )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        images = {
            "gm_crop": prepared.gray_matter,
            "wm_fill_crop": prepared.white_matter_fill,
            "wm_mask_crop": prepared.white_matter_mask,
            "wm_fill_topology_corrected": prepared.topology_corrected_white_matter_fill,
        }
        for name, image in images.items():
            _save_nifti_gz(image, output_paths[name])
            print(f"Wrote {output_paths[name]} (RAS)")
        parameters = {
            "remove_cerebellum": args.remove_cerebellum,
            "removed_voxel_counts": prepared.removed_voxel_counts,
            "crop_edge_width": args.edge_width,
            "crop_bounds": [list(axis) for axis in prepared.crop_parameters.bounds],
            "topology": {
                "foreground_background_connectivity": "(6,18)",
                "paint_threshold": 0.5,
                "lowest_level": 0.0,
                "highest_level": 1.0,
                "regularize": 1e-4,
                "multi_atlas_refinement": False,
            },
            "stage_order": prepared.stage_order,
        }
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={
                "gm": args.gm,
                "wm_fill": args.wm_fill,
                "wm_mask": args.wm_mask,
                "segmentation": args.segmentation,
            },
            parameters=parameters,
            outputs=output_paths,
        )
    elif args.command == "toads-cruise-core":
        core_parameters = TGDMBaselineParameters(
            white_matter_iso=0.6,
            initial_iterations=args.tgdm_initial_iterations,
            inner_iterations=args.tgdm_inner_iterations,
            central_iterations=args.tgdm_central_iterations,
            # The selected MATLAB command wrapper passes IterCentralSurf to
            # the Java outer-iteration parameter as well.
            outer_iterations=args.tgdm_central_iterations,
        )
        result = run_toads_cruise_core(
            _load_nifti_gz(args.gm),
            _load_nifti_gz(args.wm_fill),
            _load_nifti_gz(args.wm_mask),
            segmentation=(
                None if args.segmentation is None else _load_nifti_gz(args.segmentation)
            ),
            remove_cerebellum=args.remove_cerebellum,
            crop_edge_width=args.edge_width,
            ace_diffusion_iterations=args.ace_diffusion_iterations,
            gvf_lambda=args.gvf_lambda,
            gvf_iterations=args.gvf_iterations,
            gvf_solver_tolerance=args.gvf_solver_tolerance,
            tgdm_parameters=core_parameters,
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        images = {
            "gm_cropped": result.preparation.gray_matter,
            "wmfill_cropped": result.preparation.white_matter_fill,
            "wmmask_cropped": result.preparation.white_matter_mask,
            "wmfill_topology_corrected": result.preparation.topology_corrected_white_matter_fill,
            "ace_enhanced_gm": result.ace.enhanced_gray_matter,
            "ace_skeleton": result.ace.skeleton,
            "ace_thinned_skeleton": result.ace.thinned_skeleton,
            "gvf": result.gvf.vector_field,
            "inner": result.inner_native,
            "central": result.central_native,
            "outer": result.outer_native,
        }
        output_paths = {name: args.output_dir / f"{name}.nii.gz" for name in images}
        collisions = [path for path in output_paths.values() if path.exists()]
        if collisions:
            raise FileExistsError(
                "Refusing to overwrite CRUISE core output(s): "
                + ", ".join(str(path) for path in collisions)
            )
        for name, image in images.items():
            _save_nifti_gz(image, output_paths[name])
            print(f"Wrote {output_paths[name]} (RAS)")
        core_options = {
            "remove_cerebellum": args.remove_cerebellum,
            "removed_voxel_counts": result.preparation.removed_voxel_counts,
            "crop_edge_width": args.edge_width,
            "crop_bounds": [list(axis) for axis in result.preparation.crop_parameters.bounds],
            "ace": {
                "intensity_threshold": 126.9,
                "tensor_alpha": 1.6,
                "gaussian_alpha": 1.8,
                "diffusion_iterations": args.ace_diffusion_iterations,
                "updated_multi_atlas_csf_mode": 0,
            },
            "gvf": {
                "lambda": args.gvf_lambda,
                "max_iterations": args.gvf_iterations,
                "solver_tolerance": args.gvf_solver_tolerance,
            },
            "tgdm": asdict(core_parameters),
            "uncrop_background": 20.0,
            "stage_order": [
                *result.preparation.stage_order,
                "ACE",
                "generalizedGVF",
                "nestedTGDM",
                "inverseCrop",
            ],
            "visualization_surfaces": False,
        }
        _write_provenance(
            args.provenance,
            command=args.command,
            inputs={
                "gm": args.gm,
                "wm_fill": args.wm_fill,
                "wm_mask": args.wm_mask,
                "segmentation": args.segmentation,
            },
            parameters=core_options,
            outputs=output_paths,
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
