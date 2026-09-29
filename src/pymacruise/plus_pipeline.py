"""Selected MaCRUISEplus entry-point options, separate from base defaults.

This module records the option values passed by the MATLAB plus entry point.
It does not implement ACE, topology correction, CRUISE, or LesionTOADS stages.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import nibabel as nib

from pymacruise.ace import ACEImageResult, ace_enhance
from pymacruise.boundary import clean_csf_boundary_image, multi_atlas_csf_boundary_image
from pymacruise.gvf import GVFImageResult, generalized_gvf_image
from pymacruise.tgdm import TGDMBaselineParameters
from pymacruise.topology import TopologyCorrectionResult, topology_correct


@dataclass(frozen=True)
class PlusPipelineOptions:
    """Configuration overrides selected by ``maCRUISEpluspipeline.m``."""

    lesion_toads_output: str = "dura removal inputs"
    lesion_toads_correct: bool = False
    wm_iso_level: float = 0.5
    gvf_lambda: float = 0.2
    update_enhance: int = 0
    generate_intermediate: bool = True
    use_updated_csf_in_ace: int = 3
    boundary_erosion_method: str = "3D"
    boundary_overlook_labels: tuple[int, ...] = (120, 121, 202, 203, 156, 157)
    refine_topology_by_maseg: bool = True

    def __post_init__(self) -> None:
        if self.lesion_toads_output != "dura removal inputs":
            raise ValueError("selected plus pipeline requires LesionTOADS dura removal inputs")
        if not isinstance(self.lesion_toads_correct, bool):
            raise ValueError("lesion_toads_correct must be boolean")
        if not 0.0 < self.wm_iso_level < 1.0:
            raise ValueError("wm_iso_level must be strictly between zero and one")
        if self.gvf_lambda < 0.0:
            raise ValueError("gvf_lambda must be nonnegative")
        if self.update_enhance not in (0, 1):
            raise ValueError("update_enhance must be 0 or 1")
        if not isinstance(self.generate_intermediate, bool):
            raise ValueError("generate_intermediate must be boolean")
        if self.use_updated_csf_in_ace not in (0, 1, 2, 3):
            raise ValueError("use_updated_csf_in_ace must be one of 0, 1, 2, or 3")
        if self.boundary_erosion_method not in ("3D", "3Dmatlab"):
            raise ValueError("unsupported boundary erosion method")
        if any(not isinstance(label, int) for label in self.boundary_overlook_labels):
            raise ValueError("boundary_overlook_labels must contain integers")
        if not isinstance(self.refine_topology_by_maseg, bool):
            raise ValueError("refine_topology_by_maseg must be boolean")

    def as_dict(self) -> dict[str, Any]:
        """Return JSON-ready Python names/values for run provenance."""
        result = asdict(self)
        result["boundary_overlook_labels"] = list(self.boundary_overlook_labels)
        return result


@dataclass(frozen=True)
class PlusBoundaryResult:
    """Raw and TOADS-WM-cleaned RAS boundary products for MaCRUISEplus."""

    raw_boundary: nib.spatialimages.SpatialImage
    cleaned_boundary: nib.spatialimages.SpatialImage


def run_plus_gvf_stage(
    scalar_volume: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> GVFImageResult:
    """Run the existing RAS GVF adapter with selected MaCRUISEplus settings.

    The source wrapper passes ``GVF_Lambda`` and fixes maximum iterations at
    200 with output normalization disabled. Solver implementation remains in
    :mod:`pymacruise.gvf`; this function only maps plus configuration to it.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    return generalized_gvf_image(
        scalar_volume,
        lambda_value=options.gvf_lambda,
        max_iterations=200,
    )


def plus_tgdm_parameters(
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> TGDMBaselineParameters:
    """Map the plus white-matter isolevel into the existing TGDM settings.

    The MATLAB plus runner overrides `WMIsoLevel` to 0.5 before calling its
    surface pipeline. All other selected TGDM values retain the baseline
    Python defaults; this adapter does not implement surface evolution.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    return TGDMBaselineParameters(
        white_matter_iso=options.wm_iso_level,
    )


def run_plus_ace_stage(
    gray_matter: nib.spatialimages.SpatialImage,
    white_matter: nib.spatialimages.SpatialImage,
    multi_atlas_csf: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> ACEImageResult:
    """Run the selected MaCRUISEplus ACE mode from typed plus settings.

    The selected entry point uses a positive-CSF merge after the speed
    Gaussian (`UseUpdatedCSFinACE=3`) and disables the later enhanced-GM /
    skeleton update branch. The boundary input is supplied on the shared grid.
    """
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    if options.update_enhance != 0:
        raise ValueError("plus ACE currently supports only update_enhance=0")
    if options.use_updated_csf_in_ace != 3:
        raise ValueError("selected plus ACE requires use_updated_csf_in_ace=3")
    return ace_enhance(
        gray_matter,
        white_matter,
        multi_atlas_csf=multi_atlas_csf,
        use_updated_csf_in_ace=options.use_updated_csf_in_ace,
        merge_method="Merge",
    )


def run_plus_topology_stage(
    volume: nib.spatialimages.SpatialImage,
    paint_mask: nib.spatialimages.SpatialImage,
    ma_wm: nib.spatialimages.SpatialImage | None,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> TopologyCorrectionResult:
    """Run selected plus topology correction and optional MaWM refinement."""
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    if options.refine_topology_by_maseg and ma_wm is None:
        raise ValueError("plus topology refinement requires the selected MaWM input")
    return topology_correct(
        volume,
        paint_mask,
        multi_atlas_wm=ma_wm if options.refine_topology_by_maseg else None,
    )


def run_plus_boundary_stage(
    segmentation: nib.spatialimages.SpatialImage,
    toads_wm_mask: nib.spatialimages.SpatialImage,
    options: PlusPipelineOptions = PlusPipelineOptions(),
) -> PlusBoundaryResult:
    """Create the selected plus 3D boundary and clean it with TOADS WM."""
    if not isinstance(options, PlusPipelineOptions):
        raise ValueError("options must be PlusPipelineOptions")
    if options.boundary_erosion_method != "3D":
        raise ValueError("selected plus boundary stage supports only ErodeMethod='3D'")
    raw = multi_atlas_csf_boundary_image(
        segmentation,
        overlooked_labels=options.boundary_overlook_labels,
    )
    cleaned = clean_csf_boundary_image(raw, toads_wm_mask, mode="toads_wm")
    return PlusBoundaryResult(raw_boundary=raw, cleaned_boundary=cleaned)
