# Complete command reference

Generated from the current command parser with `python scripts/generate_cli_reference.py`. All commands below run independently; supply the required inputs shown in their usage. Use new output directories, retain provenance, and keep images on the physical grids required by [the detailed I/O guide](cli-help.md).

For full PLUS testing and a stage-by-stage recipe, start with [the PLUS walkthrough](testing/full-plus-walkthrough.md) and [the individual tools guide](testing/individual-tools.md). Ordinary argparse defaults are listed below; the PLUS testing kit deliberately overrides initial atlas alignment to 50 iterations/four levels.

## `prepare-raw`

```text
usage: pymacruise prepare-raw [-h] --t1 T1 --maseg MASEG
                              --brain-mask BRAIN_MASK --output-dir OUTPUT_DIR
                              [--overwrite]

options:
  -h, --help            show this help message and exit
  --t1 T1
  --maseg MASEG
  --brain-mask BRAIN_MASK
  --output-dir OUTPUT_DIR
  --overwrite
```

Defaults: `--overwrite` = `False`.

## `cortical-separation`

```text
usage: pymacruise cortical-separation [-h] --maseg MASEG
                                      --output-dir OUTPUT_DIR

options:
  -h, --help            show this help message and exit
  --maseg MASEG
  --output-dir OUTPUT_DIR
```

## `calculate`

```text
usage: pymacruise calculate [-h] --input-a INPUT_A --input-b INPUT_B
                            --operation OPERATION --output OUTPUT
                            [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --input-a INPUT_A
  --input-b INPUT_B
  --operation OPERATION
  --output OUTPUT
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--provenance` = `None`.

## `topology-correct`

```text
usage: pymacruise topology-correct [-h] --volume VOLUME
                                   --paint-mask PAINT_MASK
                                   [--multi-atlas-wm MULTI_ATLAS_WM]
                                   --output OUTPUT
                                   [--atlas-refined-output ATLAS_REFINED_OUTPUT]
                                   [--paint-threshold PAINT_THRESHOLD]
                                   [--lowest-level LOWEST_LEVEL]
                                   [--highest-level HIGHEST_LEVEL]
                                   [--regularize REGULARIZE]
                                   [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --volume VOLUME       scalar volume .nii.gz
  --paint-mask PAINT_MASK
                        object paint mask .nii.gz
  --multi-atlas-wm MULTI_ATLAS_WM
                        optional multi-atlas WM mask .nii.gz
  --output OUTPUT       topology-corrected volume .nii.gz
  --atlas-refined-output ATLAS_REFINED_OUTPUT
                        optional V3 multi-atlas-refined volume .nii.gz;
                        requires --multi-atlas-wm
  --paint-threshold PAINT_THRESHOLD
  --lowest-level LOWEST_LEVEL
  --highest-level HIGHEST_LEVEL
  --regularize REGULARIZE
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--multi-atlas-wm` = `None`; `--atlas-refined-output` = `None`; `--paint-threshold` = `0.5`; `--lowest-level` = `0.0`; `--highest-level` = `1.0`; `--regularize` = `0.0001`; `--provenance` = `None`.

## `segmentation-refine`

```text
usage: pymacruise segmentation-refine [-h] --segmentation SEGMENTATION
                                      --inner-levelset INNER_LEVELSET
                                      --outer-levelset OUTER_LEVELSET
                                      --output OUTPUT --summary SUMMARY
                                      [--inner-coefficient-mm INNER_COEFFICIENT_MM]
                                      [--outer-coefficient-mm OUTER_COEFFICIENT_MM]
                                      [--inner-central-coefficient-mm INNER_CENTRAL_COEFFICIENT_MM]
                                      [--outer-central-coefficient-mm OUTER_CENTRAL_COEFFICIENT_MM]
                                      [--deal-with-central | --no-deal-with-central]
                                      [--neighborhood-radius NEIGHBORHOOD_RADIUS]
                                      [--max-sweeps MAX_SWEEPS]
                                      [--source-faithful | --no-source-faithful]
                                      [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --segmentation SEGMENTATION
                        MaSeg labels .nii.gz
  --inner-levelset INNER_LEVELSET
                        inner level set .nii.gz
  --outer-levelset OUTER_LEVELSET
                        outer level set .nii.gz
  --output OUTPUT       refined labels .nii.gz
  --summary SUMMARY     JSON per-rule result summary
  --inner-coefficient-mm INNER_COEFFICIENT_MM
  --outer-coefficient-mm OUTER_COEFFICIENT_MM
  --inner-central-coefficient-mm INNER_CENTRAL_COEFFICIENT_MM
  --outer-central-coefficient-mm OUTER_CENTRAL_COEFFICIENT_MM
  --deal-with-central, --no-deal-with-central
                        enable central-label handling (default: enabled)
  --neighborhood-radius NEIGHBORHOOD_RADIUS
  --max-sweeps MAX_SWEEPS
  --source-faithful, --no-source-faithful
                        match MATLAB positive-integer tabulate/maxind labels
                        (default: enabled; both modes agree)
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--inner-coefficient-mm` = `0.2`; `--outer-coefficient-mm` = `0.2`; `--inner-central-coefficient-mm` = `0.2`; `--outer-central-coefficient-mm` = `0.2`; `--deal-with-central` = `True`; `--neighborhood-radius` = `2`; `--max-sweeps` = `50`; `--source-faithful` = `True`; `--provenance` = `None`.

## `ace`

```text
usage: pymacruise ace [-h] --gray GRAY --white WHITE [--atlas-csf ATLAS_CSF]
                      --enhanced ENHANCED --skeleton SKELETON
                      --thinned-skeleton THINNED_SKELETON
                      [--diffusion-iterations DIFFUSION_ITERATIONS]
                      [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --gray GRAY
  --white WHITE
  --atlas-csf ATLAS_CSF
  --enhanced ENHANCED   enhanced GM .nii.gz
  --skeleton SKELETON   ACE skeleton .nii.gz
  --thinned-skeleton THINNED_SKELETON
                        thinned skeleton .nii.gz
  --diffusion-iterations DIFFUSION_ITERATIONS
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--atlas-csf` = `None`; `--diffusion-iterations` = `25`; `--provenance` = `None`.

## `gvf`

```text
usage: pymacruise gvf [-h] --input INPUT --output OUTPUT
                      [--lambda LAMBDA_VALUE]
                      [--max-iterations MAX_ITERATIONS]
                      [--solver-tolerance SOLVER_TOLERANCE]
                      [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --input INPUT
  --output OUTPUT       component-last vector .nii.gz
  --lambda LAMBDA_VALUE
  --max-iterations MAX_ITERATIONS
  --solver-tolerance SOLVER_TOLERANCE
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--lambda` = `1.0`; `--max-iterations` = `200`; `--solver-tolerance` = `1e-06`; `--provenance` = `None`.

## `n4`

```text
usage: pymacruise n4 [-h] --input INPUT --output OUTPUT
                     [--working-reference WORKING_REFERENCE]
                     [--orientation-only-working-grid]
                     [--shrink-factor SHRINK_FACTOR]
                     [--iterations ITERATIONS [ITERATIONS ...]]
                     [--convergence-tolerance CONVERGENCE_TOLERANCE]
                     [--spline-mesh SPLINE_MESH SPLINE_MESH SPLINE_MESH]
                     [--estimation-mask {all,legacy-otsu}]
                     [--rescale-intensities | --no-rescale-intensities]
                     [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --input INPUT
  --output OUTPUT       corrected image .nii.gz
  --working-reference WORKING_REFERENCE
                        map corrected native-grid data to this working grid
                        after N4
  --orientation-only-working-grid
                        map to the working grid by axis order only, without
                        interpolation
  --shrink-factor SHRINK_FACTOR
  --iterations ITERATIONS [ITERATIONS ...]
  --convergence-tolerance CONVERGENCE_TOLERANCE
  --spline-mesh SPLINE_MESH SPLINE_MESH SPLINE_MESH
  --estimation-mask {all,legacy-otsu}
                        N4 fit region; BASE selects historical ANTs 200-bin
                        automatic Otsu
  --rescale-intensities, --no-rescale-intensities
                        rescale corrected values to the input range (default:
                        disabled without an explicit mask)
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--working-reference` = `None`; `--orientation-only-working-grid` = `False`; `--shrink-factor` = `4`; `--iterations` = `[50, 50, 50, 50]`; `--convergence-tolerance` = `0.0`; `--spline-mesh` = `[1, 1, 1]`; `--estimation-mask` = `all`; `--rescale-intensities` = `False`; `--provenance` = `None`.

## `base-preprocessing`

```text
usage: pymacruise base-preprocessing [-h] --t1-n4 T1_N4
                                     --brain-mask BRAIN_MASK
                                     --output-dir OUTPUT_DIR
                                     [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --t1-n4 T1_N4
  --brain-mask BRAIN_MASK
  --output-dir OUTPUT_DIR
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--provenance` = `None`.

## `tgdm`

```text
usage: pymacruise tgdm [-h] --initial INITIAL --white WHITE --gray GRAY
                       --ace-gray ACE_GRAY [--brain-mask BRAIN_MASK]
                       [--brain-mask-level-set BRAIN_MASK_LEVEL_SET]
                       [--ace-skeleton ACE_SKELETON] [--gvf GVF]
                       --output-dir OUTPUT_DIR [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --initial INITIAL
  --white WHITE
  --gray GRAY
  --ace-gray ACE_GRAY
  --brain-mask BRAIN_MASK
  --brain-mask-level-set BRAIN_MASK_LEVEL_SET
                        pass an existing signed Brain Mask Level Set with its
                        numeric magnitudes preserved
  --ace-skeleton ACE_SKELETON
  --gvf GVF
  --output-dir OUTPUT_DIR
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--brain-mask` = `None`; `--brain-mask-level-set` = `None`; `--ace-skeleton` = `None`; `--gvf` = `None`; `--provenance` = `None`.

## `dura`

```text
usage: pymacruise dura [-h] --t1 T1 --white WHITE --gray GRAY
                       --white-mask WHITE_MASK --output-dir OUTPUT_DIR
                       [--max-iterations MAX_ITERATIONS]
                       [--source-axis-order {RAS,LPS}]
                       [--membership-level MEMBERSHIP_LEVEL]
                       [--wm-topology-regularization WM_TOPOLOGY_REGULARIZATION]
                       [--gm-topology-regularization GM_TOPOLOGY_REGULARIZATION]
                       [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --t1 T1
  --white WHITE
  --gray GRAY
  --white-mask WHITE_MASK
  --output-dir OUTPUT_DIR
  --max-iterations MAX_ITERATIONS
  --source-axis-order {RAS,LPS}
                        internal Dura voxel order; BASE selects the source LPS
                        order
  --membership-level MEMBERSHIP_LEVEL
                        target membership (MaCRUISE run_RemoveDura_jist.m
                        selects 0.3)
  --wm-topology-regularization WM_TOPOLOGY_REGULARIZATION
  --gm-topology-regularization GM_TOPOLOGY_REGULARIZATION
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--max-iterations` = `5`; `--source-axis-order` = `RAS`; `--membership-level` = `0.3`; `--wm-topology-regularization` = `0.05`; `--gm-topology-regularization` = `0.02`; `--provenance` = `None`.

## `atlasem-prelesion-dura`

```text
usage: pymacruise atlasem-prelesion-dura [-h] --t1 T1
                                         --memberships MEMBERSHIPS
                                         --brain-mask BRAIN_MASK
                                         --class-names CLASS_NAMES
                                         --output-dir OUTPUT_DIR
                                         [--max-iterations MAX_ITERATIONS]
                                         [--provenance PROVENANCE]

Consume an existing channels-last AtlasEM membership NIfTI and a matching
T1/mask. This runs AtlasEM output generation and the Dura stage; it does not
load an atlas, estimate registration, or run the AtlasEM EM loop.

options:
  -h, --help            show this help message and exit
  --t1 T1
  --memberships MEMBERSHIPS
                        existing float membership image with classes on the
                        last axis
  --brain-mask BRAIN_MASK
                        binary mask on the T1 physical grid
  --class-names CLASS_NAMES
                        JSON array in membership channel order, e.g.
                        '["BG","CSF","GM","WM"]'
  --output-dir OUTPUT_DIR
  --max-iterations MAX_ITERATIONS
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--max-iterations` = `5`; `--provenance` = `None`.

## `atlasem-atlas-inventory`

```text
usage: pymacruise atlasem-atlas-inventory [-h] --descriptor DESCRIPTOR

Read descriptor metadata only. Report class order, required payload paths, and
whether payload files are available as .nii.gz. Voxel data are not loaded.

options:
  -h, --help            show this help message and exit
  --descriptor DESCRIPTOR
```

## `atlasem-run`

```text
usage: pymacruise atlasem-run [-h] --t1 T1 --brain-mask BRAIN_MASK
                              --descriptor DESCRIPTOR --output-dir OUTPUT_DIR
                              [--modality MODALITY]
                              [--background-ratio BACKGROUND_RATIO]
                              [--crop-border CROP_BORDER] [--levels LEVELS]
                              [--initial-registration-iterations INITIAL_REGISTRATION_ITERATIONS]
                              [--post-recrop-registration-iterations POST_RECROP_REGISTRATION_ITERATIONS]
                              [--max-em-iterations MAX_EM_ITERATIONS]
                              [--max-membership-change MAX_MEMBERSHIP_CHANGE]
                              [--smoothing SMOOTHING]
                              [--max-dura-iterations MAX_DURA_ITERATIONS]
                              [--resume] [--overwrite]
                              [--provenance PROVENANCE]

Run the supported single-channel AtlasEM workflow from a T1 .nii.gz, matching
brain mask, and descriptor whose referenced atlas payloads are compressed
NIfTI files. This produces native-grid RAS AtlasEM and pre-LesionTOADS Dura
products; direct optimizer parity is not established.

options:
  -h, --help            show this help message and exit
  --t1 T1
  --brain-mask BRAIN_MASK
  --descriptor DESCRIPTOR
  --output-dir OUTPUT_DIR
  --modality MODALITY
  --background-ratio BACKGROUND_RATIO
  --crop-border CROP_BORDER
  --levels LEVELS
  --initial-registration-iterations INITIAL_REGISTRATION_ITERATIONS
  --post-recrop-registration-iterations POST_RECROP_REGISTRATION_ITERATIONS
  --max-em-iterations MAX_EM_ITERATIONS
  --max-membership-change MAX_MEMBERSHIP_CHANGE
  --smoothing SMOOTHING
  --max-dura-iterations MAX_DURA_ITERATIONS
  --resume
  --overwrite
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--modality` = `T1_MPRAGE`; `--background-ratio` = `0.0`; `--crop-border` = `10`; `--levels` = `4`; `--initial-registration-iterations` = `50`; `--post-recrop-registration-iterations` = `1`; `--max-em-iterations` = `99`; `--max-membership-change` = `0.01`; `--smoothing` = `10.0`; `--max-dura-iterations` = `5`; `--resume` = `False`; `--overwrite` = `False`; `--provenance` = `None`.

## `remove-cerebellum`

```text
usage: pymacruise remove-cerebellum [-h] --segmentation SEGMENTATION --gm GM
                                    --wm-fill WM_FILL --wm-mask WM_MASK
                                    --output-dir OUTPUT_DIR
                                    [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --segmentation SEGMENTATION
                        MaSeg label map .nii.gz
  --gm GM               second-pass TOADS GM map .nii.gz
  --wm-fill WM_FILL     second-pass TOADS WM-fill map .nii.gz
  --wm-mask WM_MASK     second-pass TOADS WM-mask map .nii.gz
  --output-dir OUTPUT_DIR
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--provenance` = `None`.

## `crop-tissue-maps`

```text
usage: pymacruise crop-tissue-maps [-h] --gm GM --wm-fill WM_FILL
                                   --wm-mask WM_MASK --output-dir OUTPUT_DIR
                                   [--edge-width EDGE_WIDTH]
                                   [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --gm GM
  --wm-fill WM_FILL
  --wm-mask WM_MASK
  --output-dir OUTPUT_DIR
  --edge-width EDGE_WIDTH
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--edge-width` = `5`; `--provenance` = `None`.

## `prepare-cruise-inputs`

```text
usage: pymacruise prepare-cruise-inputs [-h] --gm GM --wm-fill WM_FILL
                                        --wm-mask WM_MASK
                                        [--segmentation SEGMENTATION]
                                        [--remove-cerebellum]
                                        [--edge-width EDGE_WIDTH]
                                        --output-dir OUTPUT_DIR
                                        [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --gm GM
  --wm-fill WM_FILL
  --wm-mask WM_MASK
  --segmentation SEGMENTATION
                        required with --remove-cerebellum
  --remove-cerebellum
  --edge-width EDGE_WIDTH
  --output-dir OUTPUT_DIR
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--segmentation` = `None`; `--remove-cerebellum` = `False`; `--edge-width` = `5`; `--provenance` = `None`.

## `base-tissue-preparation`

```text
usage: pymacruise base-tissue-preparation [-h] --maseg MASEG --gm GM
                                          --wm-fill WM_FILL --wm-mask WM_MASK
                                          [--remove-cerebellum]
                                          [--separate-hemispheres]
                                          [--edge-width EDGE_WIDTH]
                                          [--gm-confident-level GM_CONFIDENT_LEVEL]
                                          --output-dir OUTPUT_DIR
                                          --provenance PROVENANCE

options:
  -h, --help            show this help message and exit
  --maseg MASEG         MaSeg labels and reference grid .nii.gz
  --gm GM               second-pass TOADS GM .nii.gz
  --wm-fill WM_FILL     second-pass TOADS WM-fill .nii.gz
  --wm-mask WM_MASK     second-pass TOADS WM-mask .nii.gz
  --remove-cerebellum
  --separate-hemispheres
                        apply the v3.2 left/right cortical boundary exclusion
                        (selected by base-run)
  --edge-width EDGE_WIDTH
  --gm-confident-level GM_CONFIDENT_LEVEL
  --output-dir OUTPUT_DIR
  --provenance PROVENANCE
```

Defaults: `--remove-cerebellum` = `False`; `--separate-hemispheres` = `False`; `--edge-width` = `5`; `--gm-confident-level` = `0.001`.

## `toads-cruise-core`

```text
usage: pymacruise toads-cruise-core [-h] --gm GM --wm-fill WM_FILL
                                    --wm-mask WM_MASK
                                    [--brain-mask BRAIN_MASK]
                                    [--brain-mask-level-set BRAIN_MASK_LEVEL_SET]
                                    [--boundary-csf BOUNDARY_CSF]
                                    [--segmentation SEGMENTATION]
                                    [--remove-cerebellum]
                                    [--base-preparation-provenance BASE_PREPARATION_PROVENANCE]
                                    [--edge-width EDGE_WIDTH]
                                    [--topology-propagation {object->background,background->object}]
                                    [--topology-source-axis-order {RAS,LPS}]
                                    [--ace-diffusion-iterations ACE_DIFFUSION_ITERATIONS]
                                    [--use-updated-csf-in-ace {0,1,2,3}]
                                    [--wm-iso-level WM_ISO_LEVEL]
                                    [--gvf-lambda GVF_LAMBDA]
                                    [--gvf-iterations GVF_ITERATIONS]
                                    [--gvf-solver-tolerance GVF_SOLVER_TOLERANCE]
                                    [--tgdm-initial-iterations TGDM_INITIAL_ITERATIONS]
                                    [--tgdm-inner-iterations TGDM_INNER_ITERATIONS]
                                    [--tgdm-central-iterations TGDM_CENTRAL_ITERATIONS]
                                    --output-dir OUTPUT_DIR
                                    [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --gm GM
  --wm-fill WM_FILL
  --wm-mask WM_MASK
  --brain-mask BRAIN_MASK
                        user-provided binary brain mask on the T1 or prepared
                        crop grid
  --brain-mask-level-set BRAIN_MASK_LEVEL_SET
                        signed Dura Brain Mask Level Set; preserves values on
                        the T1 or crop grid
  --boundary-csf BOUNDARY_CSF
                        V3 MaSeg CSF boundary cleaned by TOADS and updated WM
                        for ACE mode 2
  --segmentation SEGMENTATION
  --remove-cerebellum
  --base-preparation-provenance BASE_PREPARATION_PROVENANCE
                        reuse the recorded MaSeg-derived crop for already-
                        fused base maps
  --edge-width EDGE_WIDTH
  --topology-propagation {object->background,background->object}
  --topology-source-axis-order {RAS,LPS}
  --ace-diffusion-iterations ACE_DIFFUSION_ITERATIONS
  --use-updated-csf-in-ace {0,1,2,3}
                        ACE updated-CSF placement mode (selected MaCRUISE V3
                        profile: 2)
  --wm-iso-level WM_ISO_LEVEL
  --gvf-lambda GVF_LAMBDA
  --gvf-iterations GVF_ITERATIONS
  --gvf-solver-tolerance GVF_SOLVER_TOLERANCE
  --tgdm-initial-iterations TGDM_INITIAL_ITERATIONS
  --tgdm-inner-iterations TGDM_INNER_ITERATIONS
  --tgdm-central-iterations TGDM_CENTRAL_ITERATIONS
  --output-dir OUTPUT_DIR
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--brain-mask` = `None`; `--brain-mask-level-set` = `None`; `--boundary-csf` = `None`; `--segmentation` = `None`; `--remove-cerebellum` = `False`; `--base-preparation-provenance` = `None`; `--edge-width` = `5`; `--topology-propagation` = `object->background`; `--topology-source-axis-order` = `RAS`; `--ace-diffusion-iterations` = `25`; `--use-updated-csf-in-ace` = `0`; `--wm-iso-level` = `0.6`; `--gvf-lambda` = `1.0`; `--gvf-iterations` = `200`; `--gvf-solver-tolerance` = `1e-06`; `--tgdm-initial-iterations` = `6`; `--tgdm-inner-iterations` = `4`; `--tgdm-central-iterations` = `7`; `--provenance` = `None`.

## `toads-loop`

```text
usage: pymacruise toads-loop [-h] --t1-n4 T1_N4 --brain-mask BRAIN_MASK
                             --output OUTPUT
                             [--dura-output-dir DURA_OUTPUT_DIR]
                             [--atlas-descriptor ATLAS_DESCRIPTOR]
                             [--modality MODALITY]
                             [--registration-levels REGISTRATION_LEVELS]
                             [--registration-iterations REGISTRATION_ITERATIONS]
                             [--post-refresh-registration-iterations POST_REFRESH_REGISTRATION_ITERATIONS]
                             [--registration-max-iterations REGISTRATION_MAX_ITERATIONS]
                             [--source-faithful-nonfinite-registration-trials]
                             [--padding PADDING]
                             [--maximum-iterations MAXIMUM_ITERATIONS]
                             [--convergence-threshold CONVERGENCE_THRESHOLD]
                             [--relation-limit RELATION_LIMIT]
                             [--maximum-distance MAXIMUM_DISTANCE]
                             [--preserve-input-axis-order |
                             --source-axis-order {RAS,LPS}]
                             [--category-exponents CATEGORY_EXPONENTS CATEGORY_EXPONENTS CATEGORY_EXPONENTS]
                             [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --t1-n4 T1_N4
  --brain-mask BRAIN_MASK
  --output OUTPUT       native-grid classification .nii.gz
  --dura-output-dir DURA_OUTPUT_DIR
                        optional directory for the four source-scaled Dura-
                        removal inputs
  --atlas-descriptor ATLAS_DESCRIPTOR
  --modality MODALITY
  --registration-levels REGISTRATION_LEVELS
  --registration-iterations REGISTRATION_ITERATIONS
  --post-refresh-registration-iterations POST_REFRESH_REGISTRATION_ITERATIONS
                        rigid atlas updates after template refresh (selected
                        V3 wrapper default: 5)
  --registration-max-iterations REGISTRATION_MAX_ITERATIONS
  --source-faithful-nonfinite-registration-trials
                        allow Java-style edge interpolation for NaN
                        registration trials; degenerate inputs may produce a
                        NaN transform and fail at the next pyramid scale
  --padding PADDING
  --maximum-iterations MAXIMUM_ITERATIONS
  --convergence-threshold CONVERGENCE_THRESHOLD
  --relation-limit RELATION_LIMIT
  --maximum-distance MAXIMUM_DISTANCE
  --preserve-input-axis-order
                        source-compatibility mode: validate/warn and normalize
                        the input to RAS, reconstruct its axial RAS/LPS voxel
                        order internally, then write every image in RAS
  --source-axis-order {RAS,LPS}
                        process in an explicit source voxel order internally
                        and export RAS images
  --category-exponents CATEGORY_EXPONENTS CATEGORY_EXPONENTS CATEGORY_EXPONENTS
                        override atlas CSF/GM/WM optimized factors (default:
                        use atlas values)
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--dura-output-dir` = `None`; `--atlas-descriptor` = `None`; `--modality` = `T1_MPRAGE`; `--registration-levels` = `4`; `--registration-iterations` = `50`; `--post-refresh-registration-iterations` = `5`; `--registration-max-iterations` = `5`; `--source-faithful-nonfinite-registration-trials` = `False`; `--padding` = `10`; `--maximum-iterations` = `50`; `--convergence-threshold` = `0.0001`; `--relation-limit` = `3`; `--maximum-distance` = `0.0001`; `--preserve-input-axis-order` = `False`; `--source-axis-order` = `None`; `--category-exponents` = `None`; `--provenance` = `None`.

## `surface-extract`

```text
usage: pymacruise surface-extract [-h] --level-set LEVEL_SET --output OUTPUT
                                  [--level LEVEL] [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --level-set LEVEL_SET
  --output OUTPUT       legacy VTK .vtk or GIFTI .surf.gii output
  --level LEVEL
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--level` = `0.0`; `--provenance` = `None`.

## `surface-convert`

```text
usage: pymacruise surface-convert [-h] --input INPUT --output OUTPUT
                                  [--input-frame {RAS,LPS}]
                                  [--point-data POINT_DATA]
                                  [--data-output DATA_OUTPUT]
                                  [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --input INPUT
  --output OUTPUT       .vtk or .surf.gii
  --input-frame {RAS,LPS}
                        required for untagged VTK inputs
  --point-data POINT_DATA
                        attach dense GIFTI data with matching vertex order
                        (repeatable)
  --data-output DATA_OUTPUT
                        also export vertex data without geometry as .shape.gii
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--input-frame` = `None`; `--point-data` = `[]`; `--data-output` = `None`; `--provenance` = `None`.

## `surface-info`

```text
usage: pymacruise surface-info [-h] --input INPUT [--input-frame {RAS,LPS}]
                               [--output OUTPUT]

options:
  -h, --help            show this help message and exit
  --input INPUT
  --input-frame {RAS,LPS}
  --output OUTPUT       optional JSON file (also printed to stdout)
```

Defaults: `--input-frame` = `None`; `--output` = `None`.

## `thickness-levelset`

```text
usage: pymacruise thickness-levelset [-h] --inner INNER --outer OUTER
                                     --output OUTPUT [--iso-level ISO_LEVEL]
                                     [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --inner INNER
  --outer OUTER
  --output OUTPUT       thickness .nii.gz
  --iso-level ISO_LEVEL
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--iso-level` = `0.0`; `--provenance` = `None`.

## `thickness-surface`

```text
usage: pymacruise thickness-surface [-h] --inner INNER --outer OUTER
                                    [--central CENTRAL]
                                    --inner-output INNER_OUTPUT
                                    [--central-output CENTRAL_OUTPUT]
                                    [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --inner INNER
  --outer OUTER
  --central CENTRAL
  --inner-output INNER_OUTPUT
  --central-output CENTRAL_OUTPUT
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--central` = `None`; `--central-output` = `None`; `--provenance` = `None`.

## `base-surface-thickness`

```text
usage: pymacruise base-surface-thickness [-h] --inner-levelset INNER_LEVELSET
                                         --central-levelset CENTRAL_LEVELSET
                                         --outer-levelset OUTER_LEVELSET
                                         --output-dir OUTPUT_DIR
                                         [--surface-format {vtk,gii}]
                                         [--provenance PROVENANCE]

options:
  -h, --help            show this help message and exit
  --inner-levelset INNER_LEVELSET
  --central-levelset CENTRAL_LEVELSET
  --outer-levelset OUTER_LEVELSET
  --output-dir OUTPUT_DIR
  --surface-format {vtk,gii}
  --provenance PROVENANCE
                        optional JSON sidecar recording inputs, effective
                        parameters, and outputs
```

Defaults: `--surface-format` = `vtk`; `--provenance` = `None`.

## `workflow-preflight`

```text
usage: pymacruise workflow-preflight [-h] --config CONFIG [--overwrite]

options:
  -h, --help       show this help message and exit
  --config CONFIG
  --overwrite      replace an existing run.json manifest only; do not remove
                   other files
```

Defaults: `--overwrite` = `False`.

## `workflow-run`

```text
usage: pymacruise workflow-run [-h] --config CONFIG [--resume]

options:
  -h, --help       show this help message and exit
  --config CONFIG
  --resume         resume from verified checkpoints in the configured output
                   directory
```

Defaults: `--resume` = `False`.

## `base-run`

```text
usage: pymacruise base-run [-h] --t1 T1 --brain-mask BRAIN_MASK
                           --atlas-descriptor ATLAS_DESCRIPTOR --maseg MASEG
                           --output-dir OUTPUT_DIR [--resume]
                           [--n4-estimation-mask {all,legacy-otsu}]
                           [--regrid-inputs | --no-regrid-inputs]
                           [--registration-levels REGISTRATION_LEVELS]
                           [--registration-iterations REGISTRATION_ITERATIONS]
                           [--post-refresh-registration-iterations POST_REFRESH_REGISTRATION_ITERATIONS]
                           [--registration-max-iterations REGISTRATION_MAX_ITERATIONS]
                           [--toads-iterations TOADS_ITERATIONS]
                           [--toads-axis-order {RAS,LPS}]
                           [--dura-iterations DURA_ITERATIONS]
                           [--dura-membership-level DURA_MEMBERSHIP_LEVEL]
                           [--crop-edge-width CROP_EDGE_WIDTH]
                           [--ace-diffusion-iterations ACE_DIFFUSION_ITERATIONS]
                           [--gvf-iterations GVF_ITERATIONS]
                           [--tgdm-iterations TGDM_ITERATIONS]
                           [--tgdm-initial-iterations TGDM_INITIAL_ITERATIONS]
                           [--tgdm-inner-iterations TGDM_INNER_ITERATIONS]
                           [--tgdm-central-iterations TGDM_CENTRAL_ITERATIONS]
                           [--wm-iso-level WM_ISO_LEVEL]
                           [--ace-use-updated-csf-in-ace {0,1,2,3}]
                           [--source-faithful-refinement | --no-source-faithful-refinement]

options:
  -h, --help            show this help message and exit
  --t1 T1               input T1 .nii.gz
  --brain-mask BRAIN_MASK
                        same-grid brain mask .nii.gz
  --atlas-descriptor ATLAS_DESCRIPTOR
  --maseg MASEG         same-grid MaSeg labels for cleanup and surface
                        refinement
  --output-dir OUTPUT_DIR
  --resume              resume from verified stage outputs in run.json
  --n4-estimation-mask {all,legacy-otsu}
                        N4 backend policy: all follows the current unmasked
                        ANTs command; legacy-otsu reproduces historical
                        automatic Otsu estimation
  --regrid-inputs, --no-regrid-inputs
                        resample oblique inputs to an orthogonal working grid
                        (default: disabled); source-compatible no-regrid
                        processing permutes/flips axes without interpolation;
                        final outputs restore the original grid and oblique
                        affine in both modes
  --registration-levels REGISTRATION_LEVELS
  --registration-iterations REGISTRATION_ITERATIONS
  --post-refresh-registration-iterations POST_REFRESH_REGISTRATION_ITERATIONS
  --registration-max-iterations REGISTRATION_MAX_ITERATIONS
  --toads-iterations TOADS_ITERATIONS
  --toads-axis-order {RAS,LPS}
                        internal TOADS voxel order (default: LPS, selected
                        source wrapper)
  --dura-iterations DURA_ITERATIONS
  --dura-membership-level DURA_MEMBERSHIP_LEVEL
  --crop-edge-width CROP_EDGE_WIDTH
  --ace-diffusion-iterations ACE_DIFFUSION_ITERATIONS
  --gvf-iterations GVF_ITERATIONS
  --tgdm-iterations TGDM_ITERATIONS
                        override all TGDM phase counts; individual phase
                        options take precedence
  --tgdm-initial-iterations TGDM_INITIAL_ITERATIONS
                        initial smoothing iterations (default: 6)
  --tgdm-inner-iterations TGDM_INNER_ITERATIONS
                        inner surface iterations (default: 4)
  --tgdm-central-iterations TGDM_CENTRAL_ITERATIONS
                        central/outer surface iterations (default: 7)
  --wm-iso-level WM_ISO_LEVEL
  --ace-use-updated-csf-in-ace {0,1,2,3}
  --source-faithful-refinement, --no-source-faithful-refinement
                        use historical MaSeg refinement behavior (default:
                        enabled)
```

Defaults: `--resume` = `False`; `--n4-estimation-mask` = `all`; `--regrid-inputs` = `False`; `--registration-levels` = `4`; `--registration-iterations` = `50`; `--post-refresh-registration-iterations` = `5`; `--registration-max-iterations` = `5`; `--toads-iterations` = `50`; `--toads-axis-order` = `LPS`; `--dura-iterations` = `5`; `--dura-membership-level` = `0.3`; `--crop-edge-width` = `5`; `--ace-diffusion-iterations` = `25`; `--gvf-iterations` = `200`; `--tgdm-iterations` = `None`; `--tgdm-initial-iterations` = `6`; `--tgdm-inner-iterations` = `4`; `--tgdm-central-iterations` = `7`; `--wm-iso-level` = `0.5`; `--ace-use-updated-csf-in-ace` = `2`; `--source-faithful-refinement` = `True`.

## `base-native-export`

```text
usage: pymacruise base-native-export [-h] --reference REFERENCE
                                     --working-reference WORKING_REFERENCE
                                     --segmentation SEGMENTATION
                                     --inner-levelset INNER_LEVELSET
                                     --central-levelset CENTRAL_LEVELSET
                                     --outer-levelset OUTER_LEVELSET
                                     --surface-dir SURFACE_DIR
                                     --output-dir OUTPUT_DIR
                                     [--orientation-only]

options:
  -h, --help            show this help message and exit
  --reference REFERENCE
  --working-reference WORKING_REFERENCE
  --segmentation SEGMENTATION
  --inner-levelset INNER_LEVELSET
  --central-levelset CENTRAL_LEVELSET
  --outer-levelset OUTER_LEVELSET
  --surface-dir SURFACE_DIR
  --output-dir OUTPUT_DIR
  --orientation-only    invert a recorded no-regrid axis-only working geometry
                        instead of physical resampling
```

Defaults: `--orientation-only` = `False`.

## `base-resume-from-toads`

```text
usage: pymacruise base-resume-from-toads [-h] --source-run SOURCE_RUN
                                         --maseg MASEG --output-dir OUTPUT_DIR
                                         [--resume]

options:
  -h, --help            show this help message and exit
  --source-run SOURCE_RUN
                        completed source run.json; read only
  --maseg MASEG         same-grid MaSeg labels .nii.gz
  --output-dir OUTPUT_DIR
                        distinct tail-run directory
  --resume              resume verified tail stages from the tail run.json
```

Defaults: `--resume` = `False`.

## `plus-register-flair`

```text
usage: pymacruise plus-register-flair [-h] --t1 T1 --flair FLAIR
                                      --output-dir OUTPUT_DIR [--resume]
                                      [--metric {mattes,meansquares,GC}]
                                      [--metric-bins METRIC_BINS]
                                      [--sampling-rate SAMPLING_RATE]
                                      [--iterations ITERATIONS [ITERATIONS ...]]
                                      [--shrink-factors SHRINK_FACTORS [SHRINK_FACTORS ...]]
                                      [--smoothing-sigmas SMOOTHING_SIGMAS [SMOOTHING_SIGMAS ...]]

options:
  -h, --help            show this help message and exit
  --t1 T1               fixed T1 reference .nii.gz
  --flair FLAIR         moving FLAIR .nii.gz
  --output-dir OUTPUT_DIR
  --resume
  --metric {mattes,meansquares,GC}
  --metric-bins METRIC_BINS
  --sampling-rate SAMPLING_RATE
  --iterations ITERATIONS [ITERATIONS ...]
  --shrink-factors SHRINK_FACTORS [SHRINK_FACTORS ...]
  --smoothing-sigmas SMOOTHING_SIGMAS [SMOOTHING_SIGMAS ...]
```

Defaults: `--resume` = `False`; `--metric` = `mattes`; `--metric-bins` = `32`; `--sampling-rate` = `0.2`; `--iterations` = `[1000, 500, 250, 0]`; `--shrink-factors` = `[6, 4, 2, 1]`; `--smoothing-sigmas` = `[3.0, 2.0, 1.0, 0.0]`.

## `plus-prepare-lesion`

```text
usage: pymacruise plus-prepare-lesion [-h] --t1 T1 --flair FLAIR
                                      --brain-mask BRAIN_MASK
                                      --output-dir OUTPUT_DIR [--resume]
                                      [--metric {mattes,meansquares,GC}]
                                      [--metric-bins METRIC_BINS]
                                      [--sampling-rate SAMPLING_RATE]
                                      [--registration-iterations REGISTRATION_ITERATIONS [REGISTRATION_ITERATIONS ...]]
                                      [--shrink-factors SHRINK_FACTORS [SHRINK_FACTORS ...]]
                                      [--smoothing-sigmas SMOOTHING_SIGMAS [SMOOTHING_SIGMAS ...]]
                                      [--n4-shrink-factor N4_SHRINK_FACTOR]
                                      [--n4-iterations N4_ITERATIONS [N4_ITERATIONS ...]]
                                      [--n4-convergence-tolerance N4_CONVERGENCE_TOLERANCE]
                                      [--n4-spline-mesh N4_SPLINE_MESH N4_SPLINE_MESH N4_SPLINE_MESH]

options:
  -h, --help            show this help message and exit
  --t1 T1               input T1 .nii.gz
  --flair FLAIR         moving FLAIR .nii.gz
  --brain-mask BRAIN_MASK
                        same-grid T1 brain mask .nii.gz
  --output-dir OUTPUT_DIR
  --resume
  --metric {mattes,meansquares,GC}
  --metric-bins METRIC_BINS
  --sampling-rate SAMPLING_RATE
  --registration-iterations REGISTRATION_ITERATIONS [REGISTRATION_ITERATIONS ...]
  --shrink-factors SHRINK_FACTORS [SHRINK_FACTORS ...]
  --smoothing-sigmas SMOOTHING_SIGMAS [SMOOTHING_SIGMAS ...]
  --n4-shrink-factor N4_SHRINK_FACTOR
  --n4-iterations N4_ITERATIONS [N4_ITERATIONS ...]
  --n4-convergence-tolerance N4_CONVERGENCE_TOLERANCE
  --n4-spline-mesh N4_SPLINE_MESH N4_SPLINE_MESH N4_SPLINE_MESH
```

Defaults: `--resume` = `False`; `--metric` = `mattes`; `--metric-bins` = `32`; `--sampling-rate` = `0.2`; `--registration-iterations` = `[1000, 500, 250, 0]`; `--shrink-factors` = `[6, 4, 2, 1]`; `--smoothing-sigmas` = `[3.0, 2.0, 1.0, 0.0]`; `--n4-shrink-factor` = `4`; `--n4-iterations` = `[50, 50, 50, 50]`; `--n4-convergence-tolerance` = `0.0`; `--n4-spline-mesh` = `[1, 1, 1]`.

## `plus-crop-lesion-inputs`

```text
usage: pymacruise plus-crop-lesion-inputs [-h] --modalities MODALITIES
                                          --topology TOPOLOGY
                                          --output-dir OUTPUT_DIR
                                          [--padding PADDING] [--resume]

options:
  -h, --help            show this help message and exit
  --modalities MODALITIES
                        4-D T1-then-FLAIR .nii.gz
  --topology TOPOLOGY   LesionTOADS topology atlas .nii.gz
  --output-dir OUTPUT_DIR
  --padding PADDING
  --resume
```

Defaults: `--padding` = `10`; `--resume` = `False`.

## `plus-lesion-initial-align`

```text
usage: pymacruise plus-lesion-initial-align [-h]
                                            --cropped-modalities CROPPED_MODALITIES
                                            --crop-manifest CROP_MANIFEST
                                            --atlas-descriptor ATLAS_DESCRIPTOR
                                            --output-dir OUTPUT_DIR
                                            [--init-align-iterations INIT_ALIGN_ITERATIONS]
                                            [--pyramid-levels PYRAMID_LEVELS]
                                            [--resume]

options:
  -h, --help            show this help message and exit
  --cropped-modalities CROPPED_MODALITIES
  --crop-manifest CROP_MANIFEST
  --atlas-descriptor ATLAS_DESCRIPTOR
  --output-dir OUTPUT_DIR
  --init-align-iterations INIT_ALIGN_ITERATIONS
  --pyramid-levels PYRAMID_LEVELS
  --resume
```

Defaults: `--init-align-iterations` = `10`; `--pyramid-levels` = `2`; `--resume` = `False`.

## `plus-lesion-registered-recrop`

```text
usage: pymacruise plus-lesion-registered-recrop [-h]
                                                --normalized-modalities NORMALIZED_MODALITIES
                                                --crop-manifest CROP_MANIFEST
                                                --initial-alignment-manifest INITIAL_ALIGNMENT_MANIFEST
                                                --atlas-descriptor ATLAS_DESCRIPTOR
                                                --output-dir OUTPUT_DIR
                                                [--resume]

options:
  -h, --help            show this help message and exit
  --normalized-modalities NORMALIZED_MODALITIES
  --crop-manifest CROP_MANIFEST
  --initial-alignment-manifest INITIAL_ALIGNMENT_MANIFEST
  --atlas-descriptor ATLAS_DESCRIPTOR
  --output-dir OUTPUT_DIR
  --resume
```

Defaults: `--resume` = `False`.

## `plus-lesion-first-iteration`

```text
usage: pymacruise plus-lesion-first-iteration [-h]
                                              --registered-recrop-manifest REGISTERED_RECROP_MANIFEST
                                              --atlas-descriptor ATLAS_DESCRIPTOR
                                              --output-dir OUTPUT_DIR
                                              [--iterations ITERATIONS]
                                              [--first-thinning-limit FIRST_THINNING_LIMIT]
                                              [--first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO]
                                              [--finite-zero-variance]
                                              [--resume]

options:
  -h, --help            show this help message and exit
  --registered-recrop-manifest REGISTERED_RECROP_MANIFEST
  --atlas-descriptor ATLAS_DESCRIPTOR
  --output-dir OUTPUT_DIR
  --iterations ITERATIONS
  --first-thinning-limit FIRST_THINNING_LIMIT
  --first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO
  --finite-zero-variance
  --resume
```

Defaults: `--iterations` = `99`; `--first-thinning-limit` = `3.0`; `--first-thinning-volume-ratio` = `0.5`; `--finite-zero-variance` = `False`; `--resume` = `False`.

## `plus-lesion-registered-branch`

```text
usage: pymacruise plus-lesion-registered-branch [-h]
                                                --registered-recrop-manifest REGISTERED_RECROP_MANIFEST
                                                --atlas-descriptor ATLAS_DESCRIPTOR
                                                --output-dir OUTPUT_DIR
                                                [--iterations ITERATIONS]
                                                [--first-thinning-limit FIRST_THINNING_LIMIT]
                                                [--first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO]
                                                [--max-distance MAX_DISTANCE]
                                                [--first-iteration-manifest FIRST_ITERATION_MANIFEST]
                                                [--finite-zero-variance]
                                                [--resume]

options:
  -h, --help            show this help message and exit
  --registered-recrop-manifest REGISTERED_RECROP_MANIFEST
  --atlas-descriptor ATLAS_DESCRIPTOR
  --output-dir OUTPUT_DIR
  --iterations ITERATIONS
  --first-thinning-limit FIRST_THINNING_LIMIT
  --first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO
  --max-distance MAX_DISTANCE
  --first-iteration-manifest FIRST_ITERATION_MANIFEST
  --finite-zero-variance
  --resume
```

Defaults: `--iterations` = `99`; `--first-thinning-limit` = `3.0`; `--first-thinning-volume-ratio` = `0.5`; `--max-distance` = `0.0001`; `--first-iteration-manifest` = `None`; `--finite-zero-variance` = `False`; `--resume` = `False`.

## `plus-lesion-registered-loop`

```text
usage: pymacruise plus-lesion-registered-loop [-h]
                                              --registered-checkpoint-manifest REGISTERED_CHECKPOINT_MANIFEST
                                              --atlas-descriptor ATLAS_DESCRIPTOR
                                              --output-dir OUTPUT_DIR
                                              [--iterations ITERATIONS]
                                              [--thinning-limit THINNING_LIMIT]
                                              [--max-passes MAX_PASSES]
                                              [--correction-fields CORRECTION_FIELDS]
                                              [--correction-polynomial-degree {1,2,3}]
                                              [--initial-correction-degree INITIAL_CORRECTION_DEGREE]
                                              [--finite-zero-variance]
                                              [--resume]

options:
  -h, --help            show this help message and exit
  --registered-checkpoint-manifest REGISTERED_CHECKPOINT_MANIFEST
  --atlas-descriptor ATLAS_DESCRIPTOR
  --output-dir OUTPUT_DIR
  --iterations ITERATIONS
  --thinning-limit THINNING_LIMIT
  --max-passes MAX_PASSES
  --correction-fields CORRECTION_FIELDS
                        optional channel-last RAS .nii.gz correction fields
                        matching the registered modalities
  --correction-polynomial-degree {1,2,3}
                        refit and carry adaptive correction fields with this
                        source degree ceiling
  --initial-correction-degree INITIAL_CORRECTION_DEGREE
  --finite-zero-variance
  --resume
```

Defaults: `--iterations` = `99`; `--thinning-limit` = `3.0`; `--max-passes` = `1`; `--correction-fields` = `None`; `--correction-polynomial-degree` = `None`; `--initial-correction-degree` = `None`; `--finite-zero-variance` = `False`; `--resume` = `False`.

## `plus-lesion-loop-run`

```text
usage: pymacruise plus-lesion-loop-run [-h]
                                       --registered-recrop-manifest REGISTERED_RECROP_MANIFEST
                                       --atlas-descriptor ATLAS_DESCRIPTOR
                                       --output-dir OUTPUT_DIR
                                       [--iterations ITERATIONS]
                                       [--first-thinning-limit FIRST_THINNING_LIMIT]
                                       [--first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO]
                                       [--max-distance MAX_DISTANCE]
                                       [--loop-convergence-threshold LOOP_CONVERGENCE_THRESHOLD]
                                       [--n-iter-max N_ITER_MAX]
                                       [--thinning-limit THINNING_LIMIT]
                                       [--max-passes MAX_PASSES]
                                       [--correction-fields CORRECTION_FIELDS]
                                       [--correction-polynomial-degree {1,2,3}]
                                       [--initial-correction-degree INITIAL_CORRECTION_DEGREE]
                                       [--finite-zero-variance] [--resume]

options:
  -h, --help            show this help message and exit
  --registered-recrop-manifest REGISTERED_RECROP_MANIFEST
  --atlas-descriptor ATLAS_DESCRIPTOR
  --output-dir OUTPUT_DIR
  --iterations ITERATIONS
  --first-thinning-limit FIRST_THINNING_LIMIT
  --first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO
  --max-distance MAX_DISTANCE
  --loop-convergence-threshold LOOP_CONVERGENCE_THRESHOLD
                        optional P41m convergence threshold, separate from
                        P41k max-distance
  --n-iter-max N_ITER_MAX
  --thinning-limit THINNING_LIMIT
  --max-passes MAX_PASSES
  --correction-fields CORRECTION_FIELDS
  --correction-polynomial-degree {1,2,3}
  --initial-correction-degree INITIAL_CORRECTION_DEGREE
  --finite-zero-variance
                        use finite limiting similarities; default preserves
                        source-faithful NaN behavior
  --resume
```

Defaults: `--iterations` = `99`; `--first-thinning-limit` = `3.0`; `--first-thinning-volume-ratio` = `0.5`; `--max-distance` = `0.0001`; `--loop-convergence-threshold` = `None`; `--n-iter-max` = `99`; `--thinning-limit` = `3.0`; `--max-passes` = `99`; `--correction-fields` = `None`; `--correction-polynomial-degree` = `None`; `--initial-correction-degree` = `None`; `--finite-zero-variance` = `False`; `--resume` = `False`.

## `plus-full-run`

```text
usage: pymacruise plus-full-run [-h] --t1 T1 --flair FLAIR
                                --brain-mask BRAIN_MASK --maseg MASEG
                                --atlas-descriptor ATLAS_DESCRIPTOR
                                [--atlasem-descriptor ATLASEM_DESCRIPTOR]
                                --output-dir OUTPUT_DIR
                                [--regrid-inputs | --no-regrid-inputs]
                                [--dura-transfer-mode {physical-affine,jist-unchanged}]
                                [--registration-iterations REGISTRATION_ITERATIONS [REGISTRATION_ITERATIONS ...]]
                                [--shrink-factors SHRINK_FACTORS [SHRINK_FACTORS ...]]
                                [--smoothing-sigmas SMOOTHING_SIGMAS [SMOOTHING_SIGMAS ...]]
                                [--registration-sampling-rate REGISTRATION_SAMPLING_RATE]
                                [--n4-shrink-factor N4_SHRINK_FACTOR]
                                [--n4-iterations N4_ITERATIONS [N4_ITERATIONS ...]]
                                [--n4-convergence-tolerance N4_CONVERGENCE_TOLERANCE]
                                [--n4-spline-mesh N4_SPLINE_MESH N4_SPLINE_MESH N4_SPLINE_MESH]
                                [--crop-padding CROP_PADDING]
                                [--init-align-iterations INIT_ALIGN_ITERATIONS]
                                [--pyramid-levels PYRAMID_LEVELS]
                                [--lesion-iterations LESION_ITERATIONS]
                                [--first-thinning-limit FIRST_THINNING_LIMIT]
                                [--first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO]
                                [--max-distance MAX_DISTANCE]
                                [--loop-convergence-threshold LOOP_CONVERGENCE_THRESHOLD]
                                [--n-iter-max N_ITER_MAX]
                                [--loop-thinning-limit LOOP_THINNING_LIMIT]
                                [--max-passes MAX_PASSES]
                                [--finite-zero-variance] [--resume]
                                [--stop-after-stage {preprocessing,atlasem_dura,crop,initial_alignment,registered_recrop,lesion_workflow,dura_outputs,native_dura_transfer,maseg_cortical_separation,plus_surface_workflow,maseg_refinement}]

options:
  -h, --help            show this help message and exit
  --t1 T1
  --flair FLAIR
  --brain-mask BRAIN_MASK
  --maseg MASEG
  --atlas-descriptor ATLAS_DESCRIPTOR
  --atlasem-descriptor ATLASEM_DESCRIPTOR
                        optional AtlasEM descriptor (testing target: adapted
                        ibsr11objV3) for the source-order pre-LesionTOADS
                        AtlasEM/Dura handoff; when omitted, run the legacy
                        LesionTOADS-only chain
  --output-dir OUTPUT_DIR
  --regrid-inputs, --no-regrid-inputs
                        recommended: resample oblique RAS T1, mask, and MaSeg
                        to an orthogonal grid (default: enabled; --no-regrid-
                        inputs only permutes/flips axes to closest RAS,
                        discards oblique directions, and does not interpolate)
  --dura-transfer-mode {physical-affine,jist-unchanged}
                        Dura mapping to native T1: recommended physical-affine
                        default, or jist-unchanged for legacy source-grid
                        compatibility
  --registration-iterations REGISTRATION_ITERATIONS [REGISTRATION_ITERATIONS ...]
  --shrink-factors SHRINK_FACTORS [SHRINK_FACTORS ...]
  --smoothing-sigmas SMOOTHING_SIGMAS [SMOOTHING_SIGMAS ...]
  --registration-sampling-rate REGISTRATION_SAMPLING_RATE
  --n4-shrink-factor N4_SHRINK_FACTOR
  --n4-iterations N4_ITERATIONS [N4_ITERATIONS ...]
  --n4-convergence-tolerance N4_CONVERGENCE_TOLERANCE
  --n4-spline-mesh N4_SPLINE_MESH N4_SPLINE_MESH N4_SPLINE_MESH
  --crop-padding CROP_PADDING
  --init-align-iterations INIT_ALIGN_ITERATIONS
  --pyramid-levels PYRAMID_LEVELS
  --lesion-iterations LESION_ITERATIONS
  --first-thinning-limit FIRST_THINNING_LIMIT
  --first-thinning-volume-ratio FIRST_THINNING_VOLUME_RATIO
  --max-distance MAX_DISTANCE
  --loop-convergence-threshold LOOP_CONVERGENCE_THRESHOLD
  --n-iter-max N_ITER_MAX
  --loop-thinning-limit LOOP_THINNING_LIMIT
  --max-passes MAX_PASSES
  --finite-zero-variance
  --resume
  --stop-after-stage {preprocessing,atlasem_dura,crop,initial_alignment,registered_recrop,lesion_workflow,dura_outputs,native_dura_transfer,maseg_cortical_separation,plus_surface_workflow,maseg_refinement}
                        checkpoint after this stage; repeat with --resume to
                        continue
```

Defaults: `--atlasem-descriptor` = `None`; `--regrid-inputs` = `True`; `--dura-transfer-mode` = `physical-affine`; `--registration-iterations` = `[1000, 500, 250, 0]`; `--shrink-factors` = `[6, 4, 2, 1]`; `--smoothing-sigmas` = `[3.0, 2.0, 1.0, 0.0]`; `--registration-sampling-rate` = `0.25`; `--n4-shrink-factor` = `4`; `--n4-iterations` = `[50, 50, 50, 50]`; `--n4-convergence-tolerance` = `0.0`; `--n4-spline-mesh` = `[1, 1, 1]`; `--crop-padding` = `10`; `--init-align-iterations` = `10`; `--pyramid-levels` = `2`; `--lesion-iterations` = `99`; `--first-thinning-limit` = `3.0`; `--first-thinning-volume-ratio` = `0.5`; `--max-distance` = `0.0001`; `--loop-convergence-threshold` = `None`; `--n-iter-max` = `99`; `--loop-thinning-limit` = `3.0`; `--max-passes` = `99`; `--finite-zero-variance` = `False`; `--resume` = `False`; `--stop-after-stage` = `None`.

## `plus-run`

```text
usage: pymacruise plus-run [-h] --t1-reference T1_REFERENCE
                           --segmentation SEGMENTATION --dura-csf DURA_CSF
                           --dura-gm DURA_GM --dura-wm-fill DURA_WM_FILL
                           --dura-wm-mask DURA_WM_MASK
                           --multi-atlas-gm MULTI_ATLAS_GM
                           --multi-atlas-wm MULTI_ATLAS_WM
                           --multi-atlas-mask MULTI_ATLAS_MASK
                           --brain-mask BRAIN_MASK --output-dir OUTPUT_DIR
                           [--dura-transfer-mode {physical-affine,jist-unchanged}]
                           [--resume] [--crop-edge-width CROP_EDGE_WIDTH]
                           [--gvf-lambda GVF_LAMBDA]
                           [--wm-iso-level WM_ISO_LEVEL]
                           [--refine-topology-by-maseg | --no-refine-topology-by-maseg]

options:
  -h, --help            show this help message and exit
  --t1-reference T1_REFERENCE
  --segmentation SEGMENTATION
  --dura-csf DURA_CSF
  --dura-gm DURA_GM
  --dura-wm-fill DURA_WM_FILL
  --dura-wm-mask DURA_WM_MASK
  --multi-atlas-gm MULTI_ATLAS_GM
  --multi-atlas-wm MULTI_ATLAS_WM
  --multi-atlas-mask MULTI_ATLAS_MASK
                        MaSeg Atlas_BrainTissueMask (.nii.gz; GM+WM+other
                        brain tissue), not the WM-only Atlas_Mask
  --brain-mask BRAIN_MASK
  --output-dir OUTPUT_DIR
  --dura-transfer-mode {physical-affine,jist-unchanged}
                        Dura mapping to native T1: recommended physical-affine
                        default, or jist-unchanged for legacy source-grid
                        compatibility
  --resume
  --crop-edge-width CROP_EDGE_WIDTH
  --gvf-lambda GVF_LAMBDA
  --wm-iso-level WM_ISO_LEVEL
  --refine-topology-by-maseg, --no-refine-topology-by-maseg
```

Defaults: `--dura-transfer-mode` = `physical-affine`; `--resume` = `False`; `--crop-edge-width` = `5`; `--gvf-lambda` = `0.2`; `--wm-iso-level` = `0.5`; `--refine-topology-by-maseg` = `True`.

## `pymacruise-lesion-output`

Also available as `python -m pymacruise.lesion_output_cli`. This writes already-computed arrays; it does not run segmentation.

```text
usage: python -m pymacruise.lesion_output_cli [-h] --reference REFERENCE
                                              --classification CLASSIFICATION
                                              --lesion-subtype LESION_SUBTYPE
                                              --mode {classification-only,class+member,cruise inputs,dura removal inputs}
                                              [--membership-classification MEMBERSHIP_CLASSIFICATION]
                                              [--mode-data MODE_DATA]
                                              [--inhomogeneity-field INHOMOGENEITY_FIELD]
                                              --input-name INPUT_NAME
                                              --output-dir OUTPUT_DIR

Write already-computed LesionTOADS terminal arrays as typed RAS NIfTI outputs.
This is an output adapter, not a segmentation run.

options:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        original T1 reference .nii.gz
  --classification CLASSIFICATION
                        terminal class labels .nii.gz
  --lesion-subtype LESION_SUBTYPE
                        lesion subtype labels .nii.gz
  --mode {classification-only,class+member,cruise inputs,dura removal inputs}
                        selected terminal output branch
  --membership-classification MEMBERSHIP_CLASSIFICATION
                        optional membership-derived class labels .nii.gz
  --mode-data MODE_DATA
                        mode product .nii.gz: [x,y,z,channels] for
                        class+member, or [4,x,y,z] in CSF/GM/WM-fill/WM-mask
                        order for CRUISE and Dura
  --inhomogeneity-field INHOMOGENEITY_FIELD
                        optional already-restored 3-D or 4-D field .nii.gz
  --input-name INPUT_NAME
                        input T1 filename used for output stem
  --output-dir OUTPUT_DIR
```
