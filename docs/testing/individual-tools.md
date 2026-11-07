# Running the tools separately

Updated October 27, 2026. Install and prepare a testing kit using [the full PLUS walkthrough](full-plus-walkthrough.md). The same inputs, settings, image grids and atlas identities apply to individual runs. Scientific intermediates can be reused within one pipeline for troubleshooting; copying them across pipelines changes the independence boundary and must be reported as a held-fixed experiment.

## Inspect one Python workflow stage at a time

The generated script uses one explicit full scientific profile. Stop it at any of these stage boundaries. Earlier verified stages are resumed and the next prefix is computed, so dependent filenames and checkpoint signatures remain correct:

```bash
export PLUS_KIT="$PWD/build/plus-e2e"
bash "$PLUS_KIT/run-python.sh" --stop-after-stage preprocessing \
  > "$PLUS_KIT/01-preprocessing.log" 2>&1
bash "$PLUS_KIT/run-python.sh" --resume --stop-after-stage atlasem_dura \
  > "$PLUS_KIT/02-atlasem-dura.log" 2>&1
# Continue with any next boundary from the table below, preserving --resume.
bash "$PLUS_KIT/run-python.sh" --resume --stop-after-stage crop \
  > "$PLUS_KIT/03-crop.log" 2>&1
# Finish the remaining stages:
bash "$PLUS_KIT/run-python.sh" --resume > "$PLUS_KIT/finish.log" 2>&1
```

| Boundary | Scientific work and outputs to inspect |
| --- | --- |
| `preprocessing` | Rigid FLAIR-to-T1 registration, transform files, N4 on both images, user masking, negative clamp, T1/FLAIR modality stack. |
| `atlasem_dura` | Selected AtlasEM initialization/registration/EM, GM/WM-fill/WM-mask, initial Dura evolution, stripped T1 and masked registered FLAIR. |
| `crop` | LesionTOADS normalization and first crop; normalized full grid, cropped modalities, topology and crop metadata. |
| `initial_alignment` | Descriptor-ordered memberships and initial rigid atlas alignment. |
| `registered_recrop` | Post-registration uncrop/recrop, updated topology and transformed shape maps. |
| `lesion_workflow` | First iteration and continuing topology/thinning/statistics/membership/registration loop; terminal classification and lesion subtype when available. |
| `dura_outputs` | LesionTOADS-derived CSF, GM, WM-fill and WM-mask outputs. This is output routing, not a second Dura evolution. |
| `native_dura_transfer` | Transfer those channels back to native T1 physical coordinates. |
| `maseg_cortical_separation` | MaSeg-derived GM, WM, mask and brain-tissue mask. |
| `plus_surface_workflow` | Crop/boundary cleanup, tissue fusion, topology, ACE, GVF, TGDM, uncrop, surfaces and thickness. |
| `maseg_refinement` | Source-rule final full-grid `SEG_refine.nii.gz`; reaching this final boundary produces a complete workflow. |

The parent `python-run/plus-full-workflow.json` points to each active child manifest and its hashed products. Every intermediate prefix has `status: checkpointed`, `terminal: false`; it is deliberately incomplete. A stage failure has `status: failed` and an error. Use a new output kit for a settings change; use identical settings plus `--resume` for an interruption. Do not edit manifests. After an interrupted longer attempt, a retained later child may still record `running` even though its process has exited. A verified prefix remains `checkpointed` and incomplete; the next full resume verifies/recomputes the interrupted child. Use the current process/log record rather than treating a saved `running` string as proof that a process is alive.

## Run a Python scientific tool directly

The [complete command reference](../tools-reference.md) contains exact usage, all flags and defaults for **every** current CLI tool. The [detailed I/O reference](../cli-help.md) explains output names, geometry and algorithm limits. Each invocation is independent and can use your own correctly prepared files; the commands below illustrate the main scientific tools. Names such as `GM.nii.gz` denote actual same-grid intermediates, not bundled example filenames.

```bash
# One registration checkpoint; FLAIR may have a different physical grid.
pymacruise plus-register-flair --t1 T1.nii.gz --flair FLAIR.nii.gz --output-dir registration

# Paired registration, N4 and mask application, before AtlasEM/LesionTOADS.
pymacruise plus-prepare-lesion --t1 T1.nii.gz --flair FLAIR.nii.gz \
  --brain-mask brain_mask.nii.gz --output-dir prepared

# One N4 correction, with source-style defaults (no mask/intensity rescale).
pymacruise n4 --input T1.nii.gz --output T1-n4.nii.gz --provenance n4.json

# MaSeg cortical masks, before fusion with LesionTOADS outputs.
pymacruise cortical-separation --maseg MaSeg.nii.gz --output-dir cortical-masks

# AtlasEM from an already corrected/masked T1, and a prepared NIfTI descriptor.
pymacruise atlasem-run --t1 T1-n4-masked.nii.gz --brain-mask brain_mask.nii.gz \
  --descriptor "$PLUS_KIT/python-atlas/LesionToads-11-obj-atlasV3-atlasem-python.txt" \
  --output-dir atlasem --modality T1_MPRAGE --levels 4 --initial-registration-iterations 50

# Standalone topology correction; use the selected full stage for PLUS-specific fusion.
pymacruise topology-correct --volume WM.nii.gz --paint-mask initial-membership.nii.gz \
  --output initial-corrected.nii.gz --provenance topology.json

# ACE and GVF: use aligned membership maps, not hard labels.
pymacruise ace --gray GM.nii.gz --white WM.nii.gz \
  --enhanced GM-ACE.nii.gz --skeleton ACE.nii.gz \
  --thinned-skeleton ACE-thinned.nii.gz --provenance ace.json
pymacruise gvf --input GM-ACE.nii.gz --output GVF.nii.gz --lambda 0.2 --provenance gvf.json

# Nested TGDM: INITIAL is an initial membership, not a signed level set.
# The selected PLUS wrapper supplies no optional mask bound.
pymacruise tgdm --initial initial-membership.nii.gz --white WM.nii.gz \
  --gray GM.nii.gz --ace-gray GM-ACE.nii.gz --ace-skeleton ACE-thinned.nii.gz \
  --gvf GVF.nii.gz --output-dir tgdm --provenance tgdm.json

# Dura evolution from T1 and AtlasEM (or selected TOADS) handoff maps.
pymacruise dura --t1 T1-n4-masked.nii.gz --white WM-fill.nii.gz \
  --gray GM.nii.gz --white-mask WM-mask.nii.gz --output-dir dura --provenance dura.json

# Meshes and thickness, on an agreed full physical grid.
pymacruise surface-extract --level-set inner.nii.gz --output inner.vtk --provenance inner.json
pymacruise surface-extract --level-set outer.nii.gz --output outer.vtk --provenance outer.json
pymacruise thickness-levelset --inner inner.nii.gz --outer outer.nii.gz \
  --output thickness.nii.gz --provenance thickness.json
pymacruise thickness-surface --inner inner.vtk --outer outer.vtk \
  --inner-output inner-measured.vtk --provenance surface-thickness.json
pymacruise segmentation-refine --segmentation MaSeg.nii.gz \
  --inner-levelset inner.nii.gz --outer-levelset outer.nii.gz \
  --output SEG_refine.nii.gz --summary refinement-summary.json --source-faithful

# Verify a complete Python run, including active manifests and output hashes.
python scripts/audit_python_plus_run.py "$PLUS_KIT/python-run" --report "$PLUS_KIT/python-audit.json"
```

The toolkit's `plus_surface_workflow` has PLUS-specific boundary/fusion and topology options. Running the plain standalone `ace`, `tgdm` or topology command with defaults is a tool experiment, not automatically a replay of that full stage. Consult resolved settings and the public `PlusPipelineOptions` before comparing a hand-assembled chain. The direct TGDM command interprets its initial image as membership; use the selected full-stage replay API for an existing signed initial level set.

The independent restartable commands for the LesionTOADS chain are:

| Tool | Inputs / purpose |
| --- | --- |
| `plus-crop-lesion-inputs` | Prepared T1/FLAIR modality stack and atlas topology; writes normalization/crop checkpoint. |
| `plus-lesion-initial-align` | Cropped modalities, crop manifest and descriptor; writes initial alignment. |
| `plus-lesion-registered-recrop` | Full normalized modalities, crop/alignment manifests and descriptor. |
| `plus-lesion-first-iteration` | Registered recrop manifest and descriptor; first memberships/thinning. |
| `plus-lesion-registered-branch` | First-iteration and recrop checkpoints; one continuing branch. |
| `plus-lesion-registered-loop` | Branch checkpoint and descriptor; bounded iteration segment. |
| `plus-lesion-loop-run` | Registered recrop and descriptor; runs initialization through continuing loop and terminal exports. |
| `atlasem-run` | Corrected/masked T1, prepared AtlasEM descriptor; exports AtlasEM products and manifest. |
| `atlasem-prelesion-dura` | T1, FLAIR, brain mask and existing AtlasEM products; Dura handoff. |
| `plus-run` | Persisted Dura channels plus aligned MaSeg-derived maps; downstream surface/thickness continuation. |

Run `pymacruise TOOL --help` or use the generated reference for the exact required option names and example output filenames. There is no single standalone command that reproduces every historical LesionTOADS plugin mode. The selected full LesionTOADS route is exposed through the checkpoint commands above and `plus-full-run`.

## Export computed LesionTOADS products and run tester utilities

The separate output adapter materializes already-computed arrays; it does not segment an image:

```bash
pymacruise-lesion-output --reference T1.nii.gz \
  --classification classification.nii.gz --lesion-subtype lesion-subtype.nii.gz \
  --mode classification-only --input-name T1.nii.gz --output-dir exported-lesion-products
```

`class+member`, `cruise inputs` and `dura removal inputs` additionally require `--mode-data` with their documented array layout. See the [exporter reference](../tools-reference.md#pymacruise-lesion-output). Keep the original reference geometry and terminal arrays; this adapter is not a way to turn arbitrary tissue maps into a complete LesionTOADS run.

| Tester utility | Invocation and purpose |
| --- | --- |
| Bundled asset inventory | `python -m pymacruise.asset_inventory` verifies packaged files, headers and hashes. |
| Kit preparation | `python scripts/prepare_plus_testing.py --input-dir INPUTS --output-dir NEW_KIT` prepares native/Python atlas resources and both run configurations. |
| Native-atlas converter | `python scripts/adapt_ibsr11objv3_atlas.py --descriptor NATIVE.txt --output-dir NEW_ATLAS` converts a separately unpacked native atlas; the kit normally calls this for you. |
| Python integrity audit | `python scripts/audit_python_plus_run.py RUN --report audit.json` requires a complete terminal run and verifies its active outputs. |
| Original output audit | `python scripts/audit_source_reference_run.py SOURCE_RUN --require-run-info --require-ras` inventories returned original files and NIfTI geometry. |
| Original run recording | `python scripts/record_plus_reference.py --help` lists the required runtime, input-kit, mesh-frame and source-output arguments; the full walkthrough gives the complete invocation. |
| PLUS comparison | `python scripts/compare_external_plus_run.py --help` lists the two run trees, four input paths, reference-info/log/map and report flags; use the full walkthrough to preserve the comparison boundary. |
| BASE comparison | `python scripts/compare_external_base_run.py --help` is separate from PLUS; use [the BASE walkthrough](full-base-walkthrough.md). |

Developer `generate_java_*`, replay and comparison harnesses under `scripts/` test specific held-fixed algorithms. They are supplemental validation tools and do not replace either independent scientific workflow.

## Run original tools separately

Use the same original MATLAB source and installed JIST classes as the full source driver. The [original driver](../../source-reference-run/matlab/run_plus_selected_mask_reference.m) is an executable record of input selection, stage order, output variables and per-stage options; copy the relevant calls into a separate MATLAB script for a held-fixed diagnostic. The following setup is sufficient for resolving the helpers; replace paths with the actual runtime paths:

```matlab
addpath('/path/to/separately-obtained/reference-drivers');
addpath(genpath('/path/to/MaCRUISE/pipeline'));
addpath(genpath('/path/to/full-multi-atlas/masi-fusion/src'));
opts = struct();
opts.mipav = '/path/to/MIPAV-install';
opts.plugindir = fullfile(opts.mipav, 'plugins');
opts.multiatlas_code_dir = '/path/to/full-multi-atlas';
opts.N4correctBIN_dir = '/path/to/ANTs/bin/N4BiasFieldCorrection';
opts.jvmmemval = '8G'; opts.runmethod = 'single';
opts.working_dir = '/path/to/new-tool-workdir';
mkdir(opts.working_dir);
opts.T1_bname = get_basename('/path/to/T1.nii.gz');
```

Use a new directory for each independent tool attempt. Legacy source helpers often skip existing products without verifying hashes. Keep the exact MATLAB calls and resulting Java command/log files.

| Original tool | Direct call and required handoff |
| --- | --- |
| Rigid FLAIR registration | `opts.rigOnly=1; opts.interp=[]; registered=FlairRegistration(T1,FLAIR,regdir,opts);` Create `regdir` first; Aladin binary is under `full-multi-atlas/niftyreg/bin/`. |
| N4 | `opts.N3orN4='N4'; corrected=run_N4Correct(image,opts,'N4','N4');` Set `opts.T1_bname` separately for T1 and registered FLAIR. |
| Mask / negative cleanup | Run the provided `apply_brain_mask_nifti.py --image IMAGE --mask MASK --output OUTPUT`, then original `removeNegative(masked,opts,'clamped','T1')`. |
| AtlasEM | `opts.atlasem_descriptor='/absolute/kit/native-atlas/ibsr11objV3-atlasem-java.txt'; [gm,wmfill,wmmask]=run_plus_atlasem_reference(T1masked,opts,'EM','AtlasEM');` Runs the original AtlasEM command with explicit selected atlas. |
| Initial Dura | `[stripped,mask]=run_RemoveDura(T1masked,gm,wmfill,wmmask,opts,'Dura','Dura');` Apply returned mask to prepared FLAIR using the supplied helper. |
| LesionTOADS | `opts.LesionToadsOutput='dura removal inputs'; opts.LesionToadsCorrect='false'; [gm,wmfill,wmmask]=run_LesionToads(T1dura,FLAIRdura,opts,'LesionTOADS','LesionTOADS');` Original wrapper supplies its atlas and iteration settings. |
| MaSeg separation | `[gm,wm,mask,brain]=cortical_seperation(MaSeg,opts.working_dir,'CorticalSeperation','Separation');` |
| Transfer and crop | `run_TransferToOriginSpace`, then `run_CropBaseOnGM`; use the driver call signatures and preserve the returned cell-array ordering. |
| Boundary and tissue fusion | `run_GetMABoundary`, `run_CleanBoundaryByToads`, `update_WMandGM`, `run_CleanBoundaryByUpdateWM`; use the source driver's overlook-label, erosion and working-directory settings. |
| Topology correction | `corrected=run_TopoCorrect(WM,initialMask,opts,'TopoCorrect','Topo');` Full PLUS additionally sets MaSeg refinement options; see source surface helper. |
| ACE | `[csfRaw,csfDiffused,csfGaussian,enhanced,skeleton,thin]=run_Ace(GM,WM,opts,'AceOrigin','ACE');` Set `NewCSFfile`, `UseUpdatedCSFinACE=3` and the selected surface options as in the driver. |
| GVF | `opts.GVF_Lambda=0.2; field=run_GVF(enhanced,opts,'GVF','GVF');` |
| TGDM / uncrop / meshes | `[inner,central,outer]=run_TGDMinvcrop(fullGridGM,initialMask,WM,GM,enhanced,thin,field,opts,'TGDMinvcrop','TGDM');` Uses original selected null-bound wrapper. |
| Refinement | `refined=run_SegmentationRefine(plusRoot,MaSeg,refineDir,opts);` Requires the original level-set folder layout. |
| Volume and surface thickness | `[volume,mat]=run_plus_thickness_reference(inner,outer,innerVTK,outerVTK,opts);` Runs original thickness helpers and retains nested logs/measurement mesh. |

For explicit original transfer/crop/fusion calls, first assign `AtlasGM`, `AtlasWM`, `AtlasMask`, `AtlasBrainMask` to cortical-separation outputs and `LesionGM`, `LesionWMfill`, `LesionWMmask` to the original LesionTOADS handoff. With `MaSeg` set to the same original label file:

```matlab
transferred = run_TransferToOriginSpace(MaSeg, ...
    {LesionGM, LesionWMfill, LesionWMmask}, opts, 'Transfer', 'Transfer');
[GMraw, cropped] = run_CropBaseOnGM(AtlasGM, ...
    [{AtlasWM, AtlasMask, AtlasBrainMask, MaSeg}, transferred], opts, 'Crop', 'Crop');
WMraw = cropped{1}; MaskRaw = cropped{2}; Brain = cropped{3}; Seg = cropped{4};
GMtoads = cropped{5}; WMtoads = cropped{6}; WMmaskToads = cropped{7};
opts.ErodeMethod = '3D'; opts.OverLookLabels = '120,121,202,203,156,157';
[Boundary, ~] = run_GetMABoundary(Seg, opts, 'Boundary', 'Boundary');
Boundary = run_CleanBoundaryByToads(Boundary, WMmaskToads, opts, 'Boundary', 'Clean');
[GM, WM, Mask] = update_WMandGM(GMraw, WMraw, MaskRaw, Brain, ...
    GMtoads, WMtoads, WMmaskToads, opts, 'Fusion', 'Fusion');
Boundary = run_CleanBoundaryByUpdateWM(Boundary, WM, opts, 'Boundary', 'CleanWM');
opts.update_Enhance = 0; opts.WMIsoLevel = 0.5; opts.GVF_Lambda = 0.2;
opts.IfGenerateIntermediate = 'true'; opts.NewCSFfile = Boundary;
opts.UseUpdatedCSFinACE = 3; opts.MaWM = WMraw; opts.refineTopoCorrectByMa = 1;
% Now use these exact GM/WM/Mask fields and opts in the ACE/GVF/TGDM
% calls above, or invoke the complete original surface-stage wrapper:
run_CruisePlusSurf_pipeline(T1, AtlasGM, GM, WM, Mask, opts, 'CruisePlusSurfSeperate', 'Surface');
```

For the exact JIST CLI behind each tool, its original `run_*_jist.m` constructs and returns the command. Set `opts.output_dir` to a new tool directory, call the constructor with the same inputs/options, save the returned text, and execute it with `system`; reject a nonzero status and require the expected files. `run_EM_jist` still hard-codes the historical atlas, so use `run_plus_atlasem_reference` for the selected AtlasEM target. Do not point it at the Python NIfTI descriptor.

For lower-level Python-only checkpoint APIs without standalone commands (`run_plus_lesion_dura_output_checkpoint`, `run_plus_lesion_native_dura_transfer`, `run_plus_maseg_cortical_checkpoint`, `run_plus_segmentation_refinement_checkpoint`), the stage-boundary route above provides a supported invocation with verified upstream dependencies. Their implementation signatures are in `src/pymacruise/plus_*`; passing hand-written checkpoint JSON is unsupported.
