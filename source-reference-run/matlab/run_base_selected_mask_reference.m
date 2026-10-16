function outputs = run_base_selected_mask_reference(maCRUISE_source_dir, mipav_dir, multiatlas_code_dir, n4_binary, python_executable)
% Run the selected v3 BASE algorithm stages with the supplied T1-grid mask.
%
% This run-only driver intentionally does not call run_maCRUISE_v3.m because
% that entry point invokes BET and emits visualization/MNI products. It calls
% the existing source stage functions in maCRUISEpipeline.m order, supplies
% the provided mask by zeroing N4-corrected T1 before TOADS, and leaves all
% intermediates in a new output folder. It does not modify upstream source.
%
% Example:
%   run_base_selected_mask_reference('/path/to/MaCRUISE', '/path/to/mipav', ...
%       '/path/to/full-multi-atlas', '/path/to/N4BiasFieldCorrection', 'python3')

if nargin < 5 || isempty(python_executable), python_executable = 'python3'; end
if nargin < 4 || isempty(n4_binary)
    error('Provide the ANTs N4BiasFieldCorrection executable path.');
end
if nargin < 3 || isempty(multiatlas_code_dir)
    error('Provide the MASI multi-atlas checkout path (used for fusion source/atlas assets).');
end
if nargin < 2 || isempty(mipav_dir)
    error('Provide the MIPAV install path with selected JIST plugins and atlas.');
end
if nargin < 1 || isempty(maCRUISE_source_dir)
    error('Provide the MaCRUISE MATLAB source checkout path.');
end

driver_dir = fileparts(mfilename('fullpath'));
packet_root = fileparts(driver_dir);
pipeline_dir = fullfile(maCRUISE_source_dir, 'pipeline');
input_dir = fullfile(packet_root, 'inputs', 'base-root');
t1_raw = fullfile(input_dir, 'T1.nii.gz');
seg_raw = fullfile(input_dir, 'MaSeg.nii.gz');
brain_mask = fullfile(input_dir, 'brain_mask.nii.gz');
if ~isfile(t1_raw) || ~isfile(seg_raw) || ~isfile(brain_mask)
    error('BASE packet is missing T1, MaSeg, or brain_mask NIfTI input.');
end

stamp = datestr(now, 'yyyymmdd_HHMMSS');
working_dir = fullfile(packet_root, 'outputs', ['base-selected-mask-v3-' stamp]);
if isfolder(working_dir)
    error('Refusing to reuse an existing reference output directory: %s', working_dir);
end
mkdir(working_dir);
diary(fullfile(working_dir, 'matlab-driver.log'));
diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

masi_fusion_dir = fullfile(multiatlas_code_dir, 'masi-fusion', 'src');
if ~isfolder(pipeline_dir) || ~isfolder(masi_fusion_dir)
    error('Could not find MaCRUISE pipeline or MASI fusion source directories.');
end
addpath(genpath(masi_fusion_dir));
addpath(genpath(pipeline_dir));

% Selected v3 wrapper parameters, excluding BET/SPECTRE and all visualization.
opts.runmethod = 'single';
opts.jvmmemval = '8G';
opts.mipav = mipav_dir;
opts.plugindir = fullfile(mipav_dir, 'plugins');
opts.workingrootfolder = working_dir;
opts.working_dir = working_dir;
opts.multiatlas_code_dir = multiatlas_code_dir;
opts.N4correctBIN_dir = n4_binary;
opts.N3orN4 = 'N4';
opts.topoCorrectVersion = 'V3';
opts.T1_sequence = 'T1_MPRAGE';
opts.FinalToadsMethod = 'cruise inputs';
opts.ErodeMethod = '3D';
opts.OverLookLabels = '120,121,202,203,156,157';

% Keep the original V3 raw preparation and multi-atlas tissue preparation.
[T1_file, Seg_file] = prepare_raw(t1_raw, seg_raw, working_dir, 'Raw', 'STEP0');
opts.rawseg = Seg_file;
[GM_file, WM_file, Mask_file, BrainMask_file] = cortical_seperation( ...
    Seg_file, working_dir, 'CorticalSeperation', 'STEP1');

% V3 N4 is unmasked and uses ANTs defaults. Apply only the supplied mask
% afterward, matching the selected Python contract before the first TOADS.
opts.T1_bname = get_basename(T1_file);
opts.working_dir = fullfile(working_dir, 'ToadsSeperate');
if ~isfolder(opts.working_dir), mkdir(opts.working_dir); end
n4_file = run_N4Correct(T1_file, opts, 'N4Correct', 'STEP2.1');
masked_n4_file = fullfile(working_dir, 'ToadsSeperate', 'T1_N4Corrected_masked.nii.gz');
helper = fullfile(packet_root, 'matlab', 'apply_brain_mask_nifti.py');
mask_cmd = sprintf('"%s" "%s" --image "%s" --mask "%s" --output "%s"', ...
    python_executable, helper, n4_file, brain_mask, masked_n4_file);
[mask_status, mask_output] = system(mask_cmd);
fprintf('%s\n', mask_output);
if mask_status ~= 0 || ~isfile(masked_n4_file)
    error('Applying the supplied mask to N4 output failed (status %d).', mask_status);
end

% Two selected TOADS passes with V3 Dura removal between them.
opts.T1_bname = get_basename(masked_n4_file);
opts.working_dir = fullfile(working_dir, 'ToadsSeperate');
opts.ToadsOutput = 'dura removal inputs';
[toads1_gm, toads1_wmfill, toads1_wmmask] = run_Toads( ...
    masked_n4_file, opts, 'Toads1stPass', 'STEP2.2');
[dura_t1, dura_mask] = run_RemoveDura( ...
    masked_n4_file, toads1_gm, toads1_wmfill, toads1_wmmask, ...
    opts, 'removeDura', 'STEP2.3');
opts.ToadsOutput = 'cruise inputs';
[toads2_gm, toads2_wmfill, toads2_wmmask] = run_Toads( ...
    dura_t1, opts, 'Toads2ndPass', 'STEP2.4');

% Continue the original BASE source sequence from transfer/crop through TGDM.
% The parent maCRUISEpipeline retains working_dir while run_Toads uses its
% local opts copy for ToadsSeperate; restore the same parent scope here.
opts.working_dir = working_dir;
transfer_orig = run_TransferToOriginSpace( ...
    Seg_file, {toads2_gm, toads2_wmfill, toads2_wmmask}, ...
    opts, 'TransferToOriginSpace', 'STEP3');
crop_files = [{WM_file, Mask_file, BrainMask_file, Seg_file}, transfer_orig];
[GM_crop_raw, cropped] = run_CropBaseOnGM( ...
    GM_file, crop_files, opts, 'Crop', 'STEP4');
WM_crop_raw = cropped{1};
Mask_crop_raw = cropped{2};
BrainMask_crop = cropped{3};
Seg_crop = cropped{4};
GM_toads = cropped{5};
WM_toads = cropped{6};
WMmask_toads = cropped{7};

opts.working_dir = working_dir;
[Boundary_file, ~] = run_GetMABoundary(Seg_crop, opts, 'NewCsfFiles', 'STEP5');
Boundary_file = run_CleanBoundaryByToads( ...
    Boundary_file, WMmask_toads, opts, 'NewCsfFiles', 'STEP5.5');
[GM_crop, WM_crop, Mask_crop] = update_WMandGM( ...
    GM_crop_raw, WM_crop_raw, Mask_crop_raw, BrainMask_crop, ...
    GM_toads, WM_toads, WMmask_toads, opts, 'UpdatedWMandGM', 'STEP6');
Boundary_file = run_CleanBoundaryByUpdateWM( ...
    Boundary_file, WM_crop, opts, 'NewCsfFiles', 'STEP5.6');

opts.update_Enhance = 0;
opts.WMIsoLevel = 0.5;
opts.GVF_Lambda = 1;
opts.IfGenerateIntermediate = 'true';
opts.NewCSFfile = Boundary_file;
opts.UseUpdatedCSFinACE = 2;
opts.MaWM = WM_crop_raw;
opts.refineTopoCorrectByMa = 1;
run_maCruiseSurf_pipeline(t1_raw, GM_file, GM_crop, WM_crop, Mask_crop, ...
    opts, 'maCruiseSurfSeperate', 'STEP7');

refine_dir = fullfile(working_dir, 'SegRefine');
refined_seg = run_SegmentationRefine(working_dir, seg_raw, refine_dir, opts);
outputs = struct('working_dir', working_dir, 'n4_t1', n4_file, ...
    'masked_n4_t1', masked_n4_file, 'toads1_gm', toads1_gm, ...
    'toads1_wmfill', toads1_wmfill, 'toads1_wmmask', toads1_wmmask, ...
    'dura_t1', dura_t1, 'dura_mask', dura_mask, 'toads2_gm', toads2_gm, ...
    'toads2_wmfill', toads2_wmfill, 'toads2_wmmask', toads2_wmmask, ...
    'refined_segmentation', refined_seg);
save(fullfile(working_dir, 'driver_outputs.mat'), 'outputs', 'opts');
fprintf('Selected-mask v3 BASE stages completed. Output directory: %s\n', working_dir);
end
