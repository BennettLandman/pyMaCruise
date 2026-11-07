function outputs = run_plus_selected_mask_reference(maCRUISE_source_dir, mipav_dir, multiatlas_code_dir, n4_binary, python_executable, config_file)
% Run selected v3+ PLUS stages with the supplied T1-grid brain mask.
%
% This run-only driver follows run_maCRUISEplus.m and maCRUISEpluspipeline.m
% while bypassing BET/SPECTRE and all visualization/MNI/FreeSurfer exports.
% FLAIR registration remains the original selected NiftyReg/Aladin rigid
% wrapper (rigOnly=1). It then runs original N4, EM, Dura, LesionTOADS,
% tissue-fusion, ACE/GVF/TGDM, and refinement stages with the supplied mask.
% It does not modify upstream source.

if nargin < 5 || isempty(python_executable), python_executable = 'python3'; end
if nargin < 4 || isempty(n4_binary), error('Provide ANTs N4BiasFieldCorrection executable path.'); end
if nargin < 3 || isempty(multiatlas_code_dir), error('Provide full-multi-atlas checkout path.'); end
if nargin < 2 || isempty(mipav_dir), error('Provide MIPAV install path with selected JIST plugins and atlas.'); end
if nargin < 1 || isempty(maCRUISE_source_dir), error('Provide MaCRUISE MATLAB source checkout path.'); end

driver_dir = fileparts(mfilename('fullpath'));
packet_root = fileparts(driver_dir);
pipeline_dir = fullfile(maCRUISE_source_dir, 'pipeline');
input_dir = fullfile(packet_root, 'inputs', '100260');
config = struct();
if nargin >= 6 && ~isempty(config_file)
    config = jsondecode(fileread(config_file));
    if isfield(config, 'input_dir'), input_dir = config.input_dir; end
end
t1_raw = fullfile(input_dir, 'T1.nii.gz');
flair_raw = fullfile(input_dir, 'FLAIR.nii.gz');
seg_raw = fullfile(input_dir, 'MaSeg.nii.gz');
brain_mask = fullfile(input_dir, 'brain_mask.nii.gz');
for file = {t1_raw, flair_raw, seg_raw, brain_mask}
    if ~isfile(file{1}), error('PLUS packet is missing input: %s', file{1}); end
end

stamp = datestr(now, 'yyyymmdd_HHMMSS');
working_dir = fullfile(packet_root, 'outputs', ['plus-selected-mask-v3-' stamp]);
if isfield(config, 'output_dir'), working_dir = config.output_dir; end
if isfolder(working_dir), error('Refusing to reuse output directory: %s', working_dir); end
mkdir(working_dir);
diary(fullfile(working_dir, 'matlab-driver.log'));
diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

masi_fusion_dir = fullfile(multiatlas_code_dir, 'masi-fusion', 'src');
if ~isfolder(pipeline_dir) || ~isfolder(masi_fusion_dir)
    error('Could not find MaCRUISE pipeline or MASI fusion source directories.');
end
addpath(genpath(masi_fusion_dir));
addpath(genpath(pipeline_dir));
required_functions = {'load_untouch_nii_gz', 'save_untouch_nii_gz', ...
    'tabulate', 'run_cmd_single', 'FlairRegistration', 'WholeBrainThickness'};
for function_index = 1:numel(required_functions)
    if isempty(which(required_functions{function_index}))
        error('Missing original MATLAB dependency: %s', required_functions{function_index});
    end
end
for executable = {n4_binary, fullfile(mipav_dir, 'jre', 'bin', 'java'), ...
        fullfile(multiatlas_code_dir, 'niftyreg', 'bin', 'reg_aladin')}
    if ~isfile(executable{1}), error('Missing source runtime executable: %s', executable{1}); end
end

% Original run_maCRUISEplus registration defaults and selected source backend.
opts.runmethod = 'single';
opts.jvmmemval = '8G';
opts.mipav = mipav_dir;
opts.plugindir = fullfile(mipav_dir, 'plugins');
opts.workingrootfolder = working_dir;
opts.working_dir = working_dir;
opts.multiatlas_code_dir = multiatlas_code_dir;
opts.N4correctBIN_dir = n4_binary;
opts.rigOnly = 1;
opts.interp = [];
if ~isfield(config, 'atlasem_descriptor')
    error('Provide source-config.json with the owner-selected native ibsr11objV3 atlasem_descriptor.');
end
opts.atlasem_descriptor = config.atlasem_descriptor;
required_atlases = { ...
    opts.atlasem_descriptor, ...
    fullfile(opts.plugindir, 'Atlas', 'lesionToads-atlas-2012', 'cruise-atlas-12obj-lesiontoads2012.txt'), ...
    fullfile(opts.plugindir, 'Atlas', 'lesionToads-atlas-2012', 'cruise-atlas-12obj-toads2012.txt'), ...
    fullfile(opts.plugindir, 'Atlas', 'lesionToads-atlas-2012', 'cruise-atlas-10obj-toads2010.txt')};
for atlas_index = 1:numel(required_atlases)
    if ~isfile(required_atlases{atlas_index})
        error('Missing selected PLUS atlas in MIPAV plugins: %s', required_atlases{atlas_index});
    end
end

% run_maCRUISEplus first rigidly registers FLAIR to T1 with NiftyReg/Aladin.
flair_reg_dir = fullfile(working_dir, 'FlairReg');
if ~isfolder(flair_reg_dir), mkdir(flair_reg_dir); end
flair_registered = FlairRegistration(t1_raw, flair_raw, flair_reg_dir, opts);

% Match the original plus pipeline's raw preparation and MaSeg separation.
plus_dir = fullfile(working_dir, 'MaCRUISEplus');
if ~isfolder(plus_dir), mkdir(plus_dir); end
[T1_file, Seg_file, Flair_file] = prepare_raw_plus( ...
    t1_raw, flair_registered, seg_raw, plus_dir, 'Raw', 'STEP0');
opts.rawseg = Seg_file;
[GM_file, WM_file, Mask_file, BrainMask_file] = cortical_seperation( ...
    Seg_file, plus_dir, 'CorticalSeperation', 'STEP1');

% Original selected source uses N4 on T1 and registered FLAIR before
% SPECTRE. Replace only the SPECTRE mask with the supplied T1-grid mask.
lesion_dir = fullfile(plus_dir, 'LesionToadsSeperate');
if ~isfolder(lesion_dir), mkdir(lesion_dir); end
opts.working_dir = lesion_dir;
opts.N3orN4 = 'N4';
opts.T1_bname = get_basename(T1_file);
t1_n4 = run_N4Correct(T1_file, opts, 'N4Correct', 'STEP2.1_T1');
opts.T1_bname = get_basename(Flair_file);
flair_n4 = run_N4Correct(Flair_file, opts, 'N4Correct', 'STEP2.1_FLAIR');

mask_helper = fullfile(packet_root, 'matlab', 'apply_brain_mask_nifti.py');
t1_masked = fullfile(lesion_dir, 'T1_N4Corrected_user_mask.nii.gz');
flair_masked = fullfile(lesion_dir, 'FLAIR_N4Corrected_user_mask.nii.gz');
apply_mask(python_executable, mask_helper, t1_n4, brain_mask, t1_masked);
apply_mask(python_executable, mask_helper, flair_n4, brain_mask, flair_masked);
t1_masked = removeNegative(t1_masked, opts, 'user_masked', 'T1');
flair_masked = removeNegative(flair_masked, opts, 'user_masked', 'FLAIR');

% Continue the original source sequence: AtlasEM -> Dura -> LesionTOADS.
t1_bname = get_basename(T1_file);
opts.T1_bname = t1_bname;
[EM_gm, EM_wmfill, EM_wmmask] = run_plus_atlasem_reference( ...
    t1_masked, opts, 'EMsegmentation', 'STEP2.7');
[t1_dura, dura_mask] = run_RemoveDura( ...
    t1_masked, EM_gm, EM_wmfill, EM_wmmask, opts, 'removeDura_init', 'STEP2.8');
opts.T1_bname = get_basename(flair_masked);
flair_dura_stripped = fullfile(lesion_dir, 'StripOrigin_Dura', ...
    [opts.T1_bname '_strip.nii.gz']);
apply_mask(python_executable, mask_helper, flair_masked, dura_mask, flair_dura_stripped);
opts.T1_bname = t1_bname;
t1_dura = moveFile(t1_dura, opts, 'removeDura_final', 'T1_dura_stripped');
flair_dura = moveFile(flair_dura_stripped, opts, 'removeDura_final', 'Flair_dura_stripped');
opts.LesionToadsOutput = 'dura removal inputs';
opts.LesionToadsCorrect = 'false';
[Toads_gm, Toads_wmfill, Toads_wmmask] = run_LesionToads( ...
    t1_dura, flair_dura, opts, 'LesionTOADS', 'STEP2.10');

% Original plus stages 3-7: transfer, crop, boundary, fusion, CRUISE core.
% The parent maCRUISEpluspipeline keeps working_dir at plus_dir after its
% nested LesionTOADS helper returns; restore that scope before stage 3.
opts.working_dir = plus_dir;
transfer_orig = run_TransferToOriginSpace( ...
    Seg_file, {Toads_gm, Toads_wmfill, Toads_wmmask}, ...
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
opts.working_dir = plus_dir;
opts.ErodeMethod = '3D';
opts.OverLookLabels = '120,121,202,203,156,157';
[Boundary_file, ~] = run_GetMABoundary(Seg_crop, opts, 'NewCsfFiles', 'STEP5');
Boundary_file = run_CleanBoundaryByToads( ...
    Boundary_file, WMmask_toads, opts, 'NewCsfFiles', 'STEP5.5');
[GM_crop, WM_crop, Mask_crop] = update_WMandGM( ...
    GM_crop_raw, WM_crop_raw, Mask_crop_raw, BrainMask_crop, ...
    GM_toads, WM_toads, WMmask_toads, opts, 'UpdatedWMandGM', 'STEP6');

opts.working_dir = plus_dir;
Boundary_file = run_CleanBoundaryByUpdateWM( ...
    Boundary_file, WM_crop, opts, 'NewCsfFiles', 'STEP5.6');
opts.update_Enhance = 0;
opts.WMIsoLevel = 0.5;
opts.GVF_Lambda = 0.2;
opts.IfGenerateIntermediate = 'true';
opts.NewCSFfile = Boundary_file;
opts.UseUpdatedCSFinACE = 3;
opts.MaWM = WM_crop_raw;
opts.refineTopoCorrectByMa = 1;
run_CruisePlusSurf_pipeline(t1_raw, GM_file, GM_crop, WM_crop, Mask_crop, ...
    opts, 'CruisePlusSurfSeperate', 'STEP7');

% The top-level wrapper always performs this selected original final stage.
refine_dir = fullfile(plus_dir, 'SegRefine');
refined_seg = run_SegmentationRefine(plus_dir, seg_raw, refine_dir, opts);
level_dir = fullfile(plus_dir, 'CruisePlusSurfSeperate', 'TGDMinvcrop', 'uncroped');
inner_level = fullfile(level_dir, 'target_image_GMimg_innerVol.nii.gz');
outer_level = fullfile(level_dir, 'target_image_GMimg_outerVol.nii.gz');
central_level = fullfile(level_dir, 'target_image_GMimg_centralVol.nii.gz');
inner_surface = fullfile(level_dir, 'target_image_GMimg_innerSurf.vtk');
outer_surface = fullfile(level_dir, 'target_image_GMimg_outerSurf.vtk');
central_surface = fullfile(level_dir, 'target_image_GMimg_centralSurf.vtk');
opts.working_dir = plus_dir;
[thickness_volume, thickness_mat] = run_plus_thickness_reference( ...
    inner_level, outer_level, inner_surface, outer_surface, opts);
outputs = struct('working_dir', working_dir, 'input_dir', input_dir, ...
    'atlasem_descriptor', opts.atlasem_descriptor, 'registered_flair', flair_registered, ...
    't1_n4', t1_n4, 'flair_n4', flair_n4, 't1_masked', t1_masked, ...
    'flair_masked', flair_masked, 'atlasem_gm', EM_gm, ...
    'atlasem_wmfill', EM_wmfill, 'atlasem_wmmask', EM_wmmask, ...
    't1_dura', t1_dura, 'dura_mask', dura_mask, ...
    'flair_dura', flair_dura, 'toads_gm', Toads_gm, ...
    'toads_wmfill', Toads_wmfill, 'toads_wmmask', Toads_wmmask, ...
    'refined_segmentation', refined_seg, ...
    'inner_level_set', inner_level, 'central_level_set', central_level, ...
    'outer_level_set', outer_level, 'inner_surface', inner_surface, ...
    'central_surface', central_surface, 'outer_surface', outer_surface, ...
    'thickness_volume', thickness_volume, 'surface_thickness_mat', thickness_mat);
for product = {refined_seg, inner_level, central_level, outer_level, ...
        inner_surface, central_surface, outer_surface, thickness_volume, thickness_mat}
    if ~isfile(product{1}), error('Missing final PLUS product: %s', product{1}); end
end
save(fullfile(working_dir, 'driver_outputs.mat'), 'outputs', 'opts');
fid = fopen(fullfile(working_dir, 'driver_outputs.json'), 'w');
if fid < 0, error('Could not write driver_outputs.json'); end
fwrite(fid, jsonencode(outputs), 'char'); fclose(fid);
fprintf('Selected-mask v3+ PLUS stages completed. Output directory: %s\n', working_dir);
end

function apply_mask(python_executable, helper, image, mask, output)
cmd = sprintf('"%s" "%s" --image "%s" --mask "%s" --output "%s"', ...
    python_executable, helper, image, mask, output);
[status, message] = system(cmd);
fprintf('%s\n', message);
if status ~= 0 || ~isfile(output)
    error('Applying a supplied PLUS mask failed (status %d): %s', status, output);
end
end
