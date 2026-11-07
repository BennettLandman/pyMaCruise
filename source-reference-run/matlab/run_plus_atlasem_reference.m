function [gm, wmfill, wmmask] = run_plus_atlasem_reference(t1, opts, folder, step)
% Run the original AtlasEM command with an explicit owner-selected atlas.
% Only -inAtlas and the output destination differ from run_EM_jist. Original
% scientific options and classes remain unchanged; preserve command and logs.
if ~isfield(opts, 'atlasem_descriptor') || ~isfile(opts.atlasem_descriptor)
    error('Provide the native ibsr11objV3 atlasem_descriptor in opts.');
end
opts.output_dir = fullfile(opts.working_dir, folder);
if isfolder(opts.output_dir), error('AtlasEM output already exists: %s', opts.output_dir); end
mkdir(opts.output_dir);
cmd = run_EM_jist(t1, opts);
[start_idx, end_idx] = regexp(cmd, ' -inAtlas "[^"]+"');
if numel(start_idx) ~= 1, error('Expected one -inAtlas in original run_EM_jist command.'); end
cmd = [cmd(1:start_idx-1) ' -inAtlas "' opts.atlasem_descriptor '"' cmd(end_idx+1:end)];
fid = fopen(fullfile(opts.output_dir, 'command.txt'), 'w');
if fid < 0, error('Cannot save AtlasEM command'); end
fwrite(fid, cmd, 'char'); fclose(fid);
fprintf('%s: %s\n', step, cmd);
[status, message] = system(cmd);
fprintf('%s\n', message);
fid = fopen(fullfile(opts.output_dir, 'atlasem.log'), 'w');
if fid < 0, error('Cannot save AtlasEM log'); end
fwrite(fid, message, 'char'); fclose(fid);
if status ~= 0, error('Original AtlasEM failed with exit status %d', status); end
gm = unique_product(opts.output_dir, '*AtlasEM_gm*.nii.gz');
wmfill = unique_product(opts.output_dir, '*AtlasEM_wmFill*.nii.gz');
wmmask = unique_product(opts.output_dir, '*AtlasEM_wmMask*.nii.gz');
end

function path = unique_product(root, pattern)
files = dir(fullfile(root, '**', pattern));
if numel(files) ~= 1, error('Expected exactly one %s under %s, found %d', pattern, root, numel(files)); end
path = fullfile(files(1).folder, files(1).name);
end
