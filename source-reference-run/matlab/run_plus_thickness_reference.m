function [volume, thickness_mat] = run_plus_thickness_reference(inner, outer, inner_surf, outer_surf, opts)
% Original scientific thickness operations, retaining logs and measured mesh.
% run_GetThickNessCRUISE removes those artifacts, so call its original helpers
% directly. Do not derive thickness from Python products.
root = fullfile(opts.working_dir, 'Thickness');
if isfolder(root), error('Thickness output already exists: %s', root); end
mkdir(root);
opts.output_dir = root;
opts.working_dir = root;
cmd = run_GetThickNessCRUISE_jist(inner, outer, opts);
fid = fopen(fullfile(root, 'command.txt'), 'w');
if fid < 0, error('Cannot save thickness command'); end
fwrite(fid, cmd, 'char'); fclose(fid);
[status, message] = system(cmd);
fprintf('%s\n', message);
fid = fopen(fullfile(root, 'thickness.log'), 'w');
if fid < 0, error('Cannot save thickness log'); end
fwrite(fid, message, 'char'); fclose(fid);
if status ~= 0, error('Original volume thickness failed with status %d', status); end
files = dir(fullfile(root, '*_thick.nii.gz'));
if numel(files) ~= 1, error('Missing or ambiguous original thickness volume'); end
volume = fullfile(files(1).folder, files(1).name);
measured = run_GetThicknessByClosestPoints(inner_surf, [], outer_surf, opts, 'surface');
if ~isfile(measured), error('Missing original measured thickness surface'); end
[mean_thickness, std_thickness, whole_thickness] = WholeBrainThickness(measured);
thickness_mat = fullfile(root, 'surface-thickness.mat');
save(thickness_mat, 'mean_thickness', 'std_thickness', 'whole_thickness');
end
