% Run with MATLAB from this fixture directory after adding source paths.
% Original helpers: MaCRUISE/pipeline and masimatlab/trunk/utils/matlab.
root=fileparts(mfilename('fullpath'));
opts=struct('InnerCoef',0.2,'OuterCoef',0.2,'InnerCoefCentral',0.2,'OuterCoefCentral',0.2,'IfDealwithCentral',1);
MAsegRefineBySurf(fullfile(root,'segmentation.nii.gz'),fullfile(root,'inner.nii.gz'),fullfile(root,'outer.nii.gz'),fullfile(root,'matlab-refined.nii.gz'),opts);
