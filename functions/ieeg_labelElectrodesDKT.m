function [t_new] = ieeg_labelElectrodesDKT(FSdir, electrodes_tsv_name, circleradius)
%
% this script labels electrodes based on the DKT Atlas from
% Freesurfer
%
% requires vistasoft in the path
%
%%% Preperation step: 
%
% Tn the terminal go to the subjects freesurfer folder and type: 
%
%  mri_convert mri/aparc.a2009s+aseg.mgz mri/aparc.a2009s+aseg.nii.gz -rt nearest
%
% This converts the mri to the original T1 space, using nearest neighbor
% reslicing to keep one label for each voxel (no interpolation)
%
% Edited 11/13/2020 to call ieeg_getLabelXyzDKT.m
%
%% 

if nargin < 3 || isempty(circleradius), circleradius = 3; end

% load electrode positions
loc_info = readtable(electrodes_tsv_name, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

elecmatrix = [loc_info.x loc_info.y loc_info.z];

[DKT_label_text, DKT_label] = ieeg_getLabelXyzDKT(elecmatrix, FSdir, circleradius); % get DKT labels

% append/edit DKT columns to loc_info
if ~ismember('DKT_label', loc_info.Properties.VariableNames)
    tDes = table(DKT_label, DKT_label_text);
    t_new = [loc_info tDes(:,{'DKT_label','DKT_label_text'})];
else
    t_new = loc_info;
    t_new.DKT_label = DKT_label;
    t_new.DKT_label_text = DKT_label_text;
end

[file_dir,file_name] = fileparts(electrodes_tsv_name);
electrodes_tsv_name_new = fullfile(file_dir,[file_name '_DKT.tsv']);

writetable(t_new, electrodes_tsv_name_new, 'FileType','text','Delimiter','\t'); 