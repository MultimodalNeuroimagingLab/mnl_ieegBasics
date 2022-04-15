function [t_new] = ieeg_labelElectrodesDestrieux(FSdir, electrodes_tsv_name, saveNew, circleradius)
%
% this script labels electrodes based on the Destrieux Atlas from
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
% Edited 11/13/2020 to call ieeg_getLabelXyzDestrieux.m
%
%% 

if nargin < 4 || isempty(circleradius), circleradius = 3; end
if nargin < 3 || isempty(saveNew), saveNew = false; end

% load electrode positions
loc_info = readtable(electrodes_tsv_name, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

elecmatrix = [loc_info.x loc_info.y loc_info.z];

[Destrieux_label_text, Destrieux_label] = ieeg_getLabelXyzDestrieux(elecmatrix, FSdir, circleradius); % get Destrieux labels

% append/edit Destrieux columns to loc_info
if ~ismember('Destrieux_label', loc_info.Properties.VariableNames)
    tDes = table(Destrieux_label, Destrieux_label_text);
    t_new = [loc_info tDes(:,{'Destrieux_label','Destrieux_label_text'})];
else
    t_new = loc_info;
    t_new.Destrieux_label = Destrieux_label;
    t_new.Destrieux_label_text = Destrieux_label_text;
end

if saveNew==1
    t_new = bids_tsv_nan2na(t_new);
    writetable(t_new, electrodes_tsv_name, 'FileType','text','Delimiter','\t'); 
end