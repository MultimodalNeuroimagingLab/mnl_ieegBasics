%
%   Returns DKT labels to input xyz coordinates based on volume labelling
%   Requires vistasoft
%
%   labs = ieeg_getLabelXyzDKT(xyzs, FSdir)
%   [labs, labs_val] = ieeg_getLabelXyzDKT(xyzs, FSdir, rad)
%       
%       xyzs =              nx3 array of n [x, y, z] coordinates
%       FSdir =             path to FreeSurfer directory for the current subject; e.g. */derivatives/freesurfer/sub-<subject>
%       rad (optional) =    radius in mm to search for labels. 3 by default
%
%   Returns:
%       labs =              nx1 cell array of string labels, corresponding to each xyz coordinate
%       labs_val =          nx1 array of colortable values that correspond to the labels in labs
%
%   Adapted from 'mnl_ieegBasics/functions/ieeg_labelElectrodesDKT.m' on 11/12/2020 by HH
%

function [labs, labs_val] = ieeg_getLabelXyzDKT(xyzs, FSdir, rad)
    if nargin < 3, rad = 3; end % circle radius
    
    %% load files
    niDKT = niftiRead(fullfile(FSdir, 'mri', 'aparc.DKTatlas+aseg.nii.gz')); % labelled nifti

    [~, ~, colortable_DKT] = read_annotation(fullfile(FSdir, 'label', 'lh.aparc.DKTatlas.annot')); % hemisphere-agnostic; use 'l' by default

    subcortical_labels = readtable(fullfile(FSdir,'mri',...
        'aseg.auto_noCCseg.label_intensities.txt'),...
        'FileType','text','Delimiter',' ','ReadVariableNames',false,'HeaderLines',0);
    
    
    [roiCodes,roiNames] = dkt_labels();
    
    %%
    voxelsize = niDKT.pixdim;
    
    labs_val = zeros(size(xyzs, 1),1); % colortable values corresponding to each DKT label
    labs = cell(size(xyzs, 1),1);
    
    els_ind = niDKT.sto_ijk * [xyzs ones(size(xyzs, 1),1)]'; % converts coordinates to indices
    els_ind = round(els_ind(1:3,:)');
    
    for elec = 1:size(xyzs, 1) % loop across electrodes
        
        if isnan(xyzs(elec,1)) % electrode no position
            labs_val(elec) = nan;
            labs{elec} = 'n/a';
            
        else
            temp.electrode = zeros(size(niDKT.data));

            % box with ones
            minidata.elecbox = NaN(round(abs(rad*4/voxelsize(1))),...
                round(abs(rad*4/voxelsize(2))),...
                round(abs(rad*4/voxelsize(3))));
            % xyz
            [minidata.x,minidata.y,minidata.z] = ...
                ind2sub(size(minidata.elecbox),find(isnan(minidata.elecbox)));
            minidata.x = minidata.x*voxelsize(1);
            minidata.y = minidata.y*voxelsize(2);
            minidata.z = minidata.z*voxelsize(3);
            minidata.mean = [mean(minidata.x) mean(minidata.y) mean(minidata.z)];
            minidata.straal = sqrt((minidata.x-minidata.mean(1)).^2+... 
                (minidata.y-minidata.mean(2)).^2+... 
                (minidata.z-minidata.mean(3)).^2);
            minidata.elecbox(minidata.straal<rad) = 1;
            minidata.elecbox(minidata.straal==min(minidata.straal)) = 1;

            % set size box 
            temp.xsize = length(minidata.elecbox(:,1,1));
            temp.ysize = length(minidata.elecbox(1,:,1));
            temp.zsize = length(minidata.elecbox(1,1,:));

            if mod(temp.xsize,2)==0 % even
                temp.xmindefine = floor((length(minidata.elecbox(:,1,1))-1)/2);
                temp.xmaxdefine = ceil((length(minidata.elecbox(:,1,1))-1)/2);
            else
                temp.xmindefine = floor(length(minidata.elecbox(:,1,1))/2);
                temp.xmaxdefine = floor(length(minidata.elecbox(:,1,1))/2);
            end 
            if mod(temp.ysize,2)==0 % even
                temp.ymindefine = floor((length(minidata.elecbox(1,:,1))-1)/2);
                temp.ymaxdefine = ceil((length(minidata.elecbox(1,:,1))-1)/2);
            else
                temp.ymindefine = floor(length(minidata.elecbox(1,:,1))/2);
                temp.ymaxdefine = floor(length(minidata.elecbox(1,:,1))/2);
            end 
            if mod(temp.zsize,2)==0 % even
                temp.zmindefine = floor((length(minidata.elecbox(1,1,:))-1)/2);
                temp.zmaxdefine = ceil((length(minidata.elecbox(1,1,:))-1)/2);
            else
                temp.zmindefine = floor(length(minidata.elecbox(1,1,:))/2);
                temp.zmaxdefine = floor(length(minidata.elecbox(1,1,:))/2);
            end

            diffxmin=1;
            diffxmax=0;
            diffymin=1;
            diffymax=0;
            diffzmin=1;
            diffzmax=0;

            if els_ind(elec,1)-temp.xmindefine<=0
                disp('ERROR: electrode x < minimal x of image')
                diffxmin=abs(els_ind(elec,1)-temp.xmindefine)+1;
                temp.xmindefine=temp.xmindefine-diffxmin;
            end
            if els_ind(elec,2)-temp.ymindefine<=0
                disp('ERROR: electrode y < minimal y of image')
                diffymin=abs(els_ind(elec,2)-temp.ymindefine)+1;
                temp.ymindefine=temp.ymindefine-diffymin;%return
            end
            if els_ind(elec,3)-temp.zmindefine<=0
                disp('ERROR: electrode z < minimal z of image')
                diffzmin=abs(els_ind(elec,3)-temp.zmindefine)+1;
                temp.zmindefine=temp.zmindefine-diffzmin;%return
            end
            if els_ind(elec,1)+temp.xmaxdefine>length(temp.electrode(:,1,1))
                disp('ERROR: electrode x > maximal x of image')
                diffxmax=els_ind(elec,1)+temp.xmaxdefine-length(temp.electrode(:,1,1));
                temp.xmaxdefine=temp.xmaxdefine-(els_ind(elec,1)+temp.xmaxdefine-length(temp.electrode(:,1,1)));
            end
            if els_ind(elec,2)+temp.ymaxdefine>length(temp.electrode(1,:,1))
                disp('ERROR: electrode y > maximal y of image')
                diffymax=els_ind(elec,2)+temp.ymaxdefine-length(temp.electrode(1,:,1));
                temp.ymaxdefine=temp.ymaxdefine-(els_ind(elec,2)+temp.ymaxdefine-length(temp.electrode(1,:,1)));
            end
            if els_ind(elec,3)+temp.zmaxdefine>length(temp.electrode(1,1,:))
                disp(['ERROR: for electrode ' int2str(elec) ' z > maximal z of image'])
                diffzmax=els_ind(elec,3)+temp.zmaxdefine-length(temp.electrode(1,1,:));
                temp.zmaxdefine=temp.zmaxdefine-(els_ind(elec,3)+temp.zmaxdefine-length(temp.electrode(1,1,:)));
                %return
                %break
            end  

            minidata.elecbox = minidata.elecbox(diffxmin:end-diffxmax,diffymin:end-diffymax,diffzmin:end-diffzmax);

            % for each electrode draw circle
            temp.electrode(els_ind(elec,1)-temp.xmindefine:els_ind(elec,1)+temp.xmaxdefine,...
                els_ind(elec,2)-temp.ymindefine:els_ind(elec,2)+temp.ymaxdefine,...
                els_ind(elec,3)-temp.zmindefine:els_ind(elec,3)+temp.zmaxdefine) = ...
                minidata.elecbox;

            % get the label using the mode:
            thisLabel = double(mode(niDKT.data(temp.electrode==1)));

            % if <60% white matter, use the most common label of other voxels
            % 2 is label for left WM, 41 is label for right WM
            if thisLabel == 2
                wm_fraction = length(find(niDKT.data(temp.electrode==1)==2))./...
                    length(niDKT.data(temp.electrode==1));
                if wm_fraction<.90 % close to WM, but use other label
                    temp_labels = niDKT.data(temp.electrode==1);
                    temp_labels(temp_labels==2) = [];
                    thisLabel = mode(temp_labels);
                end
            elseif thisLabel == 41
                wm_fraction = length(find(niDKT.data(temp.electrode==1)==41))./...
                    length(niDKT.data(temp.electrode==1));
                if wm_fraction<.90 % close to WM, but use other label
                    temp_labels = niDKT.data(temp.electrode==1);
                    temp_labels(temp_labels==41) = [];
                    thisLabel = mode(temp_labels);
                end
            elseif thisLabel == 0 % if CSF, use other label if 10% of voxels have another label
                other_fraction = length(find(niDKT.data(temp.electrode==1)==0))./...
                    length(niDKT.data(temp.electrode==1));
                if other_fraction<.90 % 
                    temp_labels = niDKT.data(temp.electrode==1);
                    temp_labels(ismember(temp_labels,[0 8 47])) = []; % remove zero labels
                    thisLabel = mode(temp_labels);
                end
            elseif thisLabel == 8 || thisLabel == 47  % Cerebellum
                other_fraction = length(find(ismember(niDKT.data(temp.electrode==1),[8 47])))./...
                    length(niDKT.data(temp.electrode==1));
                if other_fraction<.90 
                    temp_labels = niDKT.data(temp.electrode==1);
                    temp_labels(ismember(temp_labels,[0 8 47])) = []; % remove 8 and 47 cerebellar labels
                    thisLabel = mode(temp_labels);
                end
            end

            labs_val(elec) = double(thisLabel);

            % get the text
            if labs_val(elec) > 0 && labs_val(elec) < 200
                if ~isempty(find(subcortical_labels.Var1==thisLabel,1))
                    labs{elec} = subcortical_labels.Var2{subcortical_labels.Var1==thisLabel};
                else
                    labs{elec} = 'n/a';
                end
            elseif labs_val(elec)>1000
                labs{elec} = roiNames{roiCodes==thisLabel};
            else
                labs{elec} = 'n/a';
            end
        end
    end
    
end