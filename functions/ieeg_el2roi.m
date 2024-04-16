function [ni] = ieeg_el2roi(elecTable,fname,el_name,circleradius,outdir)
% input: 
%   els: table with electrodes .name and .x .y .z coordinates (native space)
%   fname: name of the nifti in which space to white electrodes
%   elname: electrode names to write ROIs for, matching els.name
%   circleradius: radius in milimeter of ROI
%   outdir: directory to write rois
%
% output:
%   nifti file with binary mask, only for last electrode written
%
% default: 
%   [ni] = ieeg_el2roi(elecTable,fname,el_name,circleradius,outdir)
%
% DH, GOV, LYR 2024, Mayo Clinic
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%    
%   Version 1.1.0, released 26-11-2009
%


%% select resliced image for electrodes

if isempty(fname)
    disp('select T1 nifti file')
    [FILENAME, PATHNAME] = ...
        uigetfile({'*.nii;*.nii.gz','nifti-files (*.nii, *.nii.gz)'},'select resliced image for electrodes');
    fname = [PATHNAME,FILENAME];
end
data = niftiRead(fname);

% convert electrodes from native 2 indices
els = [elecTable.x elecTable.y elecTable.z];
els_ind = data.sto_ijk * [els ones(size(els,1),1)]';
els_ind = round(els_ind(1:3,:)');

voxelsize = data.pixdim(1:3);

%%%% define box with electrode (circle of ones)
% box with ones
minidata.elecbox = zeros(round(abs(circleradius*3/voxelsize(1))),...
    round(abs(circleradius*3/voxelsize(2))),...
    round(abs(circleradius*3/voxelsize(3))));
% xyz within minibox
[minidata.x,minidata.y,minidata.z] = ...
    ind2sub(size(minidata.elecbox),find(minidata.elecbox==0));
minidata.x = minidata.x*voxelsize(1);
minidata.y = minidata.y*voxelsize(2);
minidata.z = minidata.z*voxelsize(3);
minidata.mean = [mean(minidata.x) mean(minidata.y) mean(minidata.z)];
minidata.straal = sqrt((minidata.x-minidata.mean(1)).^2+... 
    (minidata.y-minidata.mean(2)).^2+... 
    (minidata.z-minidata.mean(3)).^2);
minidata.elecbox(minidata.straal<circleradius) = 1; % circle
minidata.elecbox(minidata.straal==min(minidata.straal)) = 1; % center
% get size box 
x_size = size(minidata.elecbox,1);
y_size = size(minidata.elecbox,2);
z_size = size(minidata.elecbox,3);

els_to_write = find(ismember(elecTable.name,el_name));

for elec_ind = 1:length(els_to_write)
    elec = els_to_write(elec_ind);
    electrode_name = elecTable.name{elec};

    if ~isnan(els_ind(elec,:))
    
        temp.electrode = zeros(size(data.data,1),size(data.data,2),size(data.data,3));
        
        %%%% make sure box fits in data
        if mod(x_size,2)==0 % even
            temp.xmindefine = floor((size(minidata.elecbox,1)-1)/2);
            temp.xmaxdefine = ceil((size(minidata.elecbox,1)-1)/2);
        else
            temp.xmindefine = floor(size(minidata.elecbox,1)/2);
            temp.xmaxdefine = floor(size(minidata.elecbox,1)/2);
        end
        if mod(y_size,2)==0 % even
            temp.ymindefine = floor((size(minidata.elecbox,2)-1)/2);
            temp.ymaxdefine = ceil((size(minidata.elecbox,2)-1)/2);
        else
            temp.ymindefine = floor(size(minidata.elecbox,2)/2);
            temp.ymaxdefine = floor(size(minidata.elecbox,2)/2);
        end
        if mod(z_size,2)==0 % even
            temp.zmindefine = floor((size(minidata.elecbox,3)-1)/2);
            temp.zmaxdefine = ceil((size(minidata.elecbox,3)-1)/2);
        else
            temp.zmindefine = floor(size(minidata.elecbox,3)/2);
            temp.zmaxdefine = floor(size(minidata.elecbox,3)/2);
        end
        
        % check indices:
        diffxmin = 1;
        diffxmax = 0;
        diffymin = 1;
        diffymax = 0;
        diffzmin = 1;
        diffzmax = 0;
        
        if els_ind(elec,1)-temp.xmindefine<=0
            warning(['for electrode ' int2str(elec) ' x < minimal x of image'])
            diffxmin = abs(els_ind(elec,1)-temp.xmindefine)+1;
            temp.xmindefine = temp.xmindefine-diffxmin;
        end
        if els_ind(elec,2)-temp.ymindefine<=0
            warning(['for electrode ' int2str(elec) ' y < minimal y of image'])
            diffymin = abs(els_ind(elec,2)-temp.ymindefine)+1;
            temp.ymindefine = temp.ymindefine-diffymin;%return
        end
        if els_ind(elec,3)-temp.zmindefine<=0
            warning(['for electrode ' int2str(elec) ' z < minimal z of image'])
            diffzmin = abs(els_ind(elec,3)-temp.zmindefine)+1;
            temp.zmindefine = temp.zmindefine-diffzmin;%return
        end
        if els_ind(elec,1)+temp.xmaxdefine>size(temp.electrode,1)
            warning(['for electrode ' int2str(elec) ' x > maximal x of image'])
            diffxmax = els_ind(elec,1)+temp.xmaxdefine-length(temp.electrode(:,1,1));
            temp.xmaxdefine = temp.xmaxdefine-(els_ind(elec,1)+temp.xmaxdefine-length(temp.electrode(:,1,1)));
        end
        if els_ind(elec,2)+temp.ymaxdefine>size(temp.electrode,2)
            warning(['for electrode ' int2str(elec) ' y > maximal y of image'])
            diffymax = els_ind(elec,2)+temp.ymaxdefine-length(temp.electrode(1,:,1));
            temp.ymaxdefine = temp.ymaxdefine-(els_ind(elec,2)+temp.ymaxdefine-length(temp.electrode(1,:,1)));
        end
        if els_ind(elec,3)+temp.zmaxdefine>size(temp.electrode,3)
            warning(['for electrode ' int2str(elec) ' z > maximal z of image'])
            diffzmax = els_ind(elec,3)+temp.zmaxdefine-length(temp.electrode(1,1,:));
            temp.zmaxdefine = temp.zmaxdefine-(els_ind(elec,3)+temp.zmaxdefine-length(temp.electrode(1,1,:)));
        end
        % clip elecbox if it does not fit in the nifti (at edges)
        minidata.elecbox = minidata.elecbox(diffxmin:end-diffxmax,diffymin:end-diffymax,diffzmin:end-diffzmax);
        
        % first get values from coordinates where we want to add an
        % electrode:
        previous_values_minibox = temp.electrode(els_ind(elec,1)-temp.xmindefine:els_ind(elec,1)+temp.xmaxdefine,...
            els_ind(elec,2)-temp.ymindefine:els_ind(elec,2)+temp.ymaxdefine,...
            els_ind(elec,3)-temp.zmindefine:els_ind(elec,3)+temp.zmaxdefine);
        
        new_minibox = max(cat(4,minidata.elecbox,previous_values_minibox),[],4);
        
        % for each electrode draw circle
        temp.electrode(els_ind(elec,1)-temp.xmindefine:els_ind(elec,1)+temp.xmaxdefine,...
            els_ind(elec,2)-temp.ymindefine:els_ind(elec,2)+temp.ymaxdefine,...
            els_ind(elec,3)-temp.zmindefine:els_ind(elec,3)+temp.zmaxdefine) = ...
            new_minibox;
        clear new_minibox 
        if temp.electrode(els_ind(elec,1),els_ind(elec,2),els_ind(elec,3))~=1
            disp(['error' int2str(elec)])
        end
        
        % and save nifti with electrode
        
        % initialize output
        ni = data;
        ni.data = temp.electrode;
        ni.dim = ni.dim(1:3); ni.pixdim = ni.pixdim(1:3);
        ni.scl_inter = 0; % always set intercept to zero to have electrodes be 1 and zero background
        
        outputname = fullfile(outdir,['ElectrodeRoi_' electrode_name '.nii']);
    
        disp(strcat(['saving ' outputname]));
        % save the data
        niftiWrite(ni,outputname)
    end
end