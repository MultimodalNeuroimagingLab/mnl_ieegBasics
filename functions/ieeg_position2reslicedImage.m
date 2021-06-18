function [ni,els,els_ind] = ieeg_position2reslicedImage(els,fname)
% input: 
% els [3x] matrix with x y z coordinates of x electrodes (native space)
% default: 
%
%     Copyright (C) 2009  D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
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

% DH, adjusting to make surround 0.5 and list electrode number 1:nr
% electrodes in the center. We do this after we get the box, to make sure
% that box values are not overwritten for small electrodes.

circleradius = 2;

%% example to load positions of electrodes
% data.elecname='C:\Users\dora\Documents\AiO\gridpatients\ctAnalysisPackage\data\maan\ct\electrodes_loc2_global.mat';
% grid=load(data.elecname);
% els=grid.elecmatrix;

%% select resliced image for electrodes

if isempty(fname)
    disp('select T1 nifti file')
    [FILENAME, PATHNAME] = ...
        uigetfile({'*.nii;*.nii.gz','nifti-files (*.nii, *.nii.gz)'},'select resliced image for electrodes');
    fname = [PATHNAME,FILENAME];
end
data = niftiRead(fname);

% convert electrodes from native 2 indices
els_ind = data.sto_ijk * [els ones(size(els,1),1)]';
els_ind = round(els_ind(1:3,:)');

temp.electrode = zeros(size(data.data));

voxelsize = data.pixdim;

%%%% define box with electrode (circle of ones, 2 at centre at size mm)
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
minidata.elecbox(minidata.straal<circleradius) = .5;
minidata.elecbox(minidata.straal==min(minidata.straal)) = .7;
% get size box 
temp.xsize = size(minidata.elecbox,1);
temp.ysize = size(minidata.elecbox,2);
temp.zsize = size(minidata.elecbox,3);



for elec = 1:size(els_ind,1)
    if ~isnan(els_ind(elec,:))
        
        %%%% make sure box fits in data
        if mod(temp.xsize,2)==0 % even
            temp.xmindefine = floor((size(minidata.elecbox,1)-1)/2);
            temp.xmaxdefine = ceil((size(minidata.elecbox,1)-1)/2);
        else
            temp.xmindefine = floor(size(minidata.elecbox,1)/2);
            temp.xmaxdefine = floor(size(minidata.elecbox,1)/2);
        end
        if mod(temp.ysize,2)==0 % even
            temp.ymindefine = floor((size(minidata.elecbox,2)-1)/2);
            temp.ymaxdefine = ceil((size(minidata.elecbox,2)-1)/2);
        else
            temp.ymindefine = floor(size(minidata.elecbox,2)/2);
            temp.ymaxdefine = floor(size(minidata.elecbox,2)/2);
        end
        if mod(temp.zsize,2)==0 % even
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
        if temp.electrode(els_ind(elec,1),els_ind(elec,2),els_ind(elec,3))~=0.7
            disp(['error' int2str(elec)])
        end
        
    end
end

el_nr = 0;
for elec = 1:size(els_ind,1)
    if ~isnan(els_ind(elec,:)) 
        % if electrode fits in scan
        if els_ind(elec,1) > 0 && els_ind(elec,2)>0 && els_ind(elec,3)>0 &&...
                els_ind(elec,1)<size(temp.electrode,1) && els_ind(elec,2)<size(temp.electrode,2) && els_ind(elec,3)<size(temp.electrode,3)
            % add electrode number in the center
            temp.electrode(els_ind(elec,1),els_ind(elec,2),els_ind(elec,3)) = elec;
            el_nr = el_nr+1;
        else
            disp(['missed el nr ' int2str(elec)])
        end
    end
end

%% and save new data:

% initialize output
ni = data;
ni.data = temp.electrode;
ni.scl_inter = 0; % always set intercept to zero to have electrodes be 1 and zero background

disp('select output dir')
DIRECTORYNAME = uigetdir(pwd, 'select output directory');

for filenummer = 1:100
    outputname = fullfile(DIRECTORYNAME,['Electrodes2Image_' int2str(filenummer) '.nii']);

    if ~exist(outputname,'file')>0
        disp(strcat(['saving ' outputname]));
        % save the data
        niftiWrite(ni,outputname)
        break
    end
end
