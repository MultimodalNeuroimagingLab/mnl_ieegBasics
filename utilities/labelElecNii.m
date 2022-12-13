%% This function saves a nifti image of electrode names used to label electrode leads when looking at CT-MRI overlay (i.e., in MRIcroGL)
% The electrode names nifti file is saved in the same directory as the input electrodes.tsv
% Names are 5hx3w (x3 slices thick) pixel letters on the axial slice only, centered on the 1st, 6th, 11th, and 16th electrodes of each lead (if the electrode position exists).
% The first letter of each electrode name (R or L) is omitted for brevity. If electrode names contain 2 letters (e.g. OC), then they are both shown with 1-pixel space inbetween.
% The intensity of all pixel characters saved to the nifti image is 1.
%
%   labelElecNii(niiPath, elecPath);
%   vol = labelElecNii(niiPath, elecPath);
%       niiPath =       char, path to T1w nifti brain volume. This is used to get 1) the dimensions of the brain volume (image.dim) and 2) the transformation matrix (image.mat)
%                           to go from voxel space <-> position space
%       elecPath =      char, path to the electrodes.tsv file. Must contain columns "name", "x", "y", and "z". Electrodes with non-nan positions are assumed to 
%                           start with 'L' or 'R', which is omitted when saving.
%
%   Returns:
%       vol =           ixjxk double. The volume with same dimensions as the brain volume where the labels are plotted.
%
%   Dependency: SPM12. Needs spm_vol.m
%
%   Potential problems:
%       - If the transformation matrix (image.mat) does not have an identity diagonal, then shearing/scaling may mess up the appearance of the letters.
%       - If there are inconsistent reflections in the image across subjects, then letters may appear upside down or left-right mirrored.
%
%   HH 2022/11
%
function vol = labelElecNii(niiPath, elecPath)

    % Request path of electrodes in order to know where to save
    electrodes = readtable(elecPath, 'FileType', 'text', 'Delimiter', '\t');
    
    % load nifti volume if necessary
    if ~exist('spm_vol.m','file')
        error('make sure spm is in your path')
    else
        nii = spm_vol(niiPath);
    end
    
    %% Configure attributes of nifti image to save
    
    s = dir(elecPath); % directory and filename of electrodes
    parts = split(s.name, '.');
    
    niiLabels = nii; % use input image as template for output image
    niiLabels.fname = fullfile(s.folder, sprintf('%s_labels.nii', parts{1})); % configure output filename
    
    %% Select electrode names and positions to save
    
    electrodes(isnan(electrodes.x), :) = []; % remove electrodes without positions (and stuff like EKG)
    elecNames = strip(upper(electrodes.name)); % some cleaning
    locs = [electrodes.x, electrodes.y, electrodes.z]; % xyz locations
    
    isDigit = @(x) x > 47 & x < 58; % returns true for char array elements that are digits (0 - 9)
    
    leadNames = cell(size(elecNames)); % Does not contain hemisphere or numbers
    leadNum = zeros(size(elecNames)); % Only contains the numbers
    
    for ii = 1:length(elecNames)
        n = elecNames{ii}(2:end); % remove the first letter (always L or R)
        leadNames{ii} = n(~isDigit(n)); % discard number
        leadNum(ii) = str2double(n(isDigit(n))); % only number
    end
    
    % keep only electrodes that are spaced 5 apart, starting at 1
    leadNames = leadNames(ismember(leadNum, [1, 6, 11, 16]));
    locs = locs(ismember(leadNum, [1, 6, 11, 16]), :);
    leadNum = leadNum(ismember(leadNum, [1, 6, 11, 16]));

    %% Add labels into nifti volume
    
    vol = zeros(nii.dim);
    
    % Transform locs to inds. Same as locs * inv(nii.mat'). Reverse: inds*nii.mat' = locs
    inds = round([locs, ones(size(locs, 1), 1)] / nii.mat');
    inds(:, 4) = []; % remove affine 1s
    
    for ii = 1:size(inds, 1)
        raster = getLetterRaster(leadNames{ii}); % get the rasterized letter, upside down
        raster = raster(end:-1:1, end:-1:1); % invert letter left-right and upside down (based on empirical observation)
        raster = raster'; % transpose so that left-right direction of letter lies in the x axis
        rasterVol = repmat(raster, 1, 1, 3); % stack 3 in the z direction
        
        if length(leadNames{ii}) == 1
            vol((inds(ii, 1)-1):(inds(ii, 1)+1), (inds(ii, 2)-2):(inds(ii, 2)+2), (inds(ii, 3)-1):(inds(ii, 3)+1)) = rasterVol;
        else
            vol((inds(ii, 1)-3):(inds(ii, 1)+3), (inds(ii, 2)-2):(inds(ii, 2)+2), (inds(ii, 3)-1):(inds(ii, 3)+1)) = rasterVol;
        end
    end
    
    %% Save nifti volume
    
    spm_write_vol(niiLabels, vol);
    
end


function raster = getLetterRaster(letters)
    
    letters = upper(letters);

    if length(letters) == 2 % rasterize each letter separately and concatenate horizontally
        raster = cat(2, [getLetterRaster(letters(1)), zeros(5, 1, 1), getLetterRaster(letters(2))]); % 2-pixel space in between
        return
    elseif length(letters) > 2
        error('Length of letters to rasterize must be 1 or 2');
    end
    
    switch letters
        case 'A'
            raster = [0 1 0; 1 0 1; 1 1 1; 1 0 1; 1 0 1];
        case 'B'
            raster = [1 1 0; 1 0 1; 1 1 0; 1 0 1; 1 1 0];
        case 'C'
            raster = [0 1 1; 1 0 0; 1 0 0; 1 0 0; 0 1 1];
        case 'D'
            raster = [1 1 0; 1 0 1; 1 0 1; 1 0 1; 1 1 0];
        case 'E'
            raster = [1 1 1; 1 0 0; 1 1 0; 1 0 0; 1 1 1];
        case 'F'
            raster = [1 1 1; 1 0 0; 1 1 0; 1 0 0; 1 0 0];
        case 'G'
            raster = [0 1 1; 1 0 0; 1 0 1; 1 0 1; 0 1 1];
        case 'H'
            raster = [1 0 1; 1 0 1; 1 1 1; 1 0 1; 1 0 1];
        case 'I'
            raster = [1 1 1; 0 1 0; 0 1 0; 0 1 0; 1 1 1];
        case 'J'
            raster = [0 0 1; 0 0 1; 0 0 1; 1 0 1; 1 1 1];
        case 'K'
            raster = [1 0 1; 1 0 1; 1 1 0; 1 0 1; 1 0 1];
        case 'L'
            raster = [1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 1 1];
        case 'M'
            raster = [1 0 1; 1 1 1; 1 0 1; 1 0 1; 1 0 1];
        case 'N'
            raster = [1 1 1; 1 0 1; 1 0 1; 1 0 1; 1 0 1];
        case 'O'
            raster = [0 1 0; 1 0 1; 1 0 1; 1 0 1; 0 1 0];
        case 'P'
            raster = [1 1 0; 1 0 1; 1 1 0; 1 0 0; 1 0 0];
        case 'Q'
            raster = [0 1 0; 1 0 1; 1 0 1; 1 1 1; 0 1 1];
        case 'R'
            raster = [1 1 0; 1 0 1; 1 1 0; 1 0 1; 1 0 1];
        case 'S'
            raster = [0 1 1; 1 0 0; 1 1 1; 0 0 1; 1 1 0];
        case 'T'
            raster = [1 1 1; 0 1 0; 0 1 0; 0 1 0; 0 1 0];
        case 'U'
            raster = [1 0 1; 1 0 1; 1 0 1; 1 0 1; 1 1 1];
        case 'V'
            raster = [1 0 1; 1 0 1; 0 1 1; 0 1 1; 0 0 1];
        case 'W'
            raster = [1 0 1; 1 0 1; 1 0 1; 1 1 1; 1 0 1];
        case 'X'
            raster = [1 0 1; 1 0 1; 0 1 0; 1 0 1; 1 0 1];
        case 'Y'
            raster = [1 0 1; 1 0 1; 0 1 0; 0 1 0; 0 1 0];
        case 'Z'
            raster = [1 1 1; 0 0 1; 0 1 0; 1 0 0; 1 1 1];
        otherwise
            error('Input must be one of the 26 standard English letters');
    end
    
end