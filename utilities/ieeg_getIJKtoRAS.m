%% This function loads a nifti and writes its ijk to RAS (voxel to position) matrix out as a text file.
%
%   USAGE:
%       ieeg_getIJKtoRAS(niftiPath);
%       ijk2ras = ieeg_getIJKtoRAS(niftiPath, outputFname);
%           niftiPath =         char, path to input nifti file
%           outputFname =       (optional) char, name of file to save. If not given, output file will have same path as niftiPath
%                                   but end in _ijk2ras.txt. Always saves as a .txt text file.
%                                   YAEL requires this file to be named 'CT_IJK_to_MR_RAS.txt'
%
%   RETURNS:
%       ijk2ras =           4x4 num, the affine matrix representing the IJK to RAS coordinate transformation, which is saved to file
%
%   DEPENDENCIES:
%       Either: spm_vol from SPM12 or niftiRead from Vistasoft
%
%   HH 2023/09
%
function ijk2ras = ieeg_getIJKtoRAS(niftiPath, outputFname)

    % Get dir of niftiPath to make output IJKtoRAS output path
    listing = dir(niftiPath);
    if nargin < 2 % name it the same way as input nifti
        outputFname = sprintf('%s_ijk2ras.txt', listing.name(1:end-4));
    elseif ~endsWith(outputFname, '.txt') % add .txt if it isn't the ending
        outputFname = sprintf('%s.txt', outputFname);
    end
    outpath = fullfile(listing.folder, outputFname);

    % try using niftiRead first (VistaSoft)
    try
        nii = niftiRead(niftiPath);
        ijk2ras = nii.sto_xyz;
        writematrix(ijk2ras, outpath, 'Delimiter', ' ');
        fprintf('Saved ijk_to_ras matrix to ''%s'' using niftiRead\n', outpath);
        return

    catch ME

        % Raise error if it isn't because niftiRead is undefined
        if ~strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            rethrow(ME);
        end

    end

    % niftiRead failed, try using spm_vol from SPM12 instead
    try
        nii = spm_vol(niftiPath);
        ijk2ras = nii.mat; % vox -> pos transformation
        writematrix(ijk2ras, outpath, 'Delimiter', ' ');
        fprintf('Saved ijk_to_ras matrix to ''%s'' using spm_vol\n', outpath);

    catch ME

        % Error that isn't undefined function
        if ~strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            rethrow(ME);
        end
        
        % Enhance error message if spm_vol and niftiRead both were undefined
        ME2 = MException('MATLAB:UndefinedFunction', 'Both niftiRead (from Vistasoft) and spm_vol (from SPM12) functions are undefined. At least one must exist in your path.');
        throw(ME2);

    end

end