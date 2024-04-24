%
% Convert electrodes in native space to MNI through surface based
% registration (through Freesurfer sphere). Since it is surface based,
% it should ony be used for ECoG.
%
% mni_coords = ieeg_mni305ThroughFsSphere(elecmatrix, hemi, FSdir, FSsubjectsdir); % minimum inputs/outputs
% [mni_coords, vertIdxFsavg, minDists, surfUsed] = ieeg_mni305ThroughFsSphere(elecmatrix, hemi, FSdir, FSsubjectsdir, surfaceChoice, distLim); % maximum inputs/outputs
%
% Input:
%   elecmatrix: nElec x 3 in the original T1 space (before freesurfer)
%   hemi: hemisphere to load for each electrode (l/r)
%   FSdir: subjects freesurfer director
%   FSsubjectsdir: freesurfer subjects directory with fsaverage
%   surfaceChoice: (optional), char: 'pial' (default), 'white', or 'closest' to determine which surface rendering to snap electrodes to.
%                   'closest' choosest whichever surface is closest for each electrode
%   distLim: (optional), num, maximum distance to snap electrodes to subject brain. Skips electrodes too far from a surface. Default = inf (no liit)
%
% Returns:
%   mni_coords: nElec x 3 matrix in MNI305 space. Skipped electrodes (by distLim) are given [NaN, NaN, NaN].
%   vertIdxFsavg: nElec x 1, vertex index in fsaverage where the mni_coord came from. Useful to quickly get the matching inflated, sphere, or flat fsaverage coordinates
%   minDists: nElec x 1, distance from each electrode to the closest surface vertex (either white or pial depending on <surfaceChoice> arg. Returned also for skipped electrodes
%   surfUsed: nElec x 1 cell array of char, each entry is 'pial' or 'white', depending on which surface the electrode was snapped to
%
%
% Figure to check: get the mni sphere index in the mni pial
% load the freesurfer MNI surface
%   [mnipial_vert,mnipial_face] = read_surf(fullfile(FSsubjects,'fsaverage','surf','lh.pial'));
% 
%   figure
%   g.faces = mnipial_face+1;
%   g.vertices = mnipial_vert;
%   g = gifti(g);
%   tH = ieeg_RenderGifti(g);
%   ieeg_label(mni_coords)
%
% Dora Hermes, 2020
% Mutimodal Neuroimaging Lab
% Mayo Clinic
%
% Edited 2024/04/24 by HH
%   - added surfaceChoice, distLim optional arguments
%   - added vertIdxFsavg, minDists, surfUsed additional output arguments
%
function [mni_coords, vertIdxFsavg, minDists, surfUsed] = ieeg_mni305ThroughFsSphere(elecmatrix, hemi, FSdir, FSsubjectsdir, surfaceChoice, distLim)
    
    if isempty(distLim)
        distLim = inf;
    end
    
    if isempty(surfaceChoice)
        surfaceChoice = 'pial';
    end
    assert(any(strcmpi(surfaceChoice, {'pial', 'white', 'closest'})), 'Error: surfaceChoice must be "pial", "white", or "closest"');
    
    if isempty(FSsubjectsdir)
        disp('select freesurfer subjects directory with fsaverage dir')
        [FSsubjectsdir] = uigetdir(pwd,'select freesurfer directory');
    end
    
    if isempty(FSdir)
        disp('select this subjects freesurfer directory')
        [FSdir] = uigetdir(pwd,'select this subjects freesurfer directory');
    end
    
    % number of electrodes
    nElecs = size(elecmatrix,1);
    disp(['getting MNI305 coordinates for ' int2str(nElecs) ' electrodes'])
    
    % load mri orig header
    origName = fullfile(FSdir,'mri','orig.mgz');
    orig = MRIread(origName,'true');
    Norig = orig.vox2ras; 
    Torig = orig.tkrvox2ras;
    
    % electrodes to freesurfer space
    freeSurfer2T1 = inv(Norig*inv(Torig));
    elCoords = freeSurfer2T1*([elecmatrix'; ones(1, nElecs)]);
    elCoords = elCoords(1:3,:)';
    
    % subject surface vertices for pial and white
    Lsubpial_vert = read_surf(fullfile(FSdir,'surf','lh.pial'));
    Rsubpial_vert = read_surf(fullfile(FSdir,'surf','rh.pial'));
    Lsubwhite_vert = read_surf(fullfile(FSdir,'surf','lh.white'));
    Rsubwhite_vert = read_surf(fullfile(FSdir,'surf','rh.white'));
    
    % figure to check electrodes in freesurfer space
    % figure
    % g.faces = pial_face+1;
    % g.vertices = pial_vert;
    % tH = ieeg_RenderGifti(g);
    % ieeg_label(elCoords)
    % set(tH,'FaceAlpha',.5) % make transparent
    
    % subject sphere
    [Lsubsphere_vert] = read_surf(fullfile(FSdir,'surf','lh.sphere.reg'));
    [Rsubsphere_vert] = read_surf(fullfile(FSdir,'surf','rh.sphere.reg'));
    
    % mni305 sphere
    [Lmnisphere_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.sphere'));
    [Rmnisphere_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.sphere'));
    
    % mni305 pial
    [Lmnipial_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','lh.pial'));
    [Rmnipial_vert] = read_surf(fullfile(FSsubjectsdir,'fsaverage','surf','rh.pial'));
    
    % index for closest point to subject's pial
    s_pial_ind = nan(nElecs, 1);
    % subject sphere coords at these indices
    sphere_coords = nan(nElecs, 3);
    
    % also keep track of the minimum distance and the surface used
    minDists = nan(nElecs, 1);
    surfUsed = cell(nElecs, 1);
    
    % output mni coordinates and idx of best vertex on fsavg
    mni_coords = nan(nElecs,3);
    vertIdxFsavg = nan(nElecs, 1);
    
    for kk = 1:nElecs
        
        xyz = elCoords(kk,:);
        
        if isequal(upper(hemi{kk}),'L')
            
            % closest point between electrode and a surface (pial or white), and the associated minimum distance
            if strcmpi(surfaceChoice, 'pial')
                [min_dist, min_ind] = min(sqrt(sum((Lsubpial_vert - xyz).^2, 2)));
                surfUsed{kk} = 'pial';
            elseif strcmpi(surfaceChoice, 'white')
                [min_dist, min_ind] = min(sqrt(sum((Lsubwhite_vert - xyz).^2, 2)));
                surfUsed{kk} = 'white';
            elseif strcmpi(surfaceChoice, 'closest')
                [min_dist_pial, min_ind_pial] = min(sqrt(sum((Lsubpial_vert - xyz).^2, 2)));
                [min_dist_white, min_ind_white] = min(sqrt(sum((Lsubwhite_vert - xyz).^2, 2)));
                [min_dist, surfaceChoiceNum] = min([min_dist_pial, min_dist_white]); % choose the closest of the two, pial or white
                if surfaceChoiceNum == 1
                    min_ind = min_ind_pial;
                    surfUsed{kk} = 'pial';
                else
                    min_ind = min_ind_white;
                    surfUsed{kk} = 'white';
                end
            end
            minDists(kk) = min_dist;
    
            % did not find a vertex within distance limit, do not match to fsavg
            if min_dist > distLim, continue; end
    
            s_pial_ind(kk) = min_ind;
            
            % get the same index on subjects sphere
            sphere_coords(kk,:) = Lsubsphere_vert(min_ind,:);
            
            % closest point subjects sphere to mni sphere
            xyz_sphere = sphere_coords(kk,:);
            [~,min_ind_sphere] = min(sqrt(sum((Lmnisphere_vert-xyz_sphere).^2,2)));
            
            % get the  mni pial at the mni sphere index, also store the index
            mni_coords(kk,:) = Lmnipial_vert(min_ind_sphere,:);
            vertIdxFsavg(kk) = min_ind_sphere;
    
            
        elseif isequal(upper(hemi{kk}),'R')
    
            % closest point between electrode and a surface (pial or white), and the associated minimum distance
            if strcmpi(surfaceChoice, 'pial')
                [min_dist, min_ind] = min(sqrt(sum((Rsubpial_vert - xyz).^2, 2)));
                surfUsed{kk} = 'pial';
            elseif strcmpi(surfaceChoice, 'white')
                [min_dist, min_ind] = min(sqrt(sum((Rsubwhite_vert - xyz).^2, 2)));
                surfUsed{kk} = 'white';
            elseif strcmpi(surfaceChoice, 'closest')
                [min_dist_pial, min_ind_pial] = min(sqrt(sum((Rsubpial_vert - xyz).^2, 2)));
                [min_dist_white, min_ind_white] = min(sqrt(sum((Rsubwhite_vert - xyz).^2, 2)));
                [min_dist, surfaceChoiceNum] = min([min_dist_pial, min_dist_white]); % choose the closest of the two, pial or white
                if surfaceChoiceNum == 1
                    min_ind = min_ind_pial;
                    surfUsed{kk} = 'pial';
                else
                    min_ind = min_ind_white;
                    surfUsed{kk} = 'white';
                end
            end
            minDists(kk) = min_dist;
    
            % did not find a vertex within distance limit, do not match to fsavg
            if min_dist > distLim, continue; end
    
            s_pial_ind(kk) = min_ind;
            
            % get the same index on subjects sphere
            sphere_coords(kk,:) = Rsubsphere_vert(min_ind,:);
            
            % closest point subjects sphere to mni sphere
            xyz_sphere = sphere_coords(kk,:);
            [~,min_ind_sphere] = min(sqrt(sum((Rmnisphere_vert-xyz_sphere).^2,2)));
            
            % get the  mni pial at the mni sphere index, also store the index
            mni_coords(kk, :) = Rmnipial_vert(min_ind_sphere,:);    
            vertIdxFsavg(kk) = min_ind_sphere;
    
        end
        
    end

end

