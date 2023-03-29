function [out] = ieeg_eldist2pial2white(loc_info,gL_pial,gR_pial,gL_white,gR_white, distLim)
%
% get relative distances from electrode to white, negative values indicate
% that the electrode is in white matter
%
% 1: Test whether pial or gray/white border is closest to electrode 
% 2: After we have the closest point in pial or white, we also get the corresponding coordinate on
% the white/pial using the index.
% 3: We then calculate distances
% 4: We calculate the relative distance to white, ensuring that electrodes
%    below the white are negative. 
%
% Input
% [out] = ieeg_eldist2pial2white(loc_info,gL_pial,gR_pial,gL_white,gR_white, distLim)
%   loc_info: bids_electrodes_tsv_file needs columns with x,y,z and hemisphere.
%       Can be read with readtable(bids_electrodes_tsv_file,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
%   gL_pial: left pial surface in subject/electrode space, can load with gifti('pial.L.surf.gii'); 
%   gR_pial: right pial surface
%   gL_white: left gray/white matter surface
%   gR_white: right gray/white matter surface
%
% Output
%     out: will have several fields with information for each electrode
%           out.name: electrode name from electrodes.tsv
%       out.xyz_pial: coordinates of vertex on pial closest to electrode
%      out.xyz_white: coordinates of white on pial closest to electrode
%      out.dist_pial: distance from electrode to pial
%     out.dist_white: distance from electrode to white
% out.dist_pialwhite: distance from pial to white
%      out.surfIndex: index into freesurfer vertex
%       out.rel_dist: relative distance to gray/white matter border,
%                     negative values are in the white matter
%
% DH 2023, Multimodal Neuroimaging Lab

elecmatrix = [loc_info.x loc_info.y loc_info.z];
out.name = loc_info.name;

%% %%%% snap gray matter electrodes to surface, get index on surface, and get point on inflated brain
out.xyz_pial = NaN(size(elecmatrix));
out.xyz_white = NaN(size(elecmatrix));
out.dist_pial = NaN(size(elecmatrix,1),1);
out.dist_white = NaN(size(elecmatrix,1),1);
out.dist_pialwhite = NaN(size(elecmatrix,1),1);
out.surfIndex = NaN(size(elecmatrix,1),1);

if isempty(distLim)
    distLim = 6;
end

for kk = 1:height(loc_info)
    if ~isnan(elecmatrix(kk,1)) % we have a position

        % get coordinates
        xyz = elecmatrix(kk,:);

        % check for left or right
        if strcmp(loc_info.hemisphere{kk},'L')
            g_pial = gL_pial;
            g_white = gL_white;
        elseif strcmp(loc_info.hemisphere{kk},'R')
            g_pial = gR_pial;
            g_white = gR_white;
        end

        % find closest point on pial surface
        dist_xyz2pial = sqrt(sum((g_pial.vertices-xyz).^2,2));
        [minDist_pial,ind_minDist_pial] = min(dist_xyz2pial);

        % find closest point on gray/white surface
        dist_xyz2white = sqrt(sum((g_white.vertices-xyz).^2,2));
        [minDist_white,ind_minDist_white] = min(dist_xyz2white);

        if minDist_pial < distLim || minDist_white  < distLim % only find points less than 6 mm away
            if minDist_pial < minDist_white  % take the index from the pial
                out.surfIndex(kk) = ind_minDist_pial;
            else % take the index from the white
                out.surfIndex(kk) = ind_minDist_white;
            end
            out.xyz_pial(kk,:) = g_pial.vertices(out.surfIndex(kk),:);
            out.xyz_white(kk,:) = g_white.vertices(out.surfIndex(kk),:);
            out.dist_pial(kk) = dist_xyz2pial(out.surfIndex(kk));
            out.dist_white(kk) = dist_xyz2white(out.surfIndex(kk));
            out.dist_pialwhite(kk) = sqrt(sum( (g_pial.vertices(out.surfIndex(kk),:) - g_white.vertices(out.surfIndex(kk),:)).^2 ));
        end

    end
end

%% now we get relative position of the electrode to gray/white border

out.rel_dist = NaN(size(out.dist_white)); % relative position of the electrode to gray/white border

for kk = 1:size(elecmatrix,1)
    if ~isnan(out.dist_pial(kk))
        if out.dist_pial(kk) <= out.dist_pialwhite(kk) && out.dist_white(kk) <= out.dist_pialwhite(kk) 
            % electrode between pial and white
            out.rel_dist(kk) = out.dist_white(kk);
        elseif out.dist_pial(kk) > out.dist_pialwhite(kk) && out.dist_pial(kk) >= out.dist_white(kk) 
            % electrode below white
            out.rel_dist(kk) = -out.dist_white(kk);
        elseif out.dist_white(kk) > out.dist_pialwhite(kk) && out.dist_white(kk) >= out.dist_pial(kk) 
            % electrode above pial 
            out.rel_dist(kk) = out.dist_white(kk);
        else
            disp('relative distance not captured - check whats going on')
        end

    end
end
