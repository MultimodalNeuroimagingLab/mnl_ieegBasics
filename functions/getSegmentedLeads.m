%%  This function takes an electrodes table as input and outputs which leads are segmented in 5 and in 6
%   Assumes that inter-electrode distances greater than 2x the median distance for that lead correspond to segmented gaps
%   For leads that are long enough to have 2 gaps (i.e., 11 for 5-segmented, 13 for 6-segmented), this function will check to make sure both are present.
%   Otherwise, a warning is returned and the lead is not considered segmented
%
%   USAGE:
%   [seg5, seg6] = getSegmentedLeads(electrodes)
%       electrodes =    table. Must have columns "name", "x", "y", "z". There cannot be NaN positions before more valid contacts. E.g., can't have NaN in RA11 if RA12 has positions.
%
%   RETURNS:
%       seg5 =          cell array of characters. Names of leads that are segmented by 5 contacts
%       seg6 =          cell array of characters. Names of leads that are segmented by 6 contacts
%
%   HH 2023
%
function [seg5, seg6] = getSegmentedLeads(electrodes)
    
    isDigit = @(x) x > 47 & x < 58; % returns true for char array elements that are digits (0 - 9)
    xyz = [electrodes.x, electrodes.y, electrodes.z];
    names = upper(strip(electrodes.name));
    
    % stores leads that are segmented in 5 and 6
    seg5 = {};
    seg6 = {};
    
    % separate electrode names into char and numerical positions
    elecPos = cellfun(@(ch) str2double(ch(isDigit(ch))), names);
    elecChar = cellfun(@(ch) ch(~isDigit(ch)), names, 'UniformOutput', false);
    leads = unique(elecChar, 'stable');
    
    for ll = 1:length(leads)
        
        % all positions at current lead, sorted
        thislead = strcmp(elecChar, leads{ll}); % log indices of elecs in this lead
        [leadPos, ix] = sort(elecPos(thislead)); % sort electrode by ascending lead positions
        
        % get the electrode xyz positions, in same order
        xyzCurr = xyz(thislead, :);
        xyzCurr = xyzCurr(ix, :);
        
        % remove nan electrodes, if any
        leadPos(any(isnan(xyzCurr), 2)) = [];
        xyzCurr(any(isnan(xyzCurr), 2), :) = [];
        
        assert(all(diff(leadPos) == 1), 'Error: Non-contiguous electrode positions in lead %s', leads{ll});
        
        if length(leadPos) < 6, continue; end % no segments can be detected for 5 or fewer contacts
        
        dists = vecnorm(diff(xyzCurr), 2, 2);
        distMed = median(dists); % median contact-to-contact distance
        
        [distSorted, distOrder] = sort(dists, 'descend');
        leadPosByDist = leadPos(distOrder);
        
        if distSorted(1) < 2*distMed % cannot be a segmented lead because distance not great enough
            continue
        end
        
        if leadPosByDist(1) == 6 || leadPosByDist(1) == 12 % largest distance is 6-7 jump or 12-13 jump
            if length(leadPos) >= 13 % enough contacts to also perform double check on second disatnce
                if distSorted(2) > 2*distMed && (leadPosByDist(2) == 6 || leadPosByDist(2) == 12) % also large jump in the second distance, which is also 6 or 12
                    seg6 = [seg6; leads{ll}];
                else % 2nd large jump does not belong to 6-segmented pattern
                    warning('2nd largest jump inconsistent with first largest jump as 6-segmented lead in %s', leads{ll});
                end
            else
                seg6 = [seg6; leads{ll}]; % assume short contact is 6-segmented
            end

        elseif leadPosByDist(1) == 5 || leadPosByDist(1) == 10 % 5-segmented
            if length(leadPos) >= 11 % enough contacts to also perform double check on second disatnce
                if distSorted(2) > 2*distMed && (leadPosByDist(2) == 5 || leadPosByDist(2) == 10) % also large jump in the second distance, which is also 6 or 12
                    seg5 = [seg5; leads{ll}];
                else
                    warning('2nd largest jump inconsistent with first largest jump as 5-segmented lead in %s', leads{ll});
                end
            else
                seg5 = [seg5; leads{ll}]; % assume short contact is 5-segmented
            end
        end
        
    end
    
end