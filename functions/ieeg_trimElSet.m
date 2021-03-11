%% This function recursively removes indices from a set until all pairwise distances are greater than a threshold
% Used to trim a set of electrodes down to a subset in which all electrodes are at least <minDist> far from each other
% helper function to ieeg_getXyzMni.m
%
%   setOut = trimElSet(setIn, dists, minDist);
%       setIn =         1xm double, set of indices for elements that form <dists>, where each element corresponds to
%                           one row/col in <dists>. m <= n, and the maximum value in <setIn> <= n.
%       dists =         nxn double, square matrix, pairwise distances between all items (e.g. electrodes) indexed by <setIn>
%       minDist =       1x1 double, minimum pairwise distance allowed between all pairs in the trimmed output set.
%
%   returns
%       setOut =        1xp double, subset of indices from setIn such that dists(setOut, setOut) is at least <minDist>. p <= m.
%
%   HH 2021
%
function setOut = ieeg_trimElSet(setIn, dists, minDist)

    dists(1:size(dists,1)+1:end) = inf; % make sure diagonal (self-self) distances are ignored
    subDists = dists(setIn, setIn);
    
    if min(subDists, [], 'all') >= minDist
        setOut = setIn; % BASE CASE: all pairwise distances are far enough
    else
        [~, toRemove] = min(min(subDists)); % remove (1st) index with closest distance to another
        setIn(toRemove) = [];
        setOut = ieeg_trimElSet(setIn, dists, minDist);
    end
    
end