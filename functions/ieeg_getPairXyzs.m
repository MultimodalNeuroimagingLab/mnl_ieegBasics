%
%   Returns the linearly-interpolated middle XYZ coordinate between pairs of electrodes
%   
%   xyzs_pair = ieeg_getPairXyzs(pairs, electrodes);
%   
%       pairs =         nx2 cell array of electrode names: n pairs, where column 1 = elec1, column 2 = elec2
%       electrodes =    table with columns 'name', 'x', 'y', 'z', with each row corresponding to the XYZ coordinate of a
%                           given electrode. All electrode names in <pairs> must be found in <electrodes>
%
%   Returns:
%       xyzs_pair =     nx3 array of XYZ pair coordinates. Rows correspond to the rows in input <pairs>.
%
%   HH 2020
%

function xyzs_pair = ieeg_getPairXyzs(pairs, electrodes)
    
    assert(size(pairs, 2) == 2, 'pairs must be an nx2 cell array: col1=elec1, col2=elec2');

    xyzs_pair = NaN(size(pairs, 1), 3);
    xyzs = [electrodes.x, electrodes.y, electrodes.z]; % ordered coordinates for the electrodes
    
    for i=1:size(pairs, 1)
        loc1 = xyzs(strcmp(pairs{i, 1}, electrodes.name), :);
        loc2 = xyzs(strcmp(pairs{i, 2}, electrodes.name), :);
        xyzs_pair(i, :) = mean([loc1; loc2], 1);
    end
end