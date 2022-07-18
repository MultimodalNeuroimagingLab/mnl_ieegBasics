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
    
    if size(pairs, 2) == 1 % assume given as single column like 'ch1-ch2'
        try
            pairs = split(pairs, '-', 2);
        catch
            error('If pairs is given as a single column, each row must be in the form ''ch1-ch2'' (hyphen-joined)');
        end
    end

    xyzs_pair = NaN(size(pairs, 1), 3);
    xyzs = [electrodes.x, electrodes.y, electrodes.z]; % ordered coordinates for the electrodes
    
    for i=1:size(pairs, 1)
        loc1 = xyzs(strcmp(pairs{i, 1}, electrodes.name), :);
        loc2 = xyzs(strcmp(pairs{i, 2}, electrodes.name), :);
        xyzs_pair(i, :) = mean([loc1; loc2], 1);
    end
end