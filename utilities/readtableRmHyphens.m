%% Create a table by reading from a file (via Matlab readtable())
% but removes hyphens from an input column (e.g. channel names) that shouldn't have hyphens (or has too many hyphens).
% Use for channels.tsv, electrodes.tsv, and events.tsv files.
%   E.g. use this to load a channels.tsv file so that channel names like 'RC-ALT2' become 'RCALT2'.
%   E.g. alternatively, use this to load an events.tsv file so that stim pair names like 'RC-ALT2-RC-ALT3' become 'RCALT2-RCALT3'.
%
%   tbl = readtableRmHyphens(path);
%   tbl = readtableRmHyphens(path, colName, numHyphens);
%   tbl = readtableRmHyphens(path, colName, numHyphens, varargin);
%       path =          char, full or relative path to table file containing file name and extension
%       colName =       char (optional), name of the column to remove hyphens from. Default = 'name' (works for channels.tsv
%                           and electrodes.tsv files).
%       numHyphens =    0 or 1 (optional), number of hyphens that SHOULD exist for each item in the column of interest.
%                           E.g. 0 for channels.tsv like 'RC1' and 1 for events.tsv like 'RC1-RC2'. Default = 0.
%                           If 1 is selected, it is assumed that the number of unwanted hyphens is balanced on both
%                           sides of the central hyphen. E.g. 'RC-ALT1-RC-ALT2' -> 'RCALT1-RCALT2', but 'RC-ALT1-RCALT2'
%                           will throw an error.
%       varargin =      variable input arguments for Matlab readtable() function (optional). Corresponds to either
%                           importoptions opts or Name, Value pairs. See Matlab readtable documentation for more details.
%                           Default = {'FileType', text', 'Delimiter', '\t'} which corresponds to tab-delineated (.tsv) files
%
% HH 2021
%
function [tbl, origCol] = readtableRmHyphens(path, colName, numHyphens, varargin)
    
    if nargin < 4
        tbl = readtable(path, 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a'); % default options for reading tsv files
    else
        tbl = readtable(path, varargin{:});
    end
    
    if nargin < 3 || isempty(numHyphens), numHyphens = 0; end
    if nargin < 2 || isempty(colName), colName = 'name'; end % which column to remote hyphens in
    assert(ismember(numHyphens, [0, 1]), 'Function is unequipped to handling columns with more than 1 deliberate hyphen');
    
    origCol = tbl.(colName);
    col = origCol;
    for cc = 1:length(col)
        
        if any(strcmp(col(cc), {'n/a', 'NaN'})), continue; end % ignore rows with NaN or n/a
        
        eles = split(col(cc), '-');
        
        switch numHyphens
            case 0 % Remove all hyphens
                col(cc) = join(eles, '');
                
            case 1
                assert(mod(length(eles), 2) == 0, "Must have even number of '-'-split elements if 1 hyphen is required");
                col(cc) = join([join(eles(1:length(eles)/2), ''), join(eles(length(eles)/2+1:end), '')], '-');
        end
        
    end
    
    tbl.(colName) = col;
    
end