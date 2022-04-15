%% Saves a sorted electrodes_sorted.tsv table by re-ordering rows in the input electrodes.tsv file according to the row 
%   order of channel names in the input channels.tsv file.
%
%   This is relevant because the channels.tsv channel names match the ephys data channels, but the rows in
%   electrodes.tsv are often arbitrarily sorted. T.f. the electrodes.tsv rows must be sorted in order to plot electrode
%   positions that correspond correctly to row indices of ephys data.
%   
%   sortElectrodes(elecPath, channelPath);
%   elecsOut = sortElectrodes(elecPath, channelPath, saveFile);
%       elecPath =      str, path to electrodes.tsv file
%       channelPath =   str, path to channels.tsv file
%       saveFile =      bool (optional). Whether to save the sorted electrodes to file (default = true)
%
%   Returns:
%       elecsOut =      table, rows of electrodes.tsv sorted in order of channel names in channels.tsv.
%                           Channel names not found in electrodes.tsv (i.e. ends of leads without coordinates) will be
%                           NaN rows in elecsOut, in order to align overall table structure with channels
%
%   HH 2021
%
function elecsOut = sortElectrodes(elecsPath, channelsPath, saveFile)

    if nargin < 3, saveFile = true; end
    
    channels = readtable(channelsPath, 'FileType', 'text', 'Delimiter', '\t'); % keep hyphens for filename
    elecs = readtable(elecsPath, 'FileType', 'text', 'Delimiter', '\t');
    
    % RM hyphens
    
    elecNames = erase(strip(elecs.name), '-'); % delete hyphens from elecs and channels to ensure they match
    channelNames = erase(strip(channels.name), '-');
    
    assert(length(unique(elecNames)) == length(elecNames), 'Repeated names in electrodes file');
    assert(length(unique(channelNames)) == length(channelNames), 'Repeated names in channels file');

    [~, locb] = ismember(channelNames, elecNames);
    
    % varTypes in input cell, removing 
    varTypes = cellfun(@class, table2cell(elecs(1, :)), 'UniformOutput', false);
    varTypes(strcmp('char', varTypes)) = {'string'}; % so that Matlab doesn't scream at me for preallocating with 'char'
    varTypes{1} = 'name'; % ignore first col for indexing purposes
    
    elecsOut = table('Size', [length(channelNames), length(elecs.Properties.VariableNames)], ...
                            'VariableNames', elecs.Properties.VariableNames, ...
                            'VariableTypes', repmat({'string'}, [1, length(elecs.Properties.VariableNames)]));
    elecsOut.name = channels.name; % save with original channel names (whether or not they had hyphens)
    elecsOut(logical(locb), ~strcmp(varTypes, 'double')) = elecs(locb(locb > 0), ~strcmp(varTypes, 'double')); % directly put in non-double rows
    
    % For double variable types, save explicitly with 8 digits of precision
    numValues = table2cell(elecs(locb(locb > 0), strcmp(varTypes, 'double')));
    numValuesStr = cellfun(@(x) num2str(x, 8), numValues, 'UniformOutput', false);
    numValuesStr(strcmp(numValuesStr, 'NaN')) = {'n/a'}; % convert NaNs to 'n/a'
    elecsOut(logical(locb), strcmp(varTypes, 'double')) = numValuesStr;
    
    %nanRow = cell(1, length(varTypes)); % what to put into rows without input electrode info
    %nanRow(strcmpi(varTypes, 'double')) = {nan};
    %nanRow(~strcmpi(varTypes, 'double')) = {'n/a'};
    nanRow = repmat({'n/a'}, [1, length(elecs.Properties.VariableNames)]); % row of all 'n/a's
    elecsOut(~logical(locb), 2:end) = repmat(nanRow(2:end), sum(~logical(locb)), 1);
    
    % convoluted way of replacing <"missing"> values (previously NaN) with 'n/a' because fillmissing does not work with strings
%     C = table2cell(elecsOut);
%     C(ismissing(elecsOut)) = {'n/a'};
%     elecsOut = cell2table(C, 'VariableNames', elecsOut.Properties.VariableNames);
        
    [saveDir, name] = fileparts(elecsPath);
    if saveFile
        outPath = fullfile(saveDir, sprintf('%s_sorted.tsv', name));
        if exist(outPath, 'file'), warning('Overwriting existing electrodes_sorted.tsv'); end
        writetable(elecsOut, outPath, 'FileType', 'text', 'Delimiter', '\t');
    end
    
end
 