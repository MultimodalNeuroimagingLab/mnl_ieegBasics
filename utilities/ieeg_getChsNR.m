%% This function returns channels that are non-recording.
%
%   chsNR = ieeg_getChsNR(channels, electrodes);
%   [chsNR, chNamesNR] = ieeg_getChsNR(channels, electrodes);
%   ieeg_getChsNR(channels, electrodes, true); % overwrites existing channels file with a copy that has NR channels
%       
%       electrodes      = table or char. Electrodes: either the table or the path to the table.
%       channels        = table or char. Channels: either the channels or the path to the channels.tsv
%       overwrite       = (optional) logical. Default = false. If true, the channels file is overwritten with a status column that corresponds to NR channels
%                               This only applies if channels is given as path (char input)
%
%   Returns:
%       chsNR           = nx1 integers. Channel indices that are non-recording.
%       chNamesNR       = nx1 char cell array. Channel names that are non-recording
%
%   Dependencies:
%       sortElectrodes.m
%
%   Harvey Huang, Multimodal Neuroimaging Lab, Mayo Clinic, 2022
%
%
function [chsNR, chNamesNR] = ieeg_getChsNR(channels, electrodes, overwrite)
    
    if nargin < 3, overwrite = false; end

    elecsSorted = ieeg_sortElectrodes(electrodes, channels, false);
    
    chsNR = find(isnan(elecsSorted.x));
    chNamesNR = elecsSorted.name(chsNR);
    
    % overwrite channels file with NR channels
    if overwrite && (ischar(channels) || isstring(channels))
        channelsPath = channels;
        channels = readtable(channels, 'Filetype', 'text', 'Delimiter', '\t');
        statusCurr = channels.status;
        
        %status = repmat({'good'}, height(channels), 1);
        % keep existing status column
        status = statusCurr;
        status(chsNR) = {'bad'};
        
        %status(ismember(channels.type, {'EKG', 'ECG', 'MISC'})) = {'good'}; % reset the EKG and MISC channels to good
        status(ismember(channels.type, {'EKG', 'ECG', 'MISC'})) = statusCurr(ismember(channels.type, {'EKG', 'ECG', 'MISC'})); % reset the EKG and MISC channels to what they were
        
        % which channels are being changed
        statusDiff = find(~strcmp(statusCurr, status));

        % check with user before overwriting if there are currently bad channels present and the status column is different than what will be saved
        if all(strcmp(status, statusCurr))
            fprintf('No new non-recording channel marked as bad. Nothing to change or save\n');
            return
            
        elseif any(strcmp(statusCurr, 'bad')) && any(strcmp(status, 'bad'))
            answer = questdlg(sprintf('There are already bad channels marked in this file. You would ADD %d bad channel(s)\nOVERWRITE?', length(statusDiff)), ...
                              'Overwrite channels', 'Yes, OVERWRITE', 'No, cancel', 'No, cancel');
            if ~strcmp(answer, 'Yes, OVERWRITE'), disp('Exited'); return; end
        end
        
        channels.status = status;
        
        fprintf('Channels file overwritten at %s\n', channelsPath);
        writetable(channels, channelsPath, 'Filetype', 'text', 'Delimiter', '\t')
    end
    
end