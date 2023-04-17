%% This function performs bipolar referencing for data recorded from SEEG channels, based on input channel names.
%
%   Bipolar signals are calculated as the signal difference between two adjacent electrode contacts on the same SEEG lead.
%   To determine which channels should be subtracted from each other, this function requires the channel names to be in
%   the format 'XXXnn', where 'XXX' = chars of any length that represent the unique name of an SEEG lead, and 'nn' is
%   the electrode position on that lead, starting from 1. E.g. 'LA12' is the 12th electrode along the 'LA' lead.
%
%   Double check that ONLY SEEG data and channels are used as input
%
%   The bipolar signals are then calculated starting from the 1st electrode in a lead upward, in the order of
%   subtracting the smaller number contact from the larger number contact. E.g. LA1 - LA2, LA2 - LA3, etc.
%
%   The user can specify a numerical list of bad channel indices. An output list of bad bipolar channels will be returned.
%   A bad bipolar channel is defined here as any bipolar pair that contains an input bad channel. The user can also specify, by channel name,
%   which leads are segmented by 5s or 6s, and the cross-segment pairs will not be returned in the bipolar output.
%
%   Usage:
%   [dataOut, bipolarChans] = ieeg_bipolarSEEG(dataIn, channelNames, []);
%   [dataOut, bipolarChans, badChans] = ieeg_bipolarSEEG(dataIn, channelNames, badChans);
%   [dataOut, bipolarChans, badChans] = ieeg_bipolarSEEG(dataIn, channelNames, badChans, seg5, seg6);
%       dataIn =            txn double. Signal data: rows are samples, columns are channels.
%       channelNames =      nx1 cell{char array}. Channel names must be in the same order as columns of <dataIn>. See
%                               the first paragraph of the documentation above for how channel names should be formatted.
%       badChans =          mx1 num, m <= n. Indices of bad channels to not use in bipolar derivation. bipolar pairs containing
%                               bad channels will be listed in <badChansOut>.
%       seg5 =              (optional) char or cell array of char. List of leads that are segmented in groups of 5. E.g., 'LA'. Bipolar pairs that cross segments
%                               will not be included in the output
%       seg6 =              (optional) char or cell array of char. List of leads that are segmented in groups of 6.
%
%   Returns:
%       dataOut =           txb double. Bipolar-referenced signal data. Rows are samples, columns are bipolar re-referenced
%                               channels.
%       bipolarNames =      bx1 cell. channel pair names corresponding to columns in <dataOut>.
%                               The names take the form 'ch1-ch2' to represent which pair of channels
%                               are used for the derivation. E.g. 'LA1-LA2' means data from 'LA1' minus data from 'LA2'
%       badChansOut =       px1 num, p <= b. Indices of bipolar rereferenced channels that are bad because they contain channels in <badChans>
%
%   Example:
%       
%       % create 48 fake channels: RA1 ... RA12, RB1 ... RB18, RC1 ... RC18
%       channelNames = [arrayfun(@(x) sprintf('RA%d', x), 1:12, 'UniformOutput', false)';
%                       arrayfun(@(x) sprintf('RB%d', x), 1:18, 'UniformOutput', false)'
%                       arrayfun(@(x) sprintf('RC%d', x), 1:18, 'UniformOutput', false)'];
%       
%       % some fake data to match
%       dataIn = rand(4800, 48);
%
%       % RA4, RB8, RB17, RB18, RC16-18 are all bad channels
%       badChans = [4, 20, 29, 30, 46:48];
%       disp('Bad channels:');
%       disp(channelNames(badChans));
%
%       % Calculate bipolar re-referenced data. RB and RC are 6-segmented leads, while RA is not a segmented lead
%       [dataOut, bipolarChans, badChansOut] = ieeg_bipolarSEEG(dataIn, channelNames, badChans, [], {'RB', 'RC'});
%
%
%   HH 2023
%
function [dataOut, bipolarNames, badChansOut] = ieeg_bipolarSEEG(dataIn, channelNames, badChans, seg5, seg6, verbose)

    if nargin < 6, verbose = true; end
    if nargin < 5 || isempty(seg6), seg6 = {}; end % no 6-segmented leads
    if nargin < 4 || isempty(seg5), seg5 = {}; end % no 5-segmented leads
    if nargin < 3 || isempty(badChans), badChans = []; end % no bad channels
    
    % Convert to cell arrays if single lead given
    if ischar(seg5), seg5 = {seg5}; end
    if ischar(seg6), seg6 = {seg6}; end
    
    assert(isempty(intersect(seg5, seg6)), 'One or more lead(s) are input as segmented by 5s AND by 6s');
        
    isDigit = @(x) x > 47 & x < 58; % returns true for char array elements that are digits (0 - 9)
    
    assert(length(channelNames) == size(dataIn, 2), 'Number of columns in dataIn must match number of channel names');
    
    emptyChs = cellfun(@isempty, channelNames);
    channelNames(emptyChs) = []; dataIn(:, emptyChs) = [];
    channelNames = strip(upper(channelNames)); % some cleaning
    
    leads = unique(cellfun(@(ch) ch(~isDigit(ch)), channelNames, 'UniformOutput', false), 'stable'); % unique lead names
    leads(cellfun(@isempty, leads)) = []; % discard any channel names that don't contain letters
    
    if any(~ismember(seg5, leads)) || any(~ismember(seg6, leads)), error('Check segmented leads input: some are not valid lead names'); end
    
    channelPos = cellfun(@(ch) str2double(ch(isDigit(ch))), channelNames); % electrode position along lead for each channel
    channelChar = cellfun(@(ch) ch(~isDigit(ch)), channelNames, 'UniformOutput', false); % just the character name of each lead
    
    minuend = []; subtrahend = []; % which indices to subtract from which
    badChansOut = [];
    if verbose, fprintf('Input leads:\n'); end
    for ll = 1:length(leads)
        
        % position of contacts (numbers) along the current lead.
        [leadPos, ix] = sort(channelPos(strcmp(channelChar, leads{ll})));
        % Indices into input channels matching leadPos
        leadInds = find(strcmp(channelChar, leads{ll})); leadInds = leadInds(ix);
        
        % Quick checks on position
        assert(length(unique(leadPos)) == length(leadPos), 'Error: Repeated lead position in %s', leads{ll});
        assert(all(diff(leadPos) == 1), 'Error: Non-consecutive electrode position in %s', leads{ll});
        
        minuendLead = leadInds(1:end-1);
        subtrahendLead = leadInds(2:end);
        
        % remove cross-segment bipolar pairs
        if ~isempty(seg5) && any(strcmpi(leads{ll}, seg5))
            if mod(leadPos(end), 5), warning('%s lead is marked as segmented in 5s, but it is length-%d (non-multiple)\n', leads{ll}, leadPos(end)); end
            minuendLead(~mod(leadPos(1:end-1), 5)) = [];
            subtrahendLead(~mod(leadPos(1:end-1), 5)) = []; % Remove positions that would be 5-6, 10-11, etc.
            if verbose, fprintf('Lead %s, contacts 1 ... %d, segmented in 5s\n', leads{ll}, max(leadPos)); end
        elseif ~isempty(seg6) && any(strcmpi(leads{ll}, seg6))
            if mod(leadPos(end), 6), warning('%s lead is marked as segmented in 6s, but it is length-%d (non-multiple)\n', leads{ll}, leadPos(end)); end
            minuendLead(~mod(leadPos(1:end-1), 6)) = [];
            subtrahendLead(~mod(leadPos(1:end-1), 6)) = []; % 6-7, 12-13, etc 
            if verbose, fprintf('Lead %s, contacts 1 ... %d, segmented in 6s\n', leads{ll}, max(leadPos)); end
        else
            if verbose, fprintf('Lead %s, contacts 1 ... %d\n', leads{ll}, max(leadPos)); end
        end
        
        % add which indices to subtract from which
        minuend = [minuend; minuendLead];
        subtrahend = [subtrahend; subtrahendLead];
        
    end
    
    dataOut = dataIn(:, minuend) - dataIn(:, subtrahend);
    badChansOut = find(ismember(minuend, badChans) | ismember(subtrahend, badChans)); % determine which channels are bad
    bipolarNames = join([channelNames(minuend), channelNames(subtrahend)], '-');
        
    if verbose, fprintf('Constructed %d bipolar channels, %d of which are bad.\n', length(minuend), length(badChansOut)); end
    

end