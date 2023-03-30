%% This function performs laplacian referencing for data recorded from SEEG channels, based on input channel names.
%
%   Each channel is re-referenced using the mean of its adjacent neighbors.
%   Does not normalize by electrode distance, see: https://www.sciencedirect.com/science/article/pii/S0165027001003302#BIB46
%   Electrodes on the ends of leads are referenced bipolar. E.g., Schalk Neuroimage 2018
%       Can consider Vaknin correction instead: https://doi.org/10.1016/0165-0270(88)90056-8 
%
%   This function requires the channel names to be in
%   the format 'XXXnn', where 'XXX' = chars of any length that represent the unique name of an SEEG lead, and 'nn' is
%   the electrode position on that lead, starting from 1. E.g. 'LA12' is the 12th electrode along the 'LA' lead.
%
%   Double check that ONLY SEEG data and channels are used as input
%
%   The user can specify a numerical list of bad channel indices. An output list of any rereferenced channels using a bad channel will be returned.
%   The user can also specify, by channel name, which leads are segmented by 5s or 6s, and electrodes along the segments will be bipolar referenced (not cross the gap)
%
%   USAGE:
%   [dataOut, badChansLap] = ieeg_lapSEEG(dataIn, channelNames);
%   [dataOut, badChansLap] = ieeg_lapSEEG(dataIn, channelNames, badChans, seg5, seg6, inclEdges, verbose);
%       dataIn =            txn double. Signal data: rows are samples, columns are channels.
%       channelNames =      nx1 cell{char array}. Channel names must be in the same order as columns of <dataIn>. See
%                               the first paragraph of the documentation above for how channel names should be formatted.
%       badChans =          mx1 num, m <= n. Indices of bad channels to not use in bipolar derivation. Rereferenced channels containing
%                               bad channels will be listed in <badChansLap>.
%       seg5 =              (optional) char or cell array of char. List of leads that are segmented in groups of 5. E.g., 'LA'. Channels at the ends of segments
%                               are bipolar-referenced if inclEdges is true, and not returned if inclEdges is false
%       seg6 =              (optional) char or cell array of char. List of leads that are segmented in groups of 6.
%       inclEdges =         (optional) logical. If true, lead ends and segment ends (within lead) are included in output but BIPOLAR referenced. If false, they are not included in output. Default = false
%
%   Returns:
%       dataOut =           txn double. Laplacian-rereferenced signals matching input data size. If inclEdges is false, channels on segment and lead edges are set as NaN.
%       badChansLap =       px1 num, p <= b. Indices of rereferenced channels that are bad because they contain channels in <badChans> input
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
%       % Calculate laplacian re-referenced data. RB and RC are 6-segmented leads, while RA is not a segmented lead
%       [dataOut, badChansLap] = ieeg_lapSEEG(dataIn, channelNames, badChans, [], {'RB', 'RC'}, false);
%
%   HH 2023
%
function [dataOut, badChansLap] = ieeg_lapSEEG(dataIn, channelNames, badChans, seg5, seg6, inclEdges, verbose)

    if nargin < 7, verbose = true; end
    if nargin < 6 || isempty(inclEdges), inclEdges = false; end % whether the lead edges should be included, bipolar referenced
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
    
    % Each ch is referenced as ch - (ref1 + ref2)/2, except for the edge cases 
    ref1 = []; ch = []; ref2 = [];
    if verbose, fprintf('Input leads:\n'); end
    for ll = 1:length(leads)
        
        % position of contacts (numbers) along the current lead.
        [leadPos, ix] = sort(channelPos(strcmp(channelChar, leads{ll})));
        % Indices into input channels matching leadPos
        leadInds = find(strcmp(channelChar, leads{ll})); leadInds = leadInds(ix);
        
        % Ignore single-contact leads (e.g., EKG)
        if length(leadPos) == 1, continue; end
        
        assert(length(unique(leadPos)) == length(leadPos), 'Repeated lead position in %s', leads{ll});
        
        if inclEdges % edges are bipolar referenced
            chLead = leadInds(1:end);
            ref1Lead = leadInds([2, 1:end-1]); % repeat pos2 as bipolar reference for the first edge
            ref2Lead = leadInds([2:end, end-1]); % repeat pos2 from the end as bipolar reference for last edge
            posReferenced = [leadPos([2, 1:end-1]), leadPos(1:end)]; % which positions along the lead are used as ref1 and ch
        else % only return non-edge channels
            chLead = leadInds(2:end-1);
            ref1Lead = leadInds(1:end-2);
            ref2Lead = leadInds(3:end);
            posReferenced = [leadPos(1:end-2), leadPos(2:end-1)];
        end
        
        % remove cross-segment bipolar pairs
        if ~isempty(seg5) && any(strcmpi(leads{ll}, seg5))
            if mod(leadPos(end), 5), warning('%s lead is marked as segmented in 5s, but it is length-%d (non-multiple)\n', leads{ll}, leadPos(end)); end
            
            % cross-gap rereferenced channels
            crossgap1 = ~mod(posReferenced(:, 1), 5); % between ref1 and ch (e.g., 5, 6, 7)
            crossgap2 = ~mod(posReferenced(:, 2), 5); % between ch and ref2 (e.g., 4, 5, 6)
            
            if inclEdges % set to bipolar
                ref1Lead(crossgap1) = ref2Lead(crossgap1); % replace ref1 with copy of ref2 for ref1-ch gaps, e.g. 7, 6, 7
                ref2Lead(crossgap2) = ref1Lead(crossgap2); % replace ref2 with copy of ref1 for ch-ref2 gaps, e.g. 4, 5, 4
            else
                chLead(crossgap1 | crossgap2) = [];
                ref1Lead(crossgap1 | crossgap2) = [];
                ref2Lead(crossgap1 | crossgap2) = [];
            end
            
            if verbose, fprintf('Lead %s, contacts 1 ... %d, segmented in 5s\n', leads{ll}, max(leadPos)); end
        elseif ~isempty(seg6) && any(strcmpi(leads{ll}, seg6))
            if mod(leadPos(end), 6), warning('%s lead is marked as segmented in 6s, but it is length-%d (non-multiple)\n', leads{ll}, leadPos(end)); end
            
            % cross-gap rereferenced channels
            crossgap1 = ~mod(posReferenced(:, 1), 6); % between ref1 and ch (e.g., 5, 6, 7)
            crossgap2 = ~mod(posReferenced(:, 2), 6); % between ch and ref2 (e.g., 4, 5, 6)
            
            if inclEdges % set to bipolar
                ref1Lead(crossgap1) = ref2Lead(crossgap1); % replace ref1 with copy of ref2 for ref1-ch gaps, e.g. 7, 6, 7
                ref2Lead(crossgap2) = ref1Lead(crossgap2); % replace ref2 with copy of ref1 for ch-ref2 gaps, e.g. 4, 5, 4
            else
                chLead(crossgap1 | crossgap2) = [];
                ref1Lead(crossgap1 | crossgap2) = [];
                ref2Lead(crossgap1 | crossgap2) = [];
            end
            
            if verbose, fprintf('Lead %s, contacts 1 ... %d, segmented in 6s\n', leads{ll}, max(leadPos)); end
        else
            if verbose, fprintf('Lead %s, contacts 1 ... %d\n', leads{ll}, max(leadPos)); end
        end
        
        % add which indices to subtract from which
        ref1 = [ref1; ref1Lead];
        ch = [ch; chLead];
        ref2 = [ref2; ref2Lead];
        
    end
    
    dataOut = nan(size(dataIn));
    dataOut(:, ch) = dataIn(:, ch) - (dataIn(:, ref1) + dataIn(:, ref2))/2;
    
    badChansLap = unique([badChans; ch(ismember(ref1, badChans)); ch(ismember(ref2, badChans))]); % bad channel if ref1 and ref2 are also inclusive of bad channels
            
    if verbose, fprintf('Constructed %d CSD channels out of %d\n', length(ch), length(channelNames)); end
    
end