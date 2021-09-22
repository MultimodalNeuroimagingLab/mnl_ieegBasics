%% This function performs bipolar referencing for data recorded from SEEG channels, based on input channel names.
%
%   Bipolar signals are calculated as the signal difference between two adjacent electrode contacts on the same SEEG lead.
%   To determine which channels should be subtracted from each other, this function requires the channel names to be in
%   the format 'XXXnn', where 'XXX' = chars of any length that represent the unique name of an SEEG lead, and 'nn' is
%   the electrode position on that lead, starting from 1. E.g. 'LA12' is the 12th electrode along the 'LA' lead.
%
%   The bipolar signals are then calculated starting from the 1st electrode in a lead upward, in the order of
%   subtracting the smaller odd number lead from the larger even number lead. E.g. LA1 - LA2, LA3 - LA4 ...
%   There is no shared component between any pair of returned bipolar signals. Any remaining singular (odd) electrode at
%   the end of a lead (e.g. LA15) is not included in the bipolar derivations. Similarly, the user can specify a
%   numerical list of bad channels, and any pairs containing a bad channel will not be included in the bipolar
%   derivations. "Good" channels that aren't included in the bipolar derivation (i.e. those at the end of leads, those 
%   paired to bad channels) will have their indices returned in <excludedChans>.
%
%   Usage:
%   [dataOut, bipolarChans] = ieeg_bipolarSEEG(dataIn, channelNames, []);
%   [dataOut, bipolarChans, excludedChans] = ieeg_bipolarSEEG(dataIn, channelNames, badChans);
%       dataIn =            txn double. Signal data: rows are samples, columns are channels. Signal should be already
%                               highpassed so that the amplitude of adjacent channels are comparable.
%       channelNames =      nx1 cell{char array}. Channel names must be in the same order as columns of <dataIn>. See
%                               the first paragraph of the documentation for how channel names should be formatted.
%       badChans =          mx1 num, m <= n. Indices of bad channels to not use in bipolar derivation. Channels paired
%                               to bad channels will be listed in <excludedChans> if they are themselves not bad.
%
%   Returns:
%       dataOut =           txb double. Bipolar-referenced signal data. Rows are samples, columns are re-referenced
%                               channels. b will be <= n/2, because original channels have been paired up to form the
%                               rereferenced data.
%       bipolarChans =      bx2 cell. 1st column contains channel pair names corresponding to columns in <dataOut>.
%                               The names take the form 'ch1-ch2', e.g. 'LA1-LA2' to represent which pair of channels
%                               are used for the derivation. 2nd column contains 1x2 array of input channel indices
%                               matching the channel pair names.
%       excludedChans =     px1 num. Indices of good channels in the original <dataIn> that were excluded during bipolar
%                               derivation. See second paragraph of docstring. [badChans; excludedChans] will give
%                               the complete list of all indices to SEEG channels not used in bipolar derivation. Note,
%                               this list will not include channels corresponding to incorrectly-formatted names, such
%                               as 'ECG' or '12'.
%
%   HH 2021
%
function [dataOut, bipolarChans, excludedChans] = ieeg_bipolarSEEG(dataIn, channelNames, badChans)
    
    isDigit = @(x) x > 47 & x < 58; % returns true for char array elements that are digits (0 - 9)
    
    assert(length(channelNames) == size(dataIn, 2), 'Number of columns in dataIn must match number of channel names');
    
    emptyChs = cellfun(@isempty, channelNames);
    channelNames(emptyChs) = []; dataIn(:, emptyChs) = [];
    channelNames = strip(upper(channelNames)); % some cleaning
    
    leads = unique(cellfun(@(ch) ch(~isDigit(ch)), channelNames, 'UniformOutput', false)); % unique lead names
    leads(cellfun(@isempty, leads)) = []; % discard any channel names that don't contain letters
    
    channelPos = cellfun(@(ch) str2double(ch(isDigit(ch))), channelNames); % electrode position along lead for each channel
    channelChar = cellfun(@(ch) ch(~isDigit(ch)), channelNames, 'UniformOutput', false); % just the character name of each lead
    
    minuend = []; subtrahend = []; % which indices to subtract from which
    bipolarChans = cell(0, 2); % 1st col = name, 2nd col = array of index pairs
    excludedChans = []; % indices of good channels skipped over
    for ll = 1:length(leads)
        
        [leadPos, ix] = sort(channelPos(strcmp(channelChar, leads{ll})));
        % Indices into input channels matching leadPos
        leadInds = find(strcmp(channelChar, leads{ll})); leadInds = leadInds(ix);
        
        fprintf('Lead %s, pos 1 ... %d\n', leads{ll}, max(leadPos));
        
        assert(length(unique(leadPos)) == length(leadPos), 'Repeated lead position in %s', leads{ll});
        
        pp = 1; % position of minuend along current lead
        while pp+1 <= max(leadPos)
            
            chPair = [leadInds(leadPos == pp); leadInds(leadPos == pp+1)]; % [minuend, subtrahend]
            
            if length(chPair) < 2 % both positions were not found
                excludedChans = [excludedChans; chPair];
                pp = pp + 2;
                continue
            end
            
            % check if either of bipolar pair is bad channel
            badBool = ismember(chPair, badChans);
            if any(badBool)
                excludedChans = [excludedChans; chPair(~badBool)]; % append non-bad channels for record keeping
                pp = pp + 2;
                continue
            end
            
            % e.g. LA1-LA2
            minuend = [minuend; chPair(1)];
            subtrahend = [subtrahend; chPair(2)];
            bipolarChans = [bipolarChans; sprintf('%s%d-%s%d', leads{ll}, pp, leads{ll}, pp+1), {chPair'}];
            
            pp = pp + 2;
        end
    end
    
    fprintf('Constructed %d bipolar channels. Excluded %d good channels.\n', length(minuend), length(excludedChans));
    dataOut = dataIn(:, minuend) - dataIn(:, subtrahend);

end