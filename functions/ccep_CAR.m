%
%   Applies common average referencing to CCEP data from blocks of 64 channels. Assumes that noise is shared between
%       contiguous 64-channel blocks, starting at channel 1. Non-ephys channels are returned unchanged
%   Excludes:
%       - Channels being stimulated in each trial
%       - Channels with high variance from 500-1000 ms after stimulation onset
%         = channels with large stimulation artifacts)
%       - Channels with high ratio of (variance from 15-100ms) / (variance from 500-1000ms)
%         = channels with significant evoked responses
%
%   signaldata = ccep_CAR64blocks(signaldata, tt, chTbl, stimNames)
%   signaldata = ccep_CAR64blocks(signaldata, tt, chTbl, stimNames, pctThresh, ratioThresh, carFirst)
%
%       signaldata =            n x T x m array, with n channels, T time periods, m trials/epochs
%       tt =                    1xT array of T timepoints
%       chTbl =                 n x __ table of channel info; columns must include "name", "type", and "status"
%       stimNames =             (optional), m x 1 cell array, names of stim sites as channels separated by hyphen. E.g. "LG1-LG2"
%       pctThresh =             (optional, default = 95) double [0 - 100], percentile cutoff for variance
%       ratioThresh =           (optional, default = 1.5) double > 1, threshold cutoff for significant evoked potential
%       carFirst =              (optional, default = false) lobical, if true, a CAR across all channels will first be applied before estimating variances. This is so 
%                                   that variance estimates are not primarily drowned out by noise.
%
%   Returns:
%       signaldata =            n x T x m array, signals in input signaldata after subtracting CARs
%
% Adapted from ccep_CAR64blocks on 2021/06/30
% DH and HH Multimodal Neuroimaging Lab, Mayo Clinic, 2020
%
function signaldata = ccep_CAR(signaldata, tt, chTbl, stimNames, pctThresh, ratioThresh, carFirst)

    if nargin < 7, carFirst = false; end
    if nargin < 6 || isempty(ratioThresh), ratioThresh = 1.5; end
    if nargin < 5 || isempty(pctThresh), pctThresh = 95; end
    
    ephys_channels = find(ismember(upper(chTbl.type), {'SEEG', 'ECOG'}));
    good_channels = find(strcmp(chTbl.status, 'good') & ismember(upper(chTbl.type), {'SEEG', 'ECOG'}));

    numVarExclude = zeros([size(signaldata, 3), 1]); % track how many channels are variance-excluded per trial
    numRespExclude = zeros([size(signaldata, 3), 1]); % track how many channels are response-excluded per trial
    for kk = 1:size(signaldata, 3)

        if nargin < 4 || isempty(stimNames)
            stimChs = [];
        else
            stimChs = find(ismember(chTbl.name, split(stimNames(kk), '-'))); % channels being stimulated
        end
        
        if carFirst
            sigdataTemp = signaldata(:, :, kk) - mean(signaldata(setdiff(good_channels, stimChs), :, kk), 1); % do a first-round prelim CAR on all good channels to more precisely evaluate variance
        else
            sigdataTemp = signaldata(:, :, kk);
        end
            
        % calculate variances on temp-CAR'ed signal data
        var_late = var(sigdataTemp(:, tt > 0.5 & tt < 1), [], 2); % measure of noise in signal
        var_early = var(sigdataTemp(:, tt > 0.015 & tt < 0.1), [], 2); % measure of response in signal

        maxVar = prctile(var_late(setdiff(good_channels, stimChs)), pctThresh); % max var to accept, calculated on good, non-stim chs
        maxEarlyVar = prctile(var_early(setdiff(good_channels, stimChs)), 75); % require at least 75% of good chs to pass response thresh (assumes true resp in <25% chs)
        
        chs_incl = setdiff(good_channels, stimChs); % exclude stimulated channels
        numVarExclude(kk) = length(intersect(find(var_late > maxVar), chs_incl)); % record how many channels are being excluded by each filter
        chs_incl = setdiff(chs_incl, find(var_late > maxVar)); % exclude noisy channels
        numRespExclude(kk) = length(intersect(find(var_early./var_late > ratioThresh & var_early > maxEarlyVar), chs_incl));
        chs_incl = setdiff(chs_incl, find(var_early./var_late > ratioThresh & var_early > maxEarlyVar)); % exclude channels with responses

        signaldata(ephys_channels, :, kk) = signaldata(ephys_channels, :, kk) - mean(signaldata(chs_incl, :, kk), 1); % subtract CAR from ephys channels

    end
    
    fprintf('Average #chs excluded per trial (total %d good channels):\n\tby variance = %.0f (%.0f%%)\n\tadditionally by response = %.0f (%.0f%%)\n', ...
                length(good_channels), mean(numVarExclude), 100*mean(numVarExclude)/length(good_channels), mean(numRespExclude), 100*mean(numRespExclude)/length(good_channels));
end