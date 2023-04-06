%
%   Applies common average referencing to CCEP data from blocks of 64 channels. Assumes that noise is shared between
%       contiguous 64-channel blocks, starting at channel 1. Don't remove bad channels before inputting to this code.
%       Non-SEEG channels are returned unchanged.
%   Excludes:
%       - Channels being stimulated in each trial
%       - Channels with high variance from 500-1000 ms after stimulation onset
%         = channels with large stimulation artifacts)
%       - Channels with high ratio of (variance from 15-100ms) / (variance from 500-1000ms)
%         = channels with significant evoked responses
%
%   signaldata = ccep_CAR64blocks(signaldata, tt, chTbl, stimNames)
%   signaldata = ccep_CAR64blocks(signaldata, tt, chTbl, stimNames, pctThresh, ratioThresh, minChs)
%
%       signaldata =            n x T x m array, with n channels, T time periods, m trials/epochs
%       tt =                    1xT array of T timepoints
%       chTbl =                 n x __ table of channel info; columns must include "name", "type", and "status"
%       stimNames =             m x 1 cell array, names of stim sites as channels separated by hyphen. E.g. "LG1-LG2"
%       pctThresh =             (optional, default = 95) double [0 - 100], percentile cutoff for variance
%       ratioThresh =           (optional, default = 1.5) double > 1, threshold cutoff for significant evoked potential
%       minChs =                (optional, default = 8) double, mininum number of channels to create CAR with.
%                                   If a block has fewer than minChs number of channels to include, the CAR will be set
%                                   to an nan signal
%       carFirst =              (optional, default = false) lobical, if true, a CAR across all channels will first be applied before estimating variances. This is so 
%                                   that variance estimates are not primarily drowned out by noise.
%
%   Returns:
%       signaldata =            n x T x m array, signals in input signaldata after subtracting CARs
%
% DH and HH Multimodal Neuroimaging Lab, Mayo Clinic, 2020
%
function signaldata = ccep_CAR64blocks(signaldata, tt, chTbl, stimNames, pctThresh, ratioThresh, minChs, carFirst)

    if nargin < 8, carFirst = false; end
    if nargin < 7 || isempty(minChs), minChs = 8; end
    if nargin < 6 || isempty(ratioThresh), ratioThresh = 1.5; end
    if nargin < 5 || isempty(pctThresh), pctThresh = 95; end
    assert(ismember(1, find(strcmp(chTbl.type, 'SEEG'))), 'SEEG channels must begin at index 1'); 
    
    nonSEEG_channels = find(~strcmp(chTbl.type, 'SEEG'));
    good_channels = find(strcmp(chTbl.status, 'good') & strcmp(chTbl.type, 'SEEG'));
    
    numVarExclude = zeros([size(signaldata, 3), 1]); % track how many channels are variance-excluded per trial
    numRespExclude = zeros([size(signaldata, 3), 1]); % track how many channels are response-excluded per trial
    for kk = 1:size(signaldata, 3)

        stimChs = find(ismember(chTbl.name, split(stimNames(kk), '-'))); % channels being stimulated
        
        if carFirst
            sigdataTemp = signaldata(:, :, kk) - mean(signaldata(setdiff(good_channels, stimChs), :, kk), 1); % do a first-round prelim CAR on all good channels to more precisely evaluate variance
        else
            sigdataTemp = signaldata(:, :, kk);
        end
        
        % calculate variances on prelim-CAR'ed data
        var_late = var(sigdataTemp(:, tt > 0.5 & tt < 1), [], 2); % measure of noise in signal
        var_early = var(sigdataTemp(:, tt > 0.015 & tt < 0.1), [], 2); % measure of response in signal

        maxVar = prctile(var_late(setdiff(good_channels, stimChs)), pctThresh); % max var to accept, calculated on good, non-stim chs
        maxEarlyVar = prctile(var_early(setdiff(good_channels, stimChs)), 75); % require at least 75% of good chs to pass response thresh (assumes true resp in <25% chs)
        
        chs_incl = setdiff(good_channels, stimChs); % exclude stimulated channels
        numVarExclude(kk) = length(intersect(find(var_late > maxVar), chs_incl)); % record how many channels are being excluded by each filter
        chs_incl = setdiff(chs_incl, find(var_late > maxVar)); % exclude noisy channels
        numRespExclude(kk) = length(intersect(find(var_early./var_late > ratioThresh & var_early > maxEarlyVar), chs_incl));
        chs_incl = setdiff(chs_incl, find(var_early./var_late > ratioThresh & var_early > maxEarlyVar)); % exclude channels with responses

        car_sets = zeros(ceil(size(signaldata, 1)/64), size(signaldata, 2)); % CAR for each 64-channel block (using only the included channels)
        for set_nr = 1:size(car_sets, 1)
            this_car_set = intersect((set_nr*64 - 63):set_nr*64, chs_incl);
            if length(this_car_set) >= minChs
                car_sets(set_nr, :) = mean(signaldata(this_car_set, :, kk), 1);
            else
                if length(intersect((set_nr*64 - 63):set_nr*64, good_channels)) >= minChs % print only if the block would've passed before filtering
                    fprintf('Trial %d: CAR block #%d set to nan for having %d channels\n', kk, set_nr, length(this_car_set));
                end
                car_sets(set_nr, :) = nan([1, size(signaldata, 2)]); % not enough channels in CAR set
            end
        end

        for ii = 1:size(signaldata, 1)
            if ismember(ii, nonSEEG_channels)
                continue % don't modify non-SEEG channel
            else
                set_nr = ceil(ii/64);
                signaldata(ii, :, kk) = signaldata(ii, :, kk) - car_sets(set_nr, :);
            end
        end

    end
    
    fprintf('Average #chs excluded per trial (total %d good channels):\n\tby variance = %.0f (%.0f%%)\n\tadditionally by response = %.0f (%.0f%%)\n', ...
                length(good_channels), mean(numVarExclude), 100*mean(numVarExclude)/length(good_channels), mean(numRespExclude), 100*mean(numRespExclude)/length(good_channels));
end