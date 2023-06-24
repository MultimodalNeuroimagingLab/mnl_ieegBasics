function [Vout, CA, stats] = ccep_CARLA(tt, V, srate, badChs, optsIn)
%
%   This function performs Common Average Re-referencing by Least Anticorrelation (CARLA) for CCEP data.
%   1. Channels are ranked in order of increasing variance or covariance across trials. 
%   2. Channels are iteratively added into a common average, and the mean correlation between each unre-referenced channel and all re-referenced channels is 
%       calculated, for each bootstrapped mean signal. zminmean denotes the mean correlation belonging, at each bootstrap, to the most anticorrelated channel.
%   3. The optimum common average size corresponds to the least anticorrelation (when zminmean takes on its least negative value)
%
%   Vout = ccep_CARLA(tt, V); % minimum inputs, NOT RECOMMENDED
%   Vout = ccep_CARLA(tt, V, srate); % recommended minimum inputs
%   [Vout, CAR, stats] = ccep_CARLA(tt, V, srate, badChs, opts); % full usage
%       tt =        1 x t num. Time points matching V, in seconds
%       V =         n x t x k num, or n x t num. Signal data to be re-referenced, with n channels by t time points by k trials (2 dimensional if k = 1).
%                       Any single time points that are nan are propagated to entire trial(s) containing those nan time points
%       srate =     num. Sampling frequency of data, used for notch filter if opts.notchfirst==true.
%                       If not given and opts.notchfirst==true, this will be estimated from the tt input (not recommended because of possible imprecision).
%       badChs =    (optional) 1 x m num, or m x 2 cell array, m < n. If array, this is a simple list of all bad channels (e.g., [36, 78, 100, 212:230]).
%                       If given as m x 2 cell array, the first column lists all channels that are either bad (for all trials) or bad only at certain trials. The
%                       second column contains arrays listing which trials are bad for that channel. These channels WILL be considered in finding the best CAR.
%                       Empty array indicates the channel is bad for all trials.
%                       E.g., {[10], [];          % This input means that channels 10 and 11 are bad for all trials, but channel 155 is only bad at trials 4, 7, and 8.
%                              [11], [];          % So the good trials of channel 155 will still be considered in determining the common average
%                              [155], [4, 7, 8]}
%       optsIn =      (optional) struct to adjust parameters, fields described below:
%           varType =       char, 'var' or 'cov'. What metric to use when ranking channels. 'var' ranks channels based on (geometric) average variance across trials.
%                               'cov' ranks channels based on average cross-trial covariance. 'var' is used when k = 1 regardless of input. Default = 'cov'.
%           notchFirst =    bool. Whether notch filters are applied at 60, 120, 180 Hz to reduce line noise on the data before ranking/optimizing by CARLA.
%                               This only affects determining which channels to include, as the rereferenced output is calculated from the unfiltered input. Default = true.
%           nBoot =         integer. Number of bootstrapped mean signal samples to generate when calculating zminmean. Ignored if k = 1. Default = 100.
%           winResp =       1 x 2 num. [start, stop] of responsive time period to calculate variance and correlations on, in seconds matching tt. Default = [0.01, 0.3].
%           alpha =         num. Alpha threshold for significance when determining whether zminmean significantly decreases from local max to next trough. Default = 0.05.
%
%   RETURNS:
%       Vout =      n x t x k num, or n x t num. Re-referenced signal data.
%       CA =        t x k num, or 1 x t num. The common average calculated for each trial
%       stats =     struct, containing fields:
%                       chsUsed =       p x 1 num, p <= n. The sorted n channel numbers used in the optimum CAR (indexing the original input)
%                       vars =          n x 1 num. Cross-trial covariance or geomean variance of each good channel (depending on optsIn.vartype)
%                       order =         n x 1 num. The good channels sorted in order of increasing (co)variance
%                       zMin =          n x n x nboot num, or n x n num if k = 1. The z (Fisher z-transformed Pearson's r) for the channel with the most negative 
%                                           average z, at each bootstrap. Dimensions are (z to each re-referenced channel) x (number of channels in CAR subset) x (bootstrap)
%                       zMinMean =      n x nboot num, or 1 x n num if k = 1. Mean zMin across the first dimension (averaged across all "target" re-referenced channels)
%                                           for each bootstrap if V was 3 dimensional
%
%   First iteration of CARLA created 2023/04/27 by Harvey Huang   
%   This function was adapted from CARLA on 2023/06/14 by Harvey Huang
%   
    %% Input checks, data formatting
    
    % Ensure ch x time points structure if k = 1
    if ismatrix(V)
        if size(V, 1) == length(tt), V = V'; end
    end
    
    % propagate single nan timepoints to entire trial(s) affected
    nanmask = repmat(any(isnan(V), 2), 1, length(tt), 1);
    V(nanmask) = nan;
    
    nTrs = size(V, 3); % number of trials
    assert(size(V, 2) == length(tt), 'Error: second dimension does not match tt. Data should be channels x timepoints x trials');
    
    if nargin < 4, badChs = []; end
    assert(isnumeric(badChs), 'Error: badChs must be given as a numeric list');
    
    % default opts
    opts.vartype = 'cov';
    opts.notchfirst = true;
    opts.nboot = 100;
    opts.winresp = [0.01, 0.3];
    opts.alpha = 0.05;
    
    % configure opts if input is given
    if exist('optsIn', 'var') && isa(optsIn, 'struct')
        inFields = fieldnames(optsIn);
        for ii = 1:length(inFields) % make case insensitive
            fieldCurr = lower(inFields{ii});
            assert(any(strcmp(fieldCurr, {'winresp', 'vartype', 'notchfirst', 'nboot', 'alpha'})), 'Error: "%s" is not a valid opts field', fieldCurr);
            
            if strcmp(fieldCurr, 'vartype'), assert(any(strcmp(optsIn.(fieldCurr), {'var', 'cov'})), 'Error: vartype has to be "var" or "cov"'); end
            if strcmp(fieldCurr, 'notchfirst'), assert(islogical(optsIn.(fieldCurr)), 'Error: notchfirst must be logical (true or false)'); end
            if strcmp(fieldCurr, 'nboot'), assert(optsIn.(fieldCurr) > 1, 'Error: nboot must be >1'); end
            if strcmp(fieldCurr, 'winresp'), assert(length(optsIn.(fieldCurr)) == 2 & all(optsIn.(fieldCurr) > 0) & all(optsIn.(fieldCurr) < tt(end)), ...
                    'Error: winresp must be given as [t1, t2] where both are greater than 0 and before tt(end)'); end
            if strcmp(fieldCurr, 'alpha'), assert(optsIn.(fieldCurr) > 0 & optsIn.(fieldCurr) < 1, 'Error: alpha must be between 0 and 1'); end
            
            opts.(fieldCurr) = optsIn.(fieldCurr);
        end
    end
    
    % Estimate sampling frequency from time points if necessary (not recommended)
    if nargin < 3 || isempty(srate) && opts.notchfirst
        srate = (length(tt) - 1) / (max(tt) - min(tt));
        warning('sampling frequency for notch filter was ESTIMATED to be %0.02f', srate);
    end
    
    % Remove channels bad for all trials, keep track of channels with selectively-bad trials
    VgoodChs = V;
    if isnumeric(badChs) % vector of badChs given (bad at all trials)
        VgoodChs(badChs, :, :) = [];
        badTrs = cell(size(VgoodChs, 1), 1); % empty cell array to match behavior of cell input
        goodChs = setdiff(1:size(V, 1), badChs); % keep track of indices of good chs
    elseif iscell(badChs)
        assert(nTrs > 1, 'For single-trial inputs, badChs should be a numeric vector');
        assert(size(badChs, 2) == 2, 'Error: expecting badChs to have 2 columns if given as cell array');
        for ii = 1:size(badChs, 1) % make each bad trial list a column vector for consistency when concatenating
            badChs{ii, 2} = badChs{ii, 2}(:);
        end
        trsListed = unique(vertcat(badChs{:, 2})); % all trials listed
        assert(all(trsListed >= 1 & trsListed <= nTrs), 'Error: Bad trials for channels are outside of trial indices');
        
        % rows in badChs with empty cells; or where all but one or all trials are mentioned. Can't calculate covariance with only 1 good trial
        badAllTrsBool = cellfun(@(x) isempty(x) | length(unique(x)) >= nTrs-1, badChs(:, 2));
        badChsAll = cell2mat(badChs(badAllTrsBool, 1));
        goodChs = setdiff(1:size(V, 1), badChsAll); % keep track of inverse of badChsAll
        
        % stores which trials are bad for the channels that are only bad at a few trials, corresponding to rows of VgoodChs
        badTrs = cell(size(V, 1), 1); 
        badTrs(cell2mat(badChs(~badAllTrsBool, 1))) = badChs(~badAllTrsBool, 2);
        
        VgoodChs(badChsAll, :, :) = [];
        badTrs(badChsAll) = [];
    else
        error('Error: badChs must be given as an integer vector or as an nx2 cell array');
    end
    
    nChs = size(VgoodChs, 1);
    
    fprintf('%d good channels out of %d total channels considered for CAR\n', nChs, size(V, 1));
    
    %% Rank channels by increasing variance or covariance
    
    stats = struct();
    
    % Notch filter signals to remove line noise, if desired. This is not used for final output.
    Vclean = VgoodChs;
    if opts.notchfirst
        fprintf('Applying notch filters: ');
        for ff = 60:60:180
            fprintf('%d Hz, ', ff);
            dNotch = designfilt('bandstopiir', 'FilterOrder', 4, ... % matches filter in ieeg_notch
                                    'DesignMethod', 'butter', ...
                                    'HalfPowerFrequency1', ff-2, ... 
                                    'HalfPowerFrequency2', ff+2, ...
                                    'SampleRate', srate);
            for ii = 1:nTrs % filter per trial
                VcleanTr = Vclean(:, :, ii); % all channels at current trial
                VcleanTr(~any(isnan(VcleanTr), 2), :) = filtfilt(dNotch, VcleanTr(~any(isnan(VcleanTr), 2), :)')'; % filter channels with non-nan time points
                Vclean(:, :, ii) = VcleanTr;
            end
        end
        fprintf('\n');
    end
    
    % extract data on response interval to perform analysis on
    Vseg = Vclean(:, tt >= opts.winresp(1) & tt <= opts.winresp(end), :);
    
    % set bad trials for each channel to nan
    for ii = 1:nChs
        Vseg(ii, :, badTrs{ii}) = nan;
        VgoodChs(ii, :, badTrs{ii}) = nan; % apply to VgoodChs too because this will be averaged to make the CAR
    end
    
    if nTrs == 1 % single trial
        stats.vars = var(Vseg, 0, 2);
    elseif strcmp(opts.vartype, 'var') % multiple trials but choose to use variance. Geomean across trials to be more prone to outliers
        stats.vars = geomean(var(Vseg, 0, 2), 3, 'omitnan'); % omitnan ignores bad trials for channels in calculation
    elseif strcmp(opts.vartype, 'cov') % multiple trials, choose to use mean covariance
        if nTrs < 6, warning('Mean covariance may be unreliable for few trials due to spurious anticorrelations between trials. Recommend using ''var'' instead'); end
        stats.vars = nan(nChs, 1);
        for ii = 1:nChs
            covCurr = cov(squeeze(Vseg(ii, :, :)));
            stats.vars(ii) = mean(covCurr(logical(triu(ones(size(covCurr)), 1))), 'all', 'omitnan'); % average the covariances above the diagonal, ignoring bad trials (nan-coded)
        end
    end
    [~, stats.order] = sort(stats.vars);
    
    %% Perform CARLA optimization
    
    % channels x CARsize x bootstrapped sample
    if nTrs == 1
        stats.zMin = nan(nChs, nChs); % no bootstrapping if only single trial
    else
        stats.zMin = nan(nChs, nChs, opts.nboot);
    end
        
    % Pull subset U (order(1:ii)) from V, growing U each time by 1, to construct CAR
    for ii = 2:nChs % min CAR size is 2
        
        % Vseg after subtracting the putative CAR on subset U
        VsegReref = Vseg - mean(Vseg(stats.order(1:ii), :, :), 1, 'omitnan'); % omitnan ignores individual channels in average for certain trials
        
        if nTrs == 1 % no bootstrapping, test correlation at the single trial level
            Useg = Vseg(stats.order(1:ii), :)'; % UNREREFERENCED channels in subset to correlate with the rereferenced signals
            UsegReref = VsegReref(stats.order(1:ii), :)'; % rereferenced signals in subset
            
            r = corr(Useg, UsegReref); % rows = which unrereferenced channel was used.
            r(1:(ii+1):end) = nan; % omit self-self correlations
            z = atanh(r); % fisher z-transform
            
            % choose "most responsive" channel with greatest anticorrelation, on average
            [~, kkMost] = min(mean(z, 2, 'omitnan'));
            stats.zMin(stats.order(1:ii), ii) = z(kkMost, :);
            
            continue;
            
        end
        
        % if nTrs > 1
        for jj = 1:opts.nboot
            
            % bootstrap to calculate mean signal
            inds = datasample(1:nTrs, nTrs);
            
            % omitnan ignores bad trials sampled in bootstrap for channels that have some bad trials. It is still possible to have ALL bad trials chosen in bootstrap, in which case ch will be nan
            Useg = mean(Vseg(stats.order(1:ii), :, inds), 3, 'omitnan')';
            UsegReref = mean(VsegReref(stats.order(1:ii), :, inds), 3, 'omitnan')';
            
            r = corr(Useg, UsegReref); % rows = which unrereferenced signal was used.
            r(1:(ii+1):end) = nan; % omit self-self correlations
            z = atanh(r); % fisher z-transform
            
            % choose "most responsive" channel with greatest anticorrelation, on average, across channels, for this bootstrapped sample
            [~, kkMost] = min(mean(z, 2, 'omitnan'));
            stats.zMin(stats.order(1:ii), ii, jj) = z(kkMost, :);
            
        end
        
    end
    
    % calculate mean across the re-referenced target channels
    stats.zMinMean = mean(stats.zMin, 1, 'omitnan');
    
    if nTrs > 1 % find first-peak optimum if it occurs before global optimum, when more than 1 trial
        
        nMin = ceil(0.1*nChs); % minimum number of channels to start at, for stability. Hardcoded at 10 %
        zMMxTrs = mean(stats.zMinMean(1, :, :), 3); % mean ZMinMean across bootstrapped samples for each n
        
        ii = nMin;
        while ii <= nChs
            
            if ii == nChs % at the end. no more comparisons needed
                nOptimum = ii;
                break
            end
            
            % loop until next local maximum found
            if zMMxTrs(ii + 1) > zMMxTrs(ii)
                ii = ii + 1;
                continue;
            end
            
            zMMxTrs(1:ii-1) = nan; % set all earlier values to nan, to avoid them being found as nextGreater or as trough
            
            % find index of the next greater zMinmean (ignores earlier nans), if it exists
            nextGreater = find(zMMxTrs > zMMxTrs(ii), 1, 'First');
            if isempty(nextGreater) % no more samples greater than current peak. Optimum is identical to global optimum. Done
                nOptimum = ii;
                break
            end
            zMMCurr = squeeze(stats.zMinMean(1, ii, :)); % all bootstrapped samples at current max
            [~, idx] = min(zMMxTrs(1:nextGreater-1)); % locate Zminmean trough before the next max (ignores earlier nans)
            zMMTrough = squeeze(stats.zMinMean(1, idx, :)); % all bootstrapped samples at trough before next greater value
            
            % all pairwise differences in the mean correlation estimate between the trough and current peak.
            % Note that zMinMean is the estimated correlation of the mean signal. We are testing the difference of zMinMean directly,
            % Not the mean of zMinMean, which would be redundant and inappropriate because the bootstrapped samples are not independent.
            [X, Y] = meshgrid(zMMTrough, zMMCurr);
            diffs = X - Y; diffs = diffs(:);
            conf = [-inf, prctile(diffs, 100*(1-opts.alpha))]; % left-tailed confidence interval (H_a = difference is less than 0)
            
            % found significant decrease. Current max is optimal.
            if conf(2) < 0
                nOptimum = ii;
                break
            end
            
            % Difference was not significant. Start at next greater sample
            ii = nextGreater;
            
        end
        
        if nOptimum == nMin, warning('Optimum CAR was detected at minimal 10% floor'); end
        
    else
        % Find optimum n at max (least negative) zminMean, averaged across trials.
        [~, nOptimum] = max(mean(stats.zMinMean, 3));
    end
    
    % delete first single dimension if 3D, after finding nOptimum
    stats.zMinMean = squeeze(stats.zMinMean);
    
    % which channels were used to make CAR (now indexing out of all channels input, including bad channels)
    stats.chsUsed = sort(goodChs(stats.order(1:nOptimum)));
    
    % calculate CAR from VgoodChs, which has nans in the bad trials. No filters applied to this.
    CA = mean(VgoodChs(stats.order(1:nOptimum), :, :), 'omitnan');
        
    Vout = V - CA;
    
    CA = squeeze(CA); % squeeze out the first dimension to make t x k, if 3 dimensional
    
end
