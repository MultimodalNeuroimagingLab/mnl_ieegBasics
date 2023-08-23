function [Vout, chsUsed] = ccep_CARVariance(tt, V, srate, badChs, grp, optsIn)
%
%   Applies adjusted common average referencing to iEEG using a subset of channels with lowest mean trial variance or cross-trial covariance on a response interval.
%   This function was created as a bridge between the old ccep_CAR functions in mnl_ieegBasics and the newer CARLA. An important change over previous methodology
%   is that this function uses the same channels as common average across all trials of the same group (e.g., stimulation site). The common average is still calculated
%   on a per-trial basis. This function was created with CCEP data in mind, but can be extrapolated to cognitive task data.
%
%   Vout = CARVariance(tt, V); % minimum inputs, NOT RECOMMENDED
%   Vout = CARVariance(tt, V, srate); % recommended minimum inputs
%   [Vout, chsUsed] = CARVariance(tt, V, srate, badChs, grp); % most common usage
%   [Vout, chsUsed] = CARVariance(tt, V, srate, badChs, grp, opts); % full usage
%       tt =        1 x t num. Time points matching V, in seconds
%       V =         n x t x k num, Signal data to be re-referenced, with n channels by t time points by k trials.
%                       If there are any singular bad trials that shouldn't be considered in CAR, set them to nan in V. Similarly, channels that shouldn't be
%                       considered only for a group (e.g. stimulated electrodes) should be set to nan for all trials corresponding to that group.
%                       Percentile calculation ignores channels that are set to nan this way. Any nans at single time points will be propagated to the entire trial(s) containing the nan time point.       
%       srate =     num. Sampling frequency of data, used for notch filter if opts.notchfirst==true.
%                       If not given and opts.notchfirst==true, this will be estimated from the tt input (not recommended because of possible imprecision).
%       badChs =    (optional) 1 x m num, m < n. List of all bad channels to not consider for CAR (e.g., [36, 78, 100, 212:230]).
%       grp =       (optional) 1 x k num, string array, or cell array(char). Grouping variable for trials. All trials in the same group must have identical values (numerical or character).
%                       If not given, assumes the entire input is one group
%       optsIn =    (optional) struct to adjust parameters, fields described below:
%           pctThresh =     num, between 0 and 100. Percentile threshold on how many channels to use for common average. Default = 25
%           winResp =       1 x 2 num. [start, stop] of responsive time period to calculate variance and correlations on, in seconds matching tt. Default = [0.01, 0.3].
%           varType =       char, 'var' or 'cov'. What metric to use when ranking channels. 'var' ranks channels based on (geometric) average variance across trials.
%                               'cov' ranks channels based on average cross-trial covariance. Default = 'cov'.
%           notchFirst =    bool. Whether notch filters are applied at 60, 120, 180 Hz to reduce line noise on the data before ranking/optimizing by CARLA.
%                               This only affects determining which channels to include, as the rereferenced output is calculated from the unfiltered input. Default = true.
%
%   Returns:
%       Vout =            n x t x k num, or n x t num. Re-referenced signal data.
%       chsUsed =         q x 1 cell array num. For each of the q unique groups in <grp>, the corresponding cell contains the list of channel numbers used to construct the
%                               common average for that group, in order of increasing <varType>.
%
% HH, DD Multimodal Neuroimaging Lab, Mayo Clinic, 2023
%
    %% Input checks, data formatting

    nTrs = size(V, 3); % number of trials
    assert(size(V, 2) == length(tt), 'Error: second dimension does not match tt. Data should be channels x timepoints x trials');
    
    % propagate any single nan timepoints in input to entire trial(s) affected
    nanmask = repmat(any(isnan(V), 2), 1, length(tt), 1);
    V(nanmask) = nan;
    
    % assume all are the same stim condition if none given
    if nargin < 5 || isempty(grp), grp = ones(nTrs, 1); end
    if nargin < 4, badChs = []; end % no bad channels
    
    assert(isnumeric(badChs), 'Error: badChs must be given as a numeric list');
    
    if isstring(grp), grp = cellstr(grp); end % convert to cell array of chars if given as string array
    assert(length(grp) == nTrs, 'Error: grp input does not match trials in length');

    % Default opts
    opts.pctthresh = 25;
    opts.winresp = [0.01, 0.3];
    opts.vartype = 'cov';
    opts.notchfirst = true;
    
    % Configure opts if input is given
    if exist('optsIn', 'var') && isa(optsIn, 'struct')
        inFields = fieldnames(optsIn);
        for ii = 1:length(inFields) % make case insensitive

            optsIn.(lower(inFields{ii})) = optsIn.(inFields{ii});
            fieldCurr = lower(inFields{ii});
            
            assert(any(strcmp(fieldCurr, {'pctthresh', 'winresp', 'vartype', 'notchfirst'})), 'Error: "%s" is not a valid opts field', fieldCurr);
            
            if strcmp(fieldCurr, 'alpha'), assert(optsIn.(fieldCurr) > 0 & optsIn.(fieldCurr) < 100, 'Error: pctthresh must be between 0 and 100'); end
            if strcmp(fieldCurr, 'winresp'), assert(length(optsIn.(fieldCurr)) == 2 & all(optsIn.(fieldCurr) > 0) & all(optsIn.(fieldCurr) < tt(end)), ...
                    'Error: winresp must be given as [t1, t2] where both are greater than 0 and before tt(end)'); end
            if strcmp(fieldCurr, 'vartype'), assert(any(strcmp(optsIn.(fieldCurr), {'var', 'cov'})), 'Error: vartype has to be "var" or "cov"'); end
            if strcmp(fieldCurr, 'notchfirst'), assert(islogical(optsIn.(fieldCurr)), 'Error: notchfirst must be logical (true or false)'); end           
            
            opts.(fieldCurr) = optsIn.(fieldCurr);
        end
    end
    
    % Estimate sampling frequency from time points if necessary (not recommended)
    if nargin < 3 || isempty(srate) && opts.notchfirst
        srate = (length(tt) - 1) / (max(tt) - min(tt));
        warning('sampling frequency for notch filter was ESTIMATED to be %0.02f', srate);
    end
    
    % Remove bad channels
    VgoodChs = V;
    VgoodChs(badChs, :, :) = [];
    goodChs = setdiff(1:size(V, 1), badChs); % keep track of indices of good chs to put back later
    
    nChs = size(VgoodChs, 1);
    fprintf('%d good channels out of %d total channels considered for CAR\n', nChs, size(V, 1));
    
    fprintf('Opts: %d percent of channels, [%0.03f - %0.03f] ms window, using %s, notchFirst=%s\n', ...
        opts.pctthresh, opts.winresp(1), opts.winresp(2), opts.vartype, string(opts.notchfirst));
        
    %% Apply notch filter to reduce line noise (if requested)
    
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
    
    %% Loop through trial conditions to determine which channels to use for each condition
    
    % extract data on response interval to perform analysis on
    Vseg = Vclean(:, tt >= opts.winresp(1) & tt <= opts.winresp(end), :);
        
    grpTypes = unique(grp);
    chsUsed = cell(length(grpTypes), 1);
    
    fprintf('Conditions: ');
    for ii = 1:length(grpTypes)
        
        if isnumeric(grp) % handle numeric and string (not char array) inputs
            VsegCurr = Vseg(:, :, grp == grpTypes(ii));
            condStr = num2str(grpTypes(ii));
        elseif iscell(grp)
            VsegCurr = Vseg(:, :, strcmp(grp, grpTypes{ii}));
            condStr = grpTypes{ii};
        end
        fprintf('%s (%d), ', condStr, size(VsegCurr, 3));
        if ~mod(ii, 5), fprintf('\n'); end
        
        if strcmp(opts.vartype, 'var') % Geomean variance across trials
            vars = geomean(var(VsegCurr, 0, 2), 3, 'omitnan'); % omitnan ignores bad trials set to nan
        elseif strcmp(opts.vartype, 'cov') % Mean covariance
            if size(VsegCurr, 3) < 6, warning('Mean covariance may be unreliable for few trials due to spurious anticorrelations between trials. Use ''var'' instead'); end
            vars = nan(nChs, 1);
            for jj = 1:nChs
                covCurr = cov(squeeze(VsegCurr(jj, :, :)));
                vars(jj) = mean(covCurr(logical(triu(ones(size(covCurr)), 1))), 'all', 'omitnan'); % average the covariances above the diagonal, ignoring bad trials (nan-coded)
            end
        end
        nanChs = find(isnan(vars));
        if ~isempty(nanChs)
            warning('All trials are nan at channels %s, not considered for CAR at current condition', join(num2str(goodChs(nanChs))));
        end
        
        % variance threshold for current condition
        thresh = prctile(vars, opts.pctthresh);
        
        % sort good channels in increasing order of variance so easier to troubleshoot
        [vars, ord] = sort(vars);
        goodChsSorted = goodChs(ord);
        chsUsed{ii} = goodChsSorted(vars < thresh); % which channels fall under threshold, by original ch number, sorted in order of increasing variance
                
    end
    fprintf('\n');
    
    %% Perform CAR on each condition
    
    Vout = V;
    
    for ii = 1:length(grpTypes)
        
        if isnumeric(grp) % handle numeric and string (not char array) inputs
            inds = grp == grpTypes(ii);
        elseif iscell(grp)
            inds = strcmp(grp, grpTypes{ii});
        end
        
        Vout(:, :, inds) = Vout(:, :, inds) - mean(Vout(chsUsed{ii}, :, inds), 'omitnan'); % omitnan ignores bad trials when calculating common average
    end
    
end