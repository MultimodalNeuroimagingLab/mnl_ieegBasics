%% This function compares time series data for trials wrt null(zero).
%   Cluster-based permutation tests, as described by Maris and Oostenveld (2007), test for significant spatiotemporal clusters in a timeseries data against a null hypothesis of zero effect,
%   using cluster mass (sum of t-values) as the statistic.
%   The null distribution was generated using sign-flip permutations for the one-sample test against zero
%   The approach supports one-sided (positive or negative effects only) or two-sided tests via tail specification.
%   
%   Please cite: Maris and Oostenveld (2007):
%       Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-and MEG-data. Journal of neuroscience methods, 164(1), 177-190.
%
%   The above paper also gives credit to Bullmore et al. (1999) for first conceiving the test statistic (tSum), where it was called the "Cluster mass test":
%       Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S., Taylor, E., & Brammer, M. J. (1999). Global, voxel, and cluster tests,
%           by theory and permutation, for a difference between two groups of structural MR images of the brain. IEEE transactions on medical imaging, 18(1), 32-42.
%
%   USAGE:
%       clusters = ieeg_clusterTest(A);    % minimal inputs/outputs
%       [clusters, pSig, tSumsSig] = ieeg_clusterTest(A, x2, alphaThresh, alphaFA, nperm, tail);   % full inputs/outputs
%           A =             t x n num. n Trials corresponding with t timepoints for the timeseries data to be compared against null
%           alphaThresh =   (optional) num, 0 < alphaThresh < 1. Alpha of the one-sided t-tests to detect samples with significant (uncorrected) difference in means between A and null.
%                               These samples are clustered on and then subjected to the non-parametric testing. alphaThresh does NOT affect the false alarm (FA)/FWER rate of the test, 
%                               but it does affect the sensitivity. Higher alphaThresh detects weaker/longer-lasting effects. Default = 0.05.
%           alphaFA =       (optional) num, 0 < alphaFA < 1. Alpha (false alarm (FA) rate) of the overall test. Sets the critical values of the one-sided permutation null distribution.
%                               This alpha controls the FA rate for ALL clusters tested (see section 3.1.1 last paragraph of Maris & Oostenveld, 2007). Default = 0.05.
%           nperm =         (optional) num > 0. Number of permutations in null distribution. Default = 1000.
%           tail  =         (optional) tail = 1(one-sided positive), -1(one-sided negative), or 0(two-sided). Default = 0.
%
%   RETURNS:
%       clusters =          b x 2 num. Each row is the [start, stop] indices in x1, x2 of a significant cluster. There are b significant clusters detected. No significant
%                               clusters if empty.
%       pSig =              b x 1 num. The permutation p-value for each significant cluster in clusters. Every element of pSig <= alphaFA.
%       tSumsSig =          b x 1 num. The cluster t-static sum for each significant cluster in clusters.
%           
% For comparision across 2 experimental conditions, see ieeg_clusterTest by Harvey Huang 08/2023 from which this function was modified to
% ZQ 01/2026
%


function [clusters, p_sig, tSums_sig] = ieeg_clusterTest1(A, alphaThresh, alphaFA, n_perm, tail)

    if nargin < 2, alphaThresh = 0.05; end
    if nargin < 3, alphaFA     = 0.05; end
    if nargin < 4, n_perm      = 1000; end
    if nargin < 5, tail        = 0;    end

    assert(abs(tail)==1 | tail==0, 'Tail can only take values {-1,0,1}!')

    % get the tstatistic sum of all clusters
    [tSums, rgs] = getTsums(A, alphaThresh, tail);

    if isempty(rgs)
        fprintf('No cluster detected at alphaThresh of %0.02f\n', alphaThresh);
        clusters = []; p_sig = [];
        return;
    end

    % sort tsums in order of decreasing absolute value (greatest differences first)
    [~, idx] = sort(abs(tSums), 'descend');
    tSums = tSums(idx);
    rgs = rgs(idx, :);

    % Permutation testing to get null distribution
    tSum_null = zeros(n_perm, 1); % tSum = 0 for permutations where no cluster above threshold is found
    
    for pp = 1:n_perm
    
        % Sign-flip permutation across trials
        signs = (rand(size(A,2),1) > 0.5)*2 - 1;  % +1 or -1
        A_perm = A .* signs';
  
        tSums_perm = getTsums(A_perm, alphaThresh, tail);
        
        tSums_perm_max = max(abs(tSums_perm));

        if ~isempty(tSums_perm_max), tSum_null(pp) = tSums_perm_max; end
    
    end

    p_val = nan(size(tSums));
    for tt = 1:length(p_val)
        p_val(tt) = (sum(tSum_null >= abs(tSums(tt))) + 1) / (n_perm + 1);
        %p_val(tt) = sum(tSum_null >= abs(tSums(tt)) | tSum_null <= -abs(tSums(tt))) / nperm; % area under null more extreme in both tails
    end
    clusters = rgs(p_val <= alphaFA, :);
    p_sig = p_val(p_val <= alphaFA);
    tSums_sig = tSums(p_val <= alphaFA);

end
%%
function [tSums, rgs] = getTsums(x, alphaThresh, tail)

    % Observed t-values (vs zero mean)
    tstat = mean(x,2) ./ std(x,0,2) * sqrt(size(x,2));

    % Suprathreshold mask
    df = size(x,2) - 1;
    switch tail
        case 1,  mask = tstat >= tinv(1-alphaThresh, df);
        case -1, mask = tstat <= tinv(alphaThresh, df);
        case 0,  mask = abs(tstat) >= tinv(1-alphaThresh/2, df);
    end

    % 1D clustering along time (contiguous suprathreshold timepoints)
    % pos_starts = find([true; diff(double(mask)) ~= 0]);
    % pos_stops = [pos_starts(2:end)-1; find(mask, 1, 'last')];
    pos_starts = find(diff([0; double(mask)]) == 1);
    pos_stops  = find(diff([double(mask); 0]) == -1);

    rgs = [pos_starts, pos_stops];

    tSums = arrayfun(@(ii) sum(tstat(pos_starts(ii):pos_stops(ii))), 1:length(pos_starts))';

end