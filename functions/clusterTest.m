%% This function compares time series data for trials across two experimenetal conditions.
%   A nonparametric test to determine whether clusters of time points exist on which there is a significant difference in means between
%   two experimental conditions. In brief, clusters of time points that are significantly different between the two conditions are detected first using
%   2-tailed t-tests on every sample. The sum of the t-statistic in each cluster is then used as a test statistic and compared to a null distribution generated
%   by permutation testing.
%
%   Please cite: Maris and Oostenveld (2007):
%       Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-and MEG-data. Journal of neuroscience methods, 164(1), 177-190.
%
%   The above paper also gives credit to Bullmore et al. (1999) for first conceiving the test statistic (tSum), where it was called the "Cluster mass test":
%       Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S., Taylor, E., & Brammer, M. J. (1999). Global, voxel, and cluster tests,
%           by theory and permutation, for a difference between two groups of structural MR images of the brain. IEEE transactions on medical imaging, 18(1), 32-42.
%
%   For a recent publication that implements this test, see:
%       Wang, Y., Luo, L., Chen, G., Luan, G., Wang, X., Wang, Q., & Fang, F. (2023). Rapid processing of invisible fearful faces in the human amygdala.
%           Journal of Neuroscience, 43(8), 1405-1413.
%
%   USAGE:
%       clusters = clusterTest(x1, x2);                                         % minimal inputs/outputs
%       [clusters, pSig] = clusterTest(x1, x2, alphaThresh, alphaPer, nperm);   % full inputs/outputs
%           x1 =            t x n num. n Trials corresponding to first experimental condition, with t time points to be compared
%           x2 =            t x m num. m Trials corresponding to second experimental condition.
%           alphaThresh =   (optional) num, 0 < alphaThresh < 1. Alpha of the 2-tailed t-tests to detect samples with significant (uncorrected) difference in means between x1 and x2.
%                               These samples are clustered on and then subjected to the non-parametric testing. alphaThresh does NOT affect the false alarm (FA)
%                               rate of the test, but it does affect the sensitivity. Higher alphaThresh detects weaker/longer-lasting effects. Default = 0.05.
%           alphaFA =       (optional) num, 0 < alphaFA < 1. Alpha (false alarm (FA) rate) of the overall test. Sets the critical values of the two-tailed permutation null distribution.
%                               This alpha controls the FA rate for ALL clusters tested (see section 3.1.1 last paragraph of Maris & Oostenveld, 2007). Default = 0.05.
%           nperm =         (optional) num > 0. Number of permutations in null distribution. Default = 1000.
%
%   RETURNS:
%       clusters =          b x 2 num. Each row is the [start, stop] indices in x1, x2 of a significant cluster. There are b significant clusters detected. No significant
%                               clusters if empty.
%       pSig =              b x 1 num. The permutation p-value for each significant cluster in clusters. Every element of pSig <= alphaFA.
%           
%   To-consider:
%       Maris & Oostenveld, 2007 discusses optional limits on the minimum cluster size. This is not currently implemented and the user can filter for minimum
%       cluster size on the function output. Note however, that the minimum cluser size can be specified as a parameter on the test statistic and permutations.
%       This could affect the null distribution. E.g. 1, a permutation where only significant clusters < minimum size would have a tSum value of 0. E.g. 2, if a
%       permutation has a max tSum corresponding to a cluster smaller than minimum size and another cluster larger than the minimum size, then setting a minimum
%       size constraint would change the tSum for that perm from the original to the lower score > minimum size.
%
% HH 08/2023
%
function [clusters, pSig] = clusterTest(x1, x2, alphaThresh, alphaFA, nperm)
    
    assert(size(x1, 1) == size(x2, 1), 'x1 and x2 must have the same number of rows (samples)');

    if nargin < 5, nperm = 1000; end
    if nargin < 4 || isempty(alphaFA), alphaFA = 0.05; end
    if nargin < 3 || isempty(alphaThresh), alphaThresh = 0.05; end

    % number of trials in each condition
    n = size(x1, 2); m = size(x2, 2);

    % Get t-statistic critical values for desired alpha1 (to detect which time points are significant)
    df = size(x1, 2) + size(x2, 2) - 2; % 2-sample ttest df
    crits = tinv([alphaThresh/2, 1 - alphaThresh/2], df); % alpha divided by 2 to reflect 2-tailed critical values on t-distribution.

    % calculate t-statistic at each time point
    tstat = calcTstat(x1, x2);
    
    % get the tstatistic sum of all clusters
    [~, ~, tSums, rgs] = getTsums(tstat, crits);

    if isempty(rgs)
        fprintf('No cluster detected at alpha1 of %0.02f\n', alphaThresh);
        clusters = []; pSig = [];
        return;
    end

    % sort tsums in order of decreasing absolute value (greatest differences first)
    [tSums, idx] = sort(abs(tSums), 'descend');
    rgs = rgs(idx, :);
    
    % Permutation testing to get null distribution
    tSumNull = zeros(nperm, 1); % tSum = 0 for permutations where no cluster above threshold is found
    for ii = 1:nperm

        % permute trials between the two conditions
        xPooled = [x1, x2];
        xPooled = xPooled(:, randperm(n + m));
        x1Perm = xPooled(:, 1:n); x2Perm = xPooled(:, n+1:end);

        % get tSumMax as statistic
        tstatPerm = calcTstat(x1Perm, x2Perm);
        tSumMaxPerm = getTsums(tstatPerm, crits);
        if ~isempty(tSumMaxPerm), tSumNull(ii) = tSumMaxPerm; end

    end

    % empirical 2-tailed p-value for each tSum
    p = nan(size(tSums));
    for ii = 1:length(p)
        p(ii) = sum(tSumNull >= abs(tSums(ii)) | tSumNull <= -abs(tSums(ii))) / nperm; % area under null more extreme in both tails
    end
    clusters = rgs(p <= alphaFA, :);
    pSig = p(p <= alphaFA);

end


%% Helper functions

% Input: a vector of t-statistics and two-tailed critical values on the t-distribution, calculated a prior with the correct degrees of freedom.
% This function calculates the sum of t-statistics (tSum) for contiguous time points with significant tstat (clusters), separately for samples with positive and negative tstat.
% Returns the tSum corresponding to the maximum cluster (tSumMax) and the cluster [start, stop]. Also returns all tSums and [start, stop] for all clusters
function [tSumMax, rgMax, tSums, rgs] = getTsums(tstat, crits)
    
    % find significant time points, for positive and negative separately
    hNeg = tstat < crits(1);
    hPos = tstat > crits(2);

    % starts and stops of the positive t-value clusters
    startsPos = find(diff([0; hPos]) == 1);
    stopsPos = find(diff([hPos; 0]) == -1);

    % calculate the tstatistic sum of each positive cluster
    tSumsPos = arrayfun(@(ii) sum(tstat(startsPos(ii):stopsPos(ii))), 1:length(startsPos))';
    
    % starts and stops of the negative t-value clusters
    startsNeg = find(diff([0; hNeg]) == 1);
    stopsNeg = find(diff([hNeg; 0]) == -1);

    % calculate the tstatistic sum of each negative cluster
    tSumsNeg = arrayfun(@(ii) sum(tstat(startsNeg(ii):stopsNeg(ii))), 1:length(startsNeg))';

    % combine positive and negative clusters into one ranges variable
    rgs = [startsPos, stopsPos; startsNeg, stopsNeg];
    tSums = [tSumsPos; tSumsNeg];

    % find cluster with maximum absolute value cluster sum
    [~, idx] = max(abs(tSums));
    rgMax = rgs(idx, :);
    tSumMax = tSums(idx);
    
end

% Calculates the t-statistic for each row of input matrices x1 and x2 (e.g. at each time point of column-wise trials
% ~5x faster than matlab ttest2 function
% Assumes equal variance between x1 and x2, which is satisfied if these are similar ERP/time series trials. This is default option for Matlab ttest2
% See documentation for Matlab ttest2 for the formulas: https://www.mathworks.com/help/stats/ttest2.html#btqjrsu-stats
function tstat = calcTstat(x1, x2)

    % Number of trials in x1 and x2
    n = size(x1, 2); m = size(x2, 2);

    % Pooled SD (assumes equal variances)
    s = sqrt(((n - 1)*var(x1, 0, 2) + (m - 1)*var(x2, 0, 2)) / (n + m - 2));

    tstat = (mean(x1, 2) - mean(x2, 2)) ./ sqrt(s.^2/n + s.^2/m);

end