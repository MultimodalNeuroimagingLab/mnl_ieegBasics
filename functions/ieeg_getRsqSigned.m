%% Computes the R-squared value for two one-dimensional distributions x1 and x2
%   R-squared in this context is the % variance explained by using the group means of x1 and x2 (compared to the null
%   model of using the global mean)
%
%   Also returns P-value and statistics of fit under F-test
%   Optionally makes plot to visualize difference in means between x1 and x2
%   Modified by HH from legacy DH/KJM code, 04/2021. Edited on 9/17 to include stats and plot
%   
% USAGE:
%   rsq = ieeg_getRsqSigned(x1, x2);
%   [rsq, p, stats] = ieeg_getRsqSigned(x1, x2, plotSamps);
%       x1 =            nx1 double. Vector of observations corresponding to condition 1
%       x2 =            mx1 double. Vector of observations corresponding to condition 2
%       plotSamps       (optional) boolean. If true, plots a visualization of the difference between x1 and x2.
%                           Default = false
%
% RETURNS:
%       rsq =           double. Signed R-squared value, [-1, 1]. % variance explained by group means of
%                           x1 and x2. Negative values indicate that mean(x1) < mean(x2).
%       p =             double. P-value, [0, 1]. P-value of fit, calculated using an F-distribution with 1 degree of
%                           freedom numerator and N-2 degrees of freedom denominator.
%       stats =         struct. Struct containing F-statistic and df corresponding to p-value.
%
% EXAMPLE:
%   x1 = randn(1000, 1);
%   x2 = 1 + randn(1000, 1);
%   [rsq, p, stats] = ieeg_getRsqSigned(x1, x2, true);
%
%   HH 2021
%
function [rsq, p, stats] = ieeg_getRsqSigned(x1, x2, plotSamps)

    if nargin < 3, plotSamps = false; end

    n1 = length(x1);
    n2 = length(x2);
    
    x1 = double(x1(:)); x2 = double(x2(:)); % ensure as double column vectors

    SStot = (n1 + n2)*var([x1; x2], 1); % total sum of squares (var(_, 1) normalized by N)
    SSres = n1*var(x1, 1) + n2*var(x2, 1); % residual sum of squares using respective means

    rsq = 1 - (SSres/SStot);

    if mean(x1) < mean(x2), rsq = -rsq; end % invert sign if mean of second vector larger than first vector
    
    %% Calculate p-value from F-test
    
    F = (SStot - SSres)/(SSres/(n1+n2-2)); % F statistic using 1 df in numerator and N-2 df in denominator
    p = 1 - fcdf(F, 1, n1+n2-2); % p-value from CDF
    
    stats.Fstat = F;
    stats.df = [1, n1+n2-2];
    
    if ~plotSamps, return; end
    
    %% Plot x1, x2, and show their means relative to overall mean
    
    p1 = plot(1:n1, x1, 'bo');
    hold on
    plot([1, n1], [mean(x1), mean(x1)], 'b-', 'LineWidth', 1);
    p2 = plot(n1+1:n1+n2, x2, 'ro');
    plot([n1+1, n1+n2], [mean(x2), mean(x2)], 'r-', 'LineWidth', 1);
    yline(mean([x1; x2]), 'k-', 'LineWidth', 1); % mean of all points (null model)
    
    legend([p1, p2], {'x1', 'x2'});
    xticks([]);
    xlabel('samples');
    
end