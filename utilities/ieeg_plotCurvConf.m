%% Plots 95% confidence interval shadings for the input curve, y
%
%   USAGE:
%       ieeg_plotCurvConf([], y);
%       ieeg_plotCurvConf(x, y);
%       h = ieeg_plotCurvConf(x, y, color, facealph)
%       [h,conf] = ieeg_plotCurvConf(x, y, color, facealph)
%           x =             1xn numeric, x values corresponding to each column in y. If empty (enter []), x is set to
%                               1 ... <number of col in y> by default
%           y =             mxn numeric, data to calculate 95% confidence interval on. Rows correspond to individual
%                               observations and columns correspond to each value in x.
%           color =         (optional) 1x3 double or char, color of confidence interval shading. Default = black ('k')
%           facealpha =     (optional) double, transparency of confidence interval shading (0 to 1). Default = 0.5.
%
%   RETURNS:
%           h =             patch object, patch for confidence interval shading, created using fill()
%           conf =          confidence interval
%
% HH 2021
%
function [h,conf] = ieeg_plotCurvConf(x, y, color, facealph, nboot)

    if nargin < 5, nboot = []; end
    if nargin < 4 || isempty(facealph), facealph = 0.5; end
    if nargin < 3 || isempty(color), color = 'k'; end
    if isempty(x), x = 1:size(y, 2); end
    x = x(:); % set vertical

    nTrials = size(y, 1);
    assert(length(x) == size(y, 2), 'x must have same length as number of columns in y');
    
    yMean = mean(y); % bootstrapped confidence interval
    
    if nboot % bootstrapped confidence interval
        assert(isa(nboot, 'double'), 'nboot input must be given as a number (of samples to bootstrap)')
        bootstat = bootstrp(nboot, @mean, y);
        conf = [prctile(bootstat, 2.5), prctile(bootstat(:, end:-1:1), 97.5)];
    else
        ySEM = std(y)/sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom
        conf = [yMean+ts(1)*ySEM, yMean(end:-1:1)+ts(2)*ySEM(end:-1:1)];
    end
    x_mirror = [x; x(end:-1:1)];
    h = fill(x_mirror, conf, color, 'EdgeColor', 'none', 'FaceAlpha', facealph);
    
end