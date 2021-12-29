% Plots multiple curves spaced out vertically by yspace
%
%   ys = plotTrials(tt, data, yspace, cm, Name, Value)
%       tt =        1xt time vector
%       data =      txn, rows are samples and columns are trials
%       yspace =    1x1, space between each trial
%       cm =        mx3, color map (optional)
%       Optional Name, Value pairs to pass to built-in plot function
%
%   Returns
%       ys =        0-value for each curve plotted
%
function ys = plotTrials(tt, data, yspace, labels, cm, varargin)
    
    if ~exist('labels', 'var') || isempty(labels)
        labels = repmat({''}, size(data, 2), 1); % no labels
    end
    if ~exist('cm', 'var') || isempty(cm)
        cm = get(0, 'DefaultAxesColorOrder'); 
    end
    assert(length(labels) == size(data, 2), 'number of labels needs to match number of columns in data');
    
    numTrials = size(data, 2);
    tt = tt(:);
    cm = repmat(cm, [ceil(size(data, 2)/size(cm, 1)), 1]); % tile cmap if not long enough
        
    ys = yspace*(0:-1:-(numTrials - 1));
    hold on
    for ii = 1:numTrials
        plot(tt, data(:, ii) + ys(ii), 'Color', cm(ii, :), varargin{:});
        yline(ys(ii), 'Color', 0.5*[1 1 1]);
    end
    hold off
    
    ylim([ys(end)-2*yspace, ys(1)+2*yspace]); % basic aesthetics
    set(gca, 'YTick', (-yspace)*(length(labels)-1:-1:0), 'YTickLabel', flip(labels), 'FontSize', 10);
    
end