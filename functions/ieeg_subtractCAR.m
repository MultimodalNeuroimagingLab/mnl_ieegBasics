%
%   Returns the time series data for a channel after subtracting a common average reference (CAR) from multiple channels.
%   The CAR is constructed from the first k channels with the least noise (least-squares from 0) over the designated
%   time interval.
%
%   V_new = ieeg_subtractCAR(V, V_all, exc_rows)
%   V_new = ieeg_subtractCAR(V, V_all, exc_rows, range)
%   V_new = ieeg_subtractCAR(V, V_all, exc_rows, range, p)
%
%       V       = 1xn time series data for the measured channel
%       V_all   = cxn time series data for all channels, c (including the measured channel)
%       exc_rows= indices of rows to remove from V_all, corresponding to the stimulated electrodes
%       range   = [optional] [start, end] indices of time series interval on which to calculate least squares
%                   default = [1, end]
%       p       = [optional] proportion, expressed between 0 and 1, of channels to use in the assembly of CAR.
%                   default = 0.75
%
%   Returns:
%       V_new   = 1xn time series equal to V after subtracting the CAR
%
function V_new = ieeg_subtractCAR(V, V_all, exc_rows, range, p)
    
    if nargin<5
        p = 0.75;
    end
    if nargin<4
        range = [1 length(V)];
    end
    
    V_all(exc_rows, :) = []; % remove stimulated channels
    [~, order] = sort(sum(V_all(:, range(1):range(2)).^2, 2)); % order of channels by least squares sum
    
    var_ind = 1:floor(p*(size(V_all, 1)));
    ref = mean(V_all(order(var_ind), :), 1); % mean time series across included channels
    V_new = V - ref;
end