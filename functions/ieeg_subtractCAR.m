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
%   Example:
%       
%       % set up time interval = [-1s, 4s]
%       epoch_length = 5; epoch_prestim_length = 1; % in s
%       tt = (1:epoch_length*srate)/srate - epoch_prestim_length;
%
%       [~, V] = readMef3(fileName, [], {'LMS2'}, 'samples', [1042583 1052823]); % Epoch from channel LMS2 matching length of tt
%       [~, V_all] = readMef3(fileName, [], [], 'samples', [1042583 1052823]); % load same epoch for all channels
%
%       [~, start_idx] = find(tt>0.1, 1, 'first'); % use t=[100ms, 2s] to calculate least-squares ordering of channels for CAR
%       [~, end_idx] = find(tt<2, 1, 'last');
%
%       exc_rows = [3, 4]; % stimulated pair matches channels 3, 4
%       V = ieeg_subtractCAR(V, V_all, exc_rows, [start_idx, end_idx]); % adjust V by subtracting CAR
%
%       Harvey Huang 2020
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