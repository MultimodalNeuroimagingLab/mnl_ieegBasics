%  Highpass filters the data to deal with DC
%
%  [signal_hp] = ieeg_highpass(signal, srate, silent)
% 
%      signal       = time X channels;
%      srate        = sampling frequency
%      silent       = [optional] flag whether be non verbose
%
%   Returns:
%      signal_hp    = time X channel
%
%
% code adapted from Brunner lab
% DH
%
function [signal_hp] = ieeg_highpass(signal, srate, silent)

    if exist('silent', 'var') == 0,  silent = 0;     end
    
    if size(signal, 1) < size(signal, 2)
        disp('data may be channels X time, transposing matrix and returning time X channels')
        signal = signal';
    end

    parameters.SamplingRate.NumericValue = srate;

    % --- highpass filter --- 
    param.highpass.Wp = 0.50; % Hz
    param.highpass.Ws = 0.05; % Hz
    param.highpass.Rp = 3;    % dB
    param.highpass.Rs = 30;   % dB

    % define passband, stopband and attenuation
    highpass{1}.Wp = param.highpass.Wp / (parameters.SamplingRate.NumericValue / 2); 
    highpass{1}.Ws = param.highpass.Ws / (parameters.SamplingRate.NumericValue / 2);
    highpass{1}.Rp = param.highpass.Rp; 
    highpass{1}.Rs = param.highpass.Rs;

    % calculate the minimum filter order for butterworth filter
    [highpass{1}.n, highpass{1}.Wn] = buttord(highpass{1}.Wp, highpass{1}.Ws, highpass{1}.Rp, highpass{1}.Rs);
    highpass{1}.n = highpass{1}.n + rem(highpass{1}.n, 2);

    % calculate the filter coefficients in Zero-Pole-Gain design
    [highpass{1}.z, highpass{1}.p,highpass{1}.k] = butter(highpass{1}.n, highpass{1}.Wn, 'high');
    [highpass{1}.sos, highpass{1}.g] = zp2sos(highpass{1}.z, highpass{1}.p, highpass{1}.k);
    highpass{1}.h = dfilt.df2sos(highpass{1}.sos, highpass{1}.g);

    signal_hp = zeros(size(signal));
    for idx_channel = 1:size(signal, 2)
        warning('off', 'signal:filtfilt:ParseSOS');
        signal_hp(:, idx_channel) = single(filtfilt(highpass{1}.sos, highpass{1}.g, double(signal(:, idx_channel))));
        if silent == 0, fprintf(1, '.'); end
    end
    
    if silent == 0, fprintf(1, '] done\n'); end
    
end