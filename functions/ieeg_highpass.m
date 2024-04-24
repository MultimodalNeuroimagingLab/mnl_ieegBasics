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

    % define the highpass filter passband, stopband, ripple and attenuation
    passFreq    = 0.50;     % Hz
    stopFreq    = 0.05;     % Hz
    passRipple  = 3;        % dB
    stopAtten   = 30;       % dB

    % normalize the passband and stopband to the Nyquist frequency
    normPassFreq = passFreq / (srate / 2); 
    normStopFreq = stopFreq / (srate / 2);

    % calculate the minimum filter order for butterworth filter
    [filtOrder, cutFreq] = buttord(normPassFreq, normStopFreq, passRipple, stopAtten);
    filtOrder = filtOrder + rem(filtOrder, 2);

    % calculate the filter coefficients in Zero-Pole-Gain design
    [filtZeros, filtPoles, filtGains] = butter(filtOrder, cutFreq, 'high');
    [filtSos, filtOverallGain] = zp2sos(filtZeros, filtPoles, filtGains);

    % filter the data
    signal_hp = zeros(size(signal));
    for idx_channel = 1:size(signal, 2)
        
        warning('off', 'signal:filtfilt:ParseSOS');
        signal_hp(:, idx_channel) = filtfilt(filtSos, filtOverallGain, double(signal(:, idx_channel)));
        if silent == 0, fprintf(1, '.'); end
        
    end
    
    if silent == 0, fprintf(1, '] done\n'); end
    
end