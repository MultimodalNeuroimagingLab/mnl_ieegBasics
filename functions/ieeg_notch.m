%% This function applies a notch (bandstop) filter to the input signal, at the fundamental frequency requested + 2 harmonics
%
%   data = ieeg_notch(data, srate, notch_freq);
%   data = ieeg_notch(data, srate, notch_freq, stopWidth);
%   data = ieeg_notch(data, srate, notch_freq, stopWidth, verbose);
%       data =          txn double. Input signal with t time points and n channels/trials. Filter will be applied along each column
%       srate =         double. Sampling frequency of signal, in Hz.
%       notch_freq =    double. Center of fundamental frequency to bandstop, in Hz. (e.g. 60Hz line noise in North America)
%       stopWidth =     (optional) double. Width of bandstop filter (between the two -3dB frequencies) in Hz. Default = 2, so bandstop at 60 Hz has -3dB
%                           frequencies at 59 Hz and 61 Hz.
%       verbose =       (optional) boolean. If true, will print info about filter. Default = false.
%
%   HH 02/2022 to replace previous ieeg_notch.m in mnl_ieegBasics.
%
function data = ieeg_notch(data, srate, notch_freq, stopWidth, verbose)

    if nargin < 5, verbose = false; end
    if nargin < 4 || isempty(stopWidth), stopWidth = 2; end

    for ii = 1:3
        fNotch = ii*notch_freq;
        fHP1 = fNotch - stopWidth/2; % -3db frequency 1
        fHP2 = fNotch + stopWidth/2; % -3db frequency 2. stopWidth is fHP2 - fHP1
        dNotch = designfilt('bandstopiir', 'FilterOrder', 4, ...
                            'DesignMethod', 'butter', ...
                            'HalfPowerFrequency1', fHP1, ... 
                            'HalfPowerFrequency2', fHP2, ...
                            'SampleRate', srate);

       if verbose, fprintf('Applying filtfilt with 4th order Butterworth notch, -3dB at %0.1f Hz and %0.1f Hz.\n', fHP1, fHP2); end
       data = filtfilt(dNotch, data);
    end
    
end