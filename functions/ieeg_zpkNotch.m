%   Passes each row of the input data through Butterworth notch filters at the specified freq + harmonics
%   Uses zero/pole/gain outputs of butter for better precision over b/a coefficient outputs
%
%   V = zpk_notch(V, srate, notch_freq);
%   V = zpk_notch(V, srate, notch_freq, butterOrd, halfWin, nharms);
%   V = zpk_notch(V, srate, notch_freq, [], [], nharms);
%       V =             n x timepoints data, where each row is a separate trial/channel
%       srate =         sample rate of data, in Hz
%       notch_freq =    fundamental frequency of notch filter, in Hz (e.g. 60 or 50)
%       butterOrd =     (optional) order of Butterworth filter to use. Default = 3
%       halfWidth =     (optional) notch stopband is from notch_freq - halfWidth to notch_freq + halfWidth (Hz).
%                                   Default = 1.
%       nharms =        (optional) number of harmonics to filter at (including fundamental frequency). Default = 3
%
%   HH 2021
%
function V = ieeg_zpkNotch(V, srate, notch_freq, butterOrd, halfWidth, nharms)

    if nargin < 6, nharms = 3; end
    if nargin < 5 || isempty(halfWidth), halfWidth = 1; end
    if nargin < 4 || isempty(butterOrd), butterOrd = 3; end

    for i=1:nharms
        [z, p, k] = butter(butterOrd, [notch_freq*i - halfWidth, notch_freq*i + halfWidth]/(srate/2), 'stop');
        [sos, g] = zp2sos(z, p, k);
        
        for j=1:size(V, 1)
            V(j, :) = filtfilt(sos, g, V(j, :));
        end
    end
end