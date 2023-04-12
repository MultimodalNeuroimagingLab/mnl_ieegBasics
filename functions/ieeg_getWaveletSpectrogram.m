%% Returns the Morlet (Gabor) wavelet transform (Spectrogram) for a signal (or multiple trials) V.
% Can choose between implementations: Matlab cwt() function, which sets wavelet scales based on default voices/octave. This is the recommended option.
% or kjm's custom function wavelet_pac(), which sets wavelet scales by user input frequencies
%
%   [S, f] = getWaveletSpectrogram(V, srate, frange);
%   [S, f] = getWaveletSpectrogram(V, srate, frange, method);
%       V =         txn or tx1 double, signal with t samples and n trials. ROWS are samples, COLS are trials.
%       srate =     num, sampling frequency of signal
%       frange =    1x2 or 1xm double, frequencies (Hz) to return for Spectrogram.
%                       Matlab cwt() implementation will only use the first and last values as frequency limits and
%                       return appropriate freq bins in between. See cwt documentation for more info.
%                       kjm wavelet_pac() implementation will return all frequencies given. If only 2 integer values are
%                       given as input, they are assumed to be the start and end frequencies; all integer frequencies in
%                       between will also be returned. E.g. [1, 200] -> 1:200.
%       method =    char (optional), 'cwt' (default) or 'kjm'; which implementation to use for transform.
%
%   Returns:
%       S =         axtxn or axt double, Spectrogram obtained from Morlet wavelet transform of each column in <V>, with a
%                       frequencies, t time points, and n trials. axt if only 1 trial given. Values are real units of
%                       POWER (calculated as the square of the absolute value of the convolution output).
%       f =         ax1 double, frequencies corresponding to S, in ascending order. a == m if kjm method is chosen and
%                       frange is not a length-2 integer array.
%
%   HH 2021, wrapper for kjm wavelet_pac and matlab cwt functions.
%

function [S, f] = ieeg_getWaveletSpectrogram(V, srate, frange, method)
    
    if nargin < 4, method = 'cwt'; end
    if size(V, 2) > size(V, 1), warning('V is wider than tall. Check that ROWS are samples'); end

    switch method
        case 'cwt'
            if length(frange) > 2, warning('Only using ends of frange for cwt implementation'); end
            [S, f] = getWavCwt(V, srate, frange([1, end]));
            
        case 'kjm'
            if length(frange) == 2 && all(floor(frange) == frange) % assume given range start and end
                warning('Returning all integer frequencies between %d and %d using kjm implementation', frange(1), frange(2));
                frange = frange(1):frange(2);
            end
            S = getWavKjm(V, srate, frange);
            f = frange;
        otherwise
            error("method must be 'cwt' (default) or 'kjm'");
    end
    S = squeeze(S); % squeeze 3rd dimension (trials) if one trial given
    
end


% Matlab cwt implementation for each column of V
function [S, f] = getWavCwt(V, srate, frange)
    
    for ii = 1:size(V, 2)
        [wt, f] = cwt(V(:, ii), 'amor', srate, 'FrequencyLimits', frange);
        S(:, :, ii) = wt;
    end
    S = abs(S).^2;
    
    % order from smallest to largest frequency
    f = flip(f);
    S = flip(S, 1);
    
end


% KJM custom wavelet implementation for each column of V
function S = getWavKjm(V, srate, frange)
    
    S = nan(length(frange), size(V, 1), size(V, 2));
    for ii = 1:size(V, 2)
        S(:, :, ii) = wavelet_pac(V(:, ii), srate, frange)';
    end
    S = abs(S).^2; % convert to power
    
end


% function [wsignal]=wavelet_pac(rawsignal,srate,frange)
% calculated wavelet fintered signal from frange (high:step:low)
% rawsignal=signal(:,44);
% frange=30:-1:1;
% srate=512;
% kjm
function wsignal = wavelet_pac(rawsignal, srate, frange, verbose)

    if nargin < 4, verbose = false; end

    wsignal = zeros(length(rawsignal), length(frange));
    fmax = frange(end);

    for kk = 1:length(frange)

        f = frange(kk);
        if mod(f, 10) == 0 && verbose, disp(['on ' num2str(f) ' of ' num2str(fmax)]), end

        %create wavelet
        t = 1:floor(5*srate/f);% 5 cycles
        wvlt = exp(1i*2*pi*f*(t-floor(max(t)/2))/srate).*gausswin(max(t))'; %gaussian envelope

        %calculate convolution
        tconv = conv(wvlt, rawsignal);
        tconv([1:(floor(length(wvlt)/2)-1) floor(length(tconv)-length(wvlt)/2+1):length(tconv)]) = []; %eliminate edges 
        wsignal(:, kk) = tconv;

    end
end