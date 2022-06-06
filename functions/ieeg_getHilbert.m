%% This function returns the instantaneous power & phase at all samples of a signal, for an input band
%   The instantaneous amplitude of a signal at any given time is its envelope, as calculated by the magnitude of the
%   analytical Hilbert transform. Similarly, its instantaneous phase is the angle of the Hilbert transform.
%   Bandpass is performed using a 3rd order Butter filter forward and backwards (filtfilt, total order = 6), per
%   ieeg_butterpass.
%
%   power = getHilbert(sigdata, band, srate);
%   [power, phase] = getHilbert(sigdata, band, srate, powerType);
%       sigdata =       txn double, input signal with rows = samples & cols = individual trials
%       band =          1x2 num, frequency band for bandpass, in Hz.
%       srate =         num, sampling frequency of signal
%       powerType =     (optional) char array, 'amplitude', 'power', or 'logpower'. The type of <power> output to
%                           return. 'amplitude' returns the magnitude of the Hilbert transform. 'power' returns the
%                           amplitude squared, and 'logpower' (default) returns the log10 transform of power.
%       verbose =       (optional). Whether to print output
%
%   Returns:
%       power =         txn double, power (see inputs) of input signal at all samples, corresponding to the input structure
%       phase =         txn double, angle of the Hilbert Transform, corresponding to instantaneous phase of signal.
%
%   HH 2021
%       
function [power, phase] = ieeg_getHilbert(sigdata, bands, srate, powerType, verbose)
    
    if nargin < 5, verbose = false; end
    if nargin < 4 || isempty(powerType), powerType = 'logpower'; end

    assert(size(sigdata, 1) > size(sigdata, 2), 'Data must have rows=samples and cols=trials');
    
    if numel(bands) == 2, bands = bands(:)'; end % ensure row vector for 1 band
    assert(size(bands, 2) == 2, 'bands must be given as 2-column array [bandstart bandstop]');
    
    powers = zeros(size(sigdata, 1), size(sigdata, 2), size(bands, 1));
    phases = zeros(size(sigdata, 1), size(sigdata, 2), size(bands, 1));
    for ii = 1:size(bands, 1)
        
        if verbose, fprintf('Calculating amplitude between %d and %d Hz\n', bands(ii, 1), bands(ii, 2)); end
        
        sigdata_pass = ieeg_butterpass(sigdata, bands(ii, :), srate, 1);
        %sigdata_pass = ieeg_butterpass2(sigdata, bands(ii, :), srate);
        sigHilb = hilbert(sigdata_pass);

        powers(:, :, ii) = abs(sigHilb).^2;
        phases(:, :, ii) = angle(sigHilb);
    end
    
    if verbose, fprintf('Returning (combined) %s and phase\n', powerType); end
    power = geomean(powers, 3); % geometric mean of powers across bands
    phase = mean(phases, 3); % arithmetic mean of phases (needs verification)

    switch lower(powerType)
        case 'amplitude'
            power = sqrt(power); % sqrt of mean power across bands
            
        case 'power'
            return
            
        case 'logpower'
            power = log10(power);

        otherwise
            error('powerType must be given as "amplitude", "power", or "logpower".');
    end
    
end