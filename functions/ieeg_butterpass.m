%  Bandpass filtering for specified band
%
%  ieeg_butterpass(signal, band, srate, silent)
% 
%      signal       = time X channels;
%      band         = [start:stop] or [start stop];
%      silent       = [optional] flag whether be non verbose
%
function [band_sig] = ieeg_butterpass(signal, band, srate, silent)
    if exist('silent', 'var') == 0,  silent = 0;     end
    assert(size(signal, 1) > size(signal, 2), 'signal must have rows = samples, columns = trials');
    num_chans = size(signal, 2);
    
    [zz, pp, ii] = butter(3, band*2/srate);
    [sos, g] = zp2sos(zz, pp, ii);

    % fband=filtfilt(bf_b, bf_a, signal(:,84)); %band pass
    band_sig=zeros(size(signal));
    for kk = 1:num_chans
        % just for nice disp:
        if silent == 0 && mod(kk,5)==0,disp(strcat(num2str(kk),'/',num2str(num_chans))),end %this is to tell us our progress as the program runs
        band_sig(:,kk) = filtfilt(sos, g, signal(:,kk)); % band pass
    end
    
end




