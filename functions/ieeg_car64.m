function [signalOut,spatfiltmatrix] = ieeg_car64(signal, chans2incl)

% This function performs Common Avregae Reference (CAR) filtering on a
% signal
%
% Inputs:
% signal = time x channels
% chans2inc = channel indiced to include in common average reference
%
% Outputs:
% signalOut=signal*spatfiltmatrix;
%
% Adapted from BCI2000 code from Joshua Fialkoff, clarified inputs/outputs
% and added option to exclude channels.
%
% DH 2010

if size(signal,1) < size(signal,2) % signal samples X electrodes
    disp('transpose signal to be samples X electrodes')
    return
end

assert(max(chans2incl) <= size(signal, 2), 'Attempting to include more channels than exist');
fprintf('Calculating common average reference by 64-channel blocks\n');

spatfiltmatrix = zeros(size(signal, 2));

% add to matrix by block
for ii = 1:64:size(signal, 2)
    
    idxEnd = min(size(signal, 2), ii + 63); % index of last channel in this block
    
    chans2inclBlock = chans2incl(chans2incl >= ii & chans2incl <= idxEnd);
    badChans = setdiff(ii:idxEnd, chans2inclBlock); % channels to not include
    
    numChansBlock = length(chans2inclBlock); % how many channels to include on this block
    if numChansBlock <= 8, warning('Only %d channels will be used as reference between chs %d and %d', numChansBlock, ii, idxEnd); end
    
    % output channel = 1*inputchannel - mean(chans2incl)
    spatfiltmatrix(ii:idxEnd, ii:idxEnd) = -1/numChansBlock;
    for jj = ii:idxEnd
        spatfiltmatrix(jj, jj) = 1 - 1/numChansBlock; % set diagonal elements
    end
    
    for bb = badChans
        spatfiltmatrix(bb, :) = 0;
        %spatfiltmatrix(:, bb) = 0; % no change to columns if CAR is also subtracted from the bad channels
        spatfiltmatrix(bb, bb) = 1;
    end
    
end

signalOut=signal*spatfiltmatrix;
