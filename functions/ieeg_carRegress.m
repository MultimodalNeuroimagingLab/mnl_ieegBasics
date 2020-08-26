%  Common average referencing flunction with regressing out the mean, rather
%  than subtracting out the mean
%
%  [signal] = ieeg_carRegress(signal, chans2incl, silent)
%
%      signal       = electrodes X samples
%      chans2incl   = channels to include in CAR
%      silent       = [optional] flag whether be non verbose
%
% dh - Oct 2010
% 
function [signal] = ieeg_carRegress(signal, chans2incl, silent)
    if exist('silent', 'var') == 0,  silent = 0;     end

    if size(signal,2) < size(signal,1) % signal samples X electrodes
        disp('transpose signal to be electrodes X samples')
        return
    end

    ca_signal = mean(signal(chans2incl,:),1);

    % regress off the mean signal
    for kk = 1:size(signal,1) % elecs
        if silent == 0, disp(['elec ' int2str(kk)]);    end
        if ismember(kk,chans2incl)
            [B,BINT,R] = regress(signal(kk,:)',ca_signal');
            signal(kk,:)=R';
        end
    end
    
end
