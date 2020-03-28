function data_epoch = ieeg_baselinesubtract(data_epoch,t_base)

% data_epoch = ecog_baselinesubtract(data_epoch,t_base)
% subtracts the mean signal during t_base in each epoch
% data_epoch = electrodes X epoch X t
% dhermes, multimodal neurimaging lab, 2020

% baseline correct
for kk = 1:size(data_epoch,1)%channels
    for mm = 1:size(data_epoch,2)%epochs
        x = squeeze(data_epoch(kk,mm,:));
        x = x-mean(x(t_base));
        data_epoch(kk,mm,:) = x;
    end
end
