function data_epoch = ecog_baselinesubtract(data_epoch,t_base)

% data_epoch = ecog_baselinesubtract(data_epoch,t_base)
% subtracts the mean signal during t_base in each epoch
% data_epoch = electrodes X epoch X t

% baseline correct
for k=1:size(data_epoch,1)%channels
    for m=1:size(data_epoch,2)%epochs
        x=squeeze(data_epoch(k,m,:));
        x=x-mean(x(t_base));
        data_epoch(k,m,:)=x;
    end
end
