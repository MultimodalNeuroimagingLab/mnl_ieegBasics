function [signaldata, out] = ccep_CAR64blocks_percent(signaldata,ttt,good_channels,perc_channels,car_timeint,notchOpt)
% 
% function [data_out,out] = ccep_CAR64blocks_percent(data,tt,good_channels,perc_channels,car_timeint);
% 
%   function applied common average referencing on CCEP data, excludes
%   channels with high variance from 10-100 ms after stimulation onset
%   these are channels with large evoked responses
%
%   function assumes that noise is shared between blocks of 64 channels.
%
%
% inputs
%   signaldata: channels X time X epochs 
%   ttt: time vector in seconds 
%   good_channels: channels numbers to include, bad channels that are noisy
%           should be excluded, stim channels should also be excluded
%   perc_channels: proportion of channels with lowest variance to include
%           in the common average: 0.1 for 10%
%   car_timeint: time interval to calculate variance, in seconds, default
%           is 15-500ms: [0.015 0.500]
%   notchOpt: structure to indicate notch of not, empty assumes no notch
%           notchOpt.do = 0; % 0: no notch, 1: do notch
%           notchOpt.freq = 60; % 50 or 60
%           notchOpt.srate = 4800; % sampling frequency
%
% outputs:
%   signaldata: channels X time X epochs 
% 
% DH and HH Multimodal Neuroimaging Lab, Mayo Clinic, 2020

if isempty(perc_channels)
    perc_channels = .30;
end

if isempty(car_timeint)
    car_timeint = [0.015 0.500];
end

if ~isempty(notchOpt) 
    if notchOpt.do==1 % do notch before getting channels used in car
        disp('notch filtering data to decide which channels go in Common Average')
        signalnotch = zeros(size(signaldata));
        for kk = 1:size(signaldata,1) % channels
            signalnotch(kk,:,:) = ieeg_notch(squeeze(signaldata(kk,:,:)), notchOpt.srate, notchOpt.freq);
        end
    end
else
    notchOpt.do = 0;
end


out = [];

% create blocks to use for 64 channel groups
set_nrs = 1:ceil(size(signaldata,1)/64);
for ss = set_nrs
    set_inds = (set_nrs(ss)*64-63):set_nrs(ss)*64; % 1:64 65:128 etc...
    set_inds = set_inds(set_inds<=max(good_channels)); % make sure to stay below last good channel
    out(ss).channels_set = set_inds;
end

% split into blocks and get car channels 
for ss = set_nrs
    
    % total set of channels that are good & in the set:
    these_channel_nrs = intersect(good_channels,out(ss).channels_set);
    
    % get data for car in this set
    if notchOpt.do==1 % only use notch filtered channels
        these_data = signalnotch(these_channel_nrs,ttt>car_timeint(1) & ttt<car_timeint(2),:);
    elseif notchOpt.do==0 % use data, no notched data
        these_data = signaldata(these_channel_nrs,ttt>car_timeint(1) & ttt<car_timeint(2),:);
    end
    
    % concatinate trials (these_data will be channels X time*trials
    these_data_cat = reshape(these_data,size(these_data,1),size(these_data,2)*size(these_data,3));

    % calculate variance for each channel across time
    chan_var = var(these_data_cat,[],2); 

    % set a threshold for which channels to reject based on variance
    var_th = quantile(chan_var,perc_channels);
       
    % include channels with smallest variance (out of good channels in set)
    chans_incl = setdiff(1:length(these_channel_nrs),find(chan_var>var_th));
    
    % original channels numbers to use for car
    out(ss).car_channels = these_channel_nrs(chans_incl);
end


% now do CAR
for ss = set_nrs % run across blocks 
    
    % average across car channels
    out(ss).car_data = mean(signaldata(out(ss).car_channels,:,:),1);

    % subtract from all other data
    signaldata(out(ss).channels_set,:,:) = signaldata(out(ss).channels_set,:,:) - repmat(out(ss).car_data,length(out(ss).channels_set),1,1);
    
end



