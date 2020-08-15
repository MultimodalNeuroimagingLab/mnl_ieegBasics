function [data_out,ref,pc1]=ieeg_pcref(data,chans2incl)
% function data_reref=pc_ref(data)
% this function will subtract out the 1st principal component of the an
% array of timeseries data (e.g. SEEG, ECoG). For equal impedances & 
% variances, it should reduce to the common average re-reference.
% 
% inputs:
% data: input data timeseries (channel x time)
% chan2incl: channel indiced to include in common average reference
% 
% outputs:
% data_reref: reconstructed data with 1st principal component removed (time x channel)
% ref: reconstructed reference (time x 1)
% pc1: first principal component (to see relative contribution of reference 
%      to each channel - can infer relative contact impedances from this)
% 
% kjm 05/2020     
% 
% dh, 08/2020, option to include good channels 

%% principal components analysis of timeseries

    tt=cov(data(chans2incl,:)'); %covariance matrix of data timeseries

    % do a principal component analysis
    [pc_vecs,pc_vals]=eig(tt); 
    [pc_vals,v_inds]=sort(sort(sum(pc_vals)),'descend'); pc_vecs=pc_vecs(:,v_inds); %reshape properly    
    
    % get timeseries of principal components
    pc_ts=data(chans2incl,:)'*pc_vecs; 


%% reconstruct data without 1st principal component
    pc_tmp=pc_vecs; pc_tmp(:,1)=0;
    data_reref=pc_ts*pinv(pc_tmp);

    % put included channels back in data:
    data_out = data;
    data_out(chans2incl,:) = data_reref';
    
%% calculate reconstructed reference and make a separate variable for 1st principal component
    ref=pc_ts(:,1)/(sqrt(size(data(chans2incl,:),1)));
    pc1=pc_vecs(:,1);

