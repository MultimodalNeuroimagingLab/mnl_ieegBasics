
function [bip_loc_info,bip_locs,hemisphere] = ieeg_Els2BipLocs(loc_info,bipolarChans)
%
% [bip_loc_info,bip_locs,hemisphere] = ieeg_Els2BipLocs(loc_info,{'LA8-LA9'})
%
% Inputs:
%   loc_info is the BIDS electrodes table, for example:
%       electrodes_tsv_name = ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv'];
%       loc_info = readtable(electrodes_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
%       Needs columns: name, x, y, z, hemisphere and Destrieux label
%   bipolarChans = bx1 cell. channel pair names take the form 'ch1-ch2' to represent which pair of channels
%       are used for the derivation e.g. {'LA1-LA2'}
%
% Outputs:
%
%   bip_loc_info: table with average location between 2 electrodes
%   bip_locs: x,y,z of the two electrodes
%   hemisphere: L or R
%
% DH, 2023

% get bipolar positions that match data
bip_locs = NaN(length(bipolarChans),3,2); % initialize
hemisphere = cell(length(bipolarChans),1);
Destrieux_label = NaN(length(bipolarChans),1);
Destrieux_label2 = NaN(length(bipolarChans),1);
for kk = 1:length(bipolarChans)
    these_els = split(bipolarChans{kk},'-');
    % check length of split - should be 2 electrodes, assuming they do not
    % have a '-' in it
    if length(these_els)==2
        el_nr1 = find(strcmpi(loc_info.name,these_els{1}));
        el_nr2 = find(strcmpi(loc_info.name,these_els{2}));
        if ~isempty(el_nr1) && ~isempty(el_nr2)
            pos12 = [loc_info.x([el_nr1 el_nr2]) loc_info.y([el_nr1 el_nr2]) loc_info.z([el_nr1 el_nr2])];
            bip_locs(kk,:,:) = pos12';

            hemisphere{kk} = loc_info.hemisphere{el_nr1};
            
            if ismember('Destrieux_label',loc_info.Properties.VariableNames)
                Destrieux_label(kk) = loc_info.Destrieux_label(el_nr1);
                Destrieux_label2(kk) = loc_info.Destrieux_label(el_nr2);
            end
        end
    else
        disp([bipolarChans{kk} ' not included - missing 2 electrodes']) % checking whether there are 2 electrodes
    end
end

for kk = 1:length(bipolarChans)
    if isempty(hemisphere{kk})
        hemisphere{kk} = 'n/a';
    end
end

% now we have locations of each electrode in the bipolar pair: bip_locs
% we can average to get the center for ease of use and plotting:
bip_centerloc = mean(bip_locs,3);

% create new loc_info structure for bipolar electrodes
name = bipolarChans;
x = bip_centerloc(:,1);
y = bip_centerloc(:,2);
z = bip_centerloc(:,3);

bip_loc_info = table(name,x,y,z,hemisphere,Destrieux_label,Destrieux_label2);
