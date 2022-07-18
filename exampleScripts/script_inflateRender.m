
% This is an example script where we inflate the cortical surface rendering
% to better show the positions of sEEG electrodes
%
% Markus Adamek (Neurotech Center) and Dora Hermes (Mayo Clinic), 2022

%% example data location

data_path = '';

% get electrode positions matrix
electrodes_tsv_name = fullfile(data_path,'sub-02_ses-ieeg01_electrodes.tsv');
loc_info = readtable(electrodes_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% load pial and inflated surfaces
gL = gifti(fullfile(data_path,'pial.L.surf.gii'));
gR = gifti(fullfile(data_path,'pial.R.surf.gii'));
gL_infl = gifti(fullfile(data_path,'inflated.L.surf.gii'));
gR_infl = gifti(fullfile(data_path,'inflated.R.surf.gii'));

% sulcal labels
sulcal_labelsL = read_curv(fullfile(data_path,'lh.sulc'));
sulcal_labelsR = read_curv(fullfile(data_path,'rh.sulc'));

%% Snap electrodes to the closest point on the surface and then move to inflated
% This defaults to closest pial surface point within 6 mm, note that this may
% include some white matter electrodes, but since Freesurfer is not always
% accurate we are happy with this for now...
% we could also consider using the white surface or so 
xyz_inflated = ieeg_snap2inflated(loc_info,gR,gL,gR_infl,gL_infl);
xyz = [loc_info.x loc_info.y loc_info.z];

%% First plot regular and inflated to check

hemi = 'r'; % choose a hemisphere to plot

if isequal(hemi,'l')
    g = gL;
    gI = gL_infl;
    sulcal_labels = sulcal_labelsL;
elseif isequal(hemi,'r')
    g = gR;
    gI = gR_infl;
    sulcal_labels = sulcal_labelsR;
end

views_plot = {[270,0],[90,0],[0 -90],[22,-10],[-22 -10],[-90 -45]};
vv = 2; % view 1

% get electrodes on this hemisphere
electrodes_thisHemi = find(ismember(loc_info.hemisphere,upper(hemi)));

% make 5 plots from pial to inflated to check
% g --> gI
perc_infl = [0:.2:1]; % steps to show in inflation process

for kk = 1:length(perc_infl)
    
    g_plot = g;
    % linear interpolation for vertices
    g_plot.vertices = (1-perc_infl(kk))*g.vertices + perc_infl(kk)*gI.vertices; 
    v_d = [views_plot{vv}(1),views_plot{vv}(2)];
    
    % linear interpolation for electrode coordinates
    els = (1-perc_infl(kk))*xyz + perc_infl(kk)*xyz_inflated;
    
    figure
    tH = ieeg_RenderGifti(g_plot,sulcal_labels);
    ieeg_elAdd(els(electrodes_thisHemi,:),[.99 .99 .99],40) % add electrode positions
    ieeg_elAdd(els(electrodes_thisHemi,:),'k',30) % add electrode positions
    ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
    tH.FaceAlpha = 0.9; % slightly transparent
    set(gcf,'PaperPositionMode','auto')
end


