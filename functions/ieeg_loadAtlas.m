function [vert_label, area_label, cmap] = loadAtlas(freesurferPath, hemi, atlasName)

% function to load Brain Atlas
% 
% input:
%     freesurferPath: local path to subject's/fsaverage freeSurfer directory
%     hemi (l/r): left/right hemisphere 
%     atlasName: Wang/BensonV/BensonE/visf/HCP/mix/NSDfloc
%     
% output:
%     vert_label: vertex labels
%     area_label: surface labels corresponding to freesurfer vertices
%     cmap: colormap for different areas
%     
%
% Example usage:
%   [vert_label, area_label, cmap] = loadAtlas(freesurferPath, 'r', 'BensonV')
% 
% ZQ 2023
%
% updated on 15 Nov 2024: (input) atlasName now case insensitive, minor fix for mix atlas (V3d) discontinuity - lines 294-295
% updated on 17 Oct 2024: included NSD-floc atlas, updated cmap for mix atlas
% updated on 17 Sep 2024: included mix atlas, updated visf to account for L/R discrepancy


if ~exist(freesurferPath, "dir")
    error('Directory not found');
end
    
assert( hemi == 'l' | hemi == 'r', 'hemisphere either l or r');

if isempty(atlasName)
    error('Choose atlas (Wang/BensonV/BensonE/visf/HCP/mix/NSDfloc)');
else
    ky = lower(atlasName);
end

switch ky

%-------------------------------
    case 'wang'
        
    % Surface labels - Wang visual ROI
    % L Wang L, et.al. Probabilistic Maps of Visual Topography in Human Cortex. Cereb Cortex (2015) 

    if ~isfile(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.wang15_mplbl.mgz']))
        error('Wang labels not found!');
    end

    surface_labels = MRIread(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.wang15_mplbl.mgz']));
    vert_label = surface_labels.vol(:);

    area_label = { 'V1v', 'V1d', ...
                   'V2v', 'V2d', ...
                   'V3v', 'V3d', ...
                   'hV4', ...
                   'VO1', 'VO2', ...
                   'PHC1','PHC2',...    % parahippocampal cortex
                   'TO2', 'TO1', ... 
                   'LO2', 'LO1', ...
                   'V3B', 'V3A', ...
                   'IPS0','IPS1','IPS2','IPS3','IPS4','IPS5', ...
                   'SPL1','FEF'};       % superior parietal lobule

    cmap = [...
        0.9438    0.3910    0.2668
        0.1974    0.5129    0.7403
        0.5978    0.8408    0.6445
        0.3686    0.3098    0.6353
        0.9955    0.8227    0.4828
        0.8417    0.9409    0.6096
        0.6196    0.0039    0.2588
        0.8000    0.8000    0.4000
        0.3539    0.7295    0.6562
        0.9877    0.6154    0.3391
        0.8206    0.2239    0.3094
        0.9000    0.6000    1.0000
        0.3438    0.3910    0.3668
        0.6974    0.5129    0.8403
        0.1978    0.8408    0.7445
        0.3955    0.8227    0.5828
        0.9686    0.3098    0.7353
        0.2417    0.9409    0.7096
        0.2196    0.0039    0.7588
        0.2000    0.8000    0.5000
        0.9539    0.7295    0.7562
        0.2877    0.6154    0.4391
        0.2206    0.2239    0.4094
        0.9000    1.0000    0.7000
        0.2000    0.5000    0.2000];

%-------------------------------
    
    case 'bensonv'

    % Surface labels - Benson visual ROI
    % NC Benson & J Winawer. Bayesian analysis of retinotopic maps. eLife (2018) 

    if ~isfile(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.benson14_varea.mgz']))
        error('Benson labels not found!');
    end

    surface_labels = MRIread(fullfile(freesurferPath, 'surf',...
        [hemi 'h.benson14_varea.mgz']));
    vert_label = surface_labels.vol(:);

    area_label = { 'V1', ...            % primary visual  
                   'V2', ...
                   'V3', ...
                   'hV4', ...           % ventral
                   'VO1', 'VO2', ...
                   'LO1', 'LO2', ...    % lateral
                   'TO1', 'TO2', ...    % temporal
                   'V3b', 'V3a'};       % dorsal

    cmap   =  [255   0   0; 
               255 128   0; 
               255 255   0; 
                 0 128 255;
                 0  76 153; 153  76   0;
                 0 102  51; 153 153   0;  
               255  51 153; 102   0 204;
                 0 255 255;   0 255   0]./255;

 %-------------------------------            
 
    case 'visf'

    % Surface labels - Probabilistic functional O-T visual (visf) atlas - KGS
    % M Rosenke, et al. A Probabilistic Functional Atlas of Human Occipito-Temporal Visual Cortex. Cerebral Cortex (2021)
    % only for fsaverage

    if ~isfile(fullfile(freesurferPath, 'label',...
                    [hemi 'h.visfAtlas.annot']))
        error('visf labels not found!');
    end

    [~, vert_label, temp] = read_annotation(fullfile(freesurferPath, 'label',...
                            [hemi 'h.visfAtlas.annot']));

    % Labels are numbered with unique ID (nothing from 1-18 so converting for rendering)
    % Also accounting for slight differences in left vs right hemi

    if hemi == 'l'
   
        for ii = 1 : 11   
            vert_label( vert_label == temp.table( ii, 5)) = ii - 1;
        end

        for ii = 12 : 18  
            vert_label( vert_label == temp.table( ii, 5)) = ii;
        end

    elseif hemi == 'r'

        for ii = 1:8  
            vert_label( vert_label == temp.table( ii, 5)) = ii - 1;
        end

        for ii = 9:17  
            vert_label( vert_label == temp.table( ii, 5)) = ii + 1;
        end

    end

    %cmap = temp.table( 2:end, 1:3) ./ 255;
    %area_label = temp.struct_names(2:end);

    cmap = [ 239  59  44; 203  24  29; 153   0	13;
               0 109  44; 199 233 192;  35 139	69; 116	196	118
             254 178  76; 254 217 118;  84  39 143;  74  20 134; 141 211 199
             198 219 239; 158 202 225; 107 174 214;  66	146	198;  33 113 181; 8	69 148]./255; 

    area_label = {'mFus-faces','pFus-faces','IOG-faces',...
                  'OTS-bodies','ITG-bodies','MTG-bodies','LOS-bodies',...
                  'pOTS-characters','IOS-characters','CoS-places','TOS-places','hMT',...
                  'v1d','v2d','v3d','v1v','v2v','v3v'};
    

%-------------------------------

    case 'bensone'
        
    % surface labels - Benson eccentricity
    % NC Benson & J Winawer. Bayesian analysis of retinotopic maps. eLife (2018) 

    if ~isfile(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.benson14_varea.mgz']))
        error('Benson labels not found!');
    end

    surface_labels = MRIread(fullfile(freesurferPath, 'surf',...
        [hemi 'h.benson14_eccen.mgz']));
    vert_label = surface_labels.vol(:);

    % load the color map -logscale
    cmap = [];
    temp = logspace(0,3,90);
    temp1 = temp(end:-1:1);
    cm = jet(10^3);
    for ii = 1:90 
        cmap(ii, :) = cm(round(temp1(ii)), :);
    end

    area_label = 1 : max(vert_label);
    
%-------------------------------

    case 'hcp'
        
    % Surface labels - HCP MMP1 Atlas
    % M Glasser, et al. A multi-modal parcellation of human cerebral cortex. Nature (2016)

    if ~isfile(fullfile(freesurferPath, 'label',...
                    [hemi 'h.HCP-MMP1.annot']))
        error('HCP labels not found!');
    end

    [~, vert_label, temp] = read_annotation(fullfile(freesurferPath, 'label',...
                              [hemi 'h.HCP-MMP1.annot']));

    % Labels are numbered in some way (nothing from 1-180).. fixing this
    for ii = 1 : size( temp.table, 1)
        vert_label( vert_label == temp.table( ii, 5)) = ii - 1;
    end

    cmap = temp.table( 2:end, 1:3) ./ 255;
    area_label = temp.struct_names(2:end);


%-------------------------------

    case 'mix'

    % Benson for early visual, Rosenke for ventral, and Wang for parietal -
    % adapted from MGYR

    % Load individual atlases
    
    % 1. Wang

    if ~isfile(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.wang15_mplbl.mgz']))
        error('Wang labels not found!');
    end

    surface_labels_W = MRIread(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.wang15_mplbl.mgz']));
    vert_label_W = surface_labels_W.vol(:);

    vert_label_W = vert_label_W - 3;    % accounting for V1-3-v&d
    vert_label_W( vert_label_W < 9 ) = 0;  % removing ROIs till PHC2
   
    
    % 2. Rosenke - note: this atlas i different from visf
    % M Rosenke, et al. A cross-validated cytoarchitectonic atlas of the human ventral visual stream. NeuroImage (2018) 
   
    if ~isfile(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.rosenke18_vcatlas.mgz']))
        error('Rosenke labels not found!');
    end

    surface_labels_R = MRIread(fullfile(freesurferPath,'surf',...
        [hemi 'h.rosenke18_vcatlas.mgz']));
    vert_label_R = surface_labels_R.vol(:);

    vert_label_R(vert_label_R < 4) = 0; % removing early visual
    

    % 3. Benson

    if ~isfile(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.benson14_varea.mgz']))
        error('Benson labels not found!');
    end

    surface_labels_B = MRIread(fullfile(freesurferPath,'surf',...
        [hemi 'h.benson14_varea.mgz']));
    vert_label_B = surface_labels_B.vol(:);

    vert_label_B( vert_label_B > 3 ) = 0; % removing everything except early visual

    % Combining now..
    vert_label = vert_label_W + vert_label_R + vert_label_B;

    % To avoid two labels for the same vertices
    % Benson V1-3 have priority
    vert_label( vert_label_B ~= 0) = vert_label_B( vert_label_B ~= 0);

    % fix for minor discontinuity in V3d
    vert_label_W = surface_labels_W.vol(:);
    vert_label( vert_label_W == 6 & vert_label == 0) = 3;

    % Colors based on ScientificColourMaps
    % F Crameri, GE Shephard, and PJ Heron. The misuse of colour in science communication. Nature Communications (2020)


    cmap = [0.8585    0.8971    0.9154; % V1    Benson  1
            0.5661    0.7364    0.8202; % V2    Benson  2
            0.2689    0.5466    0.6909; % V3    Benson  3
            0.0422    0.3632    0.5678; % hOc4v Rosenke 4
            0.9338    0.8585    0.8154; % FG1   Rosenke 5
            0.8633    0.6740    0.5647; % FG2   Rosenke 6
            0.7896    0.5017    0.3354; % FG3   Rosenke 7
            0.7033    0.3264    0.1215; % FG4   Rosenke 8
            0.4912    0.1132    0.3900; % TO2   Wang    9 (orig:12)
            0.6486    0.2615    0.5485; % TO1   Wang    10
            0.7678    0.4061    0.6748; % LO2   Wang    11
            0.8587    0.5846    0.7874; % LO1   Wang    12
            0.9267    0.7758    0.8868; % V3B   Wang    13
            0.9641    0.9105    0.9438; % V3A   Wang    14
            0.9069    0.9404    0.8523; % IPS0  Wang    15
            0.7574    0.8539    0.6341; % IPS1  Wang    16
            0.5518    0.7102    0.3970; % IPS2  Wang    17
            0.3683    0.5632    0.2393; % IPS3  Wang    18
            0.2232    0.4385    0.1280; % IPS4  Wang    19
            0.0491    0.2967         0; % IPS5  Wang    20
            1.0000    1.0000    0.4002; % SPL1  Wang    21
            0.8981    0.8731    0.4087];% FEF   Wang    22 (25)    

    % old colormaps
    % cmap = [0.9882  0.5647  0.5647; % V1    Benson  1
    %         0.9882  0.7098  0.5255  % V2    Benson  2
    %             1     1     0.6000; % V3    Benson  3
    %         0.0784  0.9490  0.9176; % hOc4v Rosenke 4
    %         0.1137  0.9647  0.0549; % FG1   Rosenke 5
    %         0.9333  0.0549  0.9647; % FG2   Rosenke 6
    %         0.0863  0.1059  0.6824; % FG3   Rosenke 7
    %         0.4784  0.0471  0.0471; % FG4   Rosenke 8
    %         0.9000  0.6000  1.0000; % TO2   Wang    9 (orig:12)
    %         0.3438  0.3910  0.3668; % TO1   Wang    10
    %         0.6974  0.5129  0.8403; % LO2   Wang    11
    %         0.1978  0.8408  0.7445; % LO1   Wang    12
    %         0.3955  0.8227  0.5828; % V3B   Wang    13
    %         0.9686  0.3098  0.7353; % V3A   Wang    14
    %         0.2417  0.9409  0.7096; % IPS0  Wang    15
    %         0.2196  0.0039  0.7588; % IPS1  Wang    16
    %         0.2000  0.8000  0.5000; % IPS2  Wang    17
    %         0.9539  0.7295  0.7562; % IPS3  Wang    18
    %         0.2877  0.6154  0.4391; % IPS4  Wang    19
    %         0.2206  0.2239  0.4094; % IPS5  Wang    20
    %         0.9000  1.0000  0.7000; % SPL1  Wang    21
    %         0.2000  0.5000  0.2000];% FEF   Wang    22 (25)
    % 
    % % pastelize colors
    % cmap(1:3,:) = (cmap(1:3,:)).^0.8;
    % cmap(4:8,:) = (cmap(4:8,:)).^0.15;
    % cmap(9:end,:) = (cmap(9:end,:)).^.5;

    area_label = { 'V1','V2','V3', ...
                   'hOc4c','FG1','FG2','FG3','FG4',...
                   'TO2','TO1','LO2','LO1','V3B','V3A'...
                   'IPS0','IPS1','IPS2','IPS3','IPS4','IPS5',...
                   'SPL1','FEF'};


%-------------------------------

    case 'nsdfloc'

    % These ROIs are based on the NSD floc task
    % The threshold is currently set at t > 2 in contrast to t > 0 used in the paper 
    % The ROIs are further decided by taking mode across subjects
    % Refer script flocROI_fsAvg.m
    % E Allen, et al. A massive 7T fMRI dataset to bridge cognitive neuroscience and artificial intelligence. Nature Neuroscience (2022)

    % Load individual atlases based on hemi

    lab = load( fullfile( freesurferPath, 'nsd_floc/vlab_modeME.mat'));

    if isempty(lab)
        error('NSD floc labels not found!');;
    end

    vert_label = lab.vlab_modeME.([hemi 'h']);

    % Load Wang atlas for early visual areas

    if ~isfile(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.wang15_mplbl.mgz']))
        error('Wang labels not found!');
    end

    surface_labels_W = MRIread(fullfile(freesurferPath, 'surf',...
                    [hemi 'h.wang15_mplbl.mgz']));
    vert_label_W = surface_labels_W.vol(:);

    vert_label_W( vert_label_W > 8 ) = 0;  % removing ROIs beyond hV4

    % Combining the two atlases

    vert_label( vert_label_W == 1 | vert_label_W == 2)  = 13; % combining v1v-d
    vert_label( vert_label_W == 3 | vert_label_W == 4)  = 14; % combining v2v-d
    vert_label( vert_label_W == 5 | vert_label_W == 6)  = 15; % combining v3v-d
    vert_label( vert_label_W == 7)  = 16; % hV4


    % Colors based on ScientificColourMaps
    % F Crameri, GE Shephard, and PJ Heron. The misuse of colour in science communication. Nature Communications (2020)

    cmap= [  1.0000    1.0000    0.4002; % 1,  bodies (EBA)
             0.8981    0.8731    0.4087; % 2,  bodies (FBA-2)
             0.8955    0.6829    0.8399; % 3,  faces  (OFA)
             0.7874    0.4377    0.6974; % 4,  faces  (FFA-1)
             0.6216    0.2348    0.5210; % 5,  faces  (FFA-2)
             0.7574    0.8539    0.6341; % 6,  places (OPA)
             0.4879    0.6607    0.3376; % 7,  places (PPA)
             0.2719    0.4811    0.1654; % 8,  places (RSC)
             0.9338    0.8585    0.8154; % 9,  word  (OVWFA)
             0.8633    0.6740    0.5647; % 10, word  (VWFA-1)
             0.7896    0.5017    0.3354; % 11, word  (VWFA-2)
             0.7033    0.3264    0.1215; % 12, word  (mfs-words)
             0.8585    0.8971    0.9154; % 13, EV    (v1)
             0.5661    0.7364    0.8202; % 14, EV    (v2)
             0.2689    0.5466    0.6909; % 15, EV    (v3)
             0.0422    0.3632    0.5678];% 16, EV    (hv4)


    area_label = { 'EBA','FBA-2', ...
                   'OFA','FFA-1','FFA-2', ...
                   'OPA','PPA','RSC', ...
                   'OVWFA','VWFA-1','VWFA-2','mfs-words', ...
                   'v1','v2','v3','hv4' };



%-------------------------------

    otherwise
        error('Unexpected input. Expected: Wang/BensonV/BensonE/visf/HCP/mix/NSDfloc')
        
end
