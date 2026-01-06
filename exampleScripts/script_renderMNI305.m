
%% example script to render electrode positions in fsaverage

% dependencies:
%   vistasoft: https://github.com/vistalab/vistasoft
%   cvncode: https://github.com/cvnlab/cvncode
%   ieegbasics: https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics

% Dora Hermes, Mayo Clinic 2025

%% Set data path


localDataPath = 'where I put the iEEG-NSD data folder';

fsaverage_folder = 'my fsaverage folder';

%% Load electrode positions in MNi305 for a subject
% these MNI305 electrode positions can be genrated by using: mnl_ieegBasics/functions/ieeg_mni305ThroughFsSphere.m

subjects = {'01'};

ss = 1;
sub_label = subjects{ss};

% get elecmatrix (in MNI305)
electrodes_tsv_name = fullfile(localDataPath,['sub-' sub_label],'ses-ieeg01','ieeg',...
    ['sub-' sub_label '_ses-ieeg01_space-MNI305_electrodes.tsv']);
loc_info = ieeg_sortElectrodes(electrodes_tsv_name, channels_table, 0);


%% pial surface

% load inflated surfaces and make a gifti
[vertex_coords, faces] = read_surf(fullfile(fsaverage_folder,'surf','lh.pial'));
g.vertices = vertex_coords; 
g.faces = faces+1;
gL_infl = gifti(g); clear g faces vertex_coords
[vertex_coords, faces] = read_surf(fullfile(fsaverage_folder,'surf','rh.pial'));
g.vertices = vertex_coords; 
g.faces = faces+1;
gR_infl = gifti(g); clear g faces vertex_coords

% curvature labels
sL_labels = read_curv(fullfile(fsaverage_folder,'surf','lh.curv'));
sR_labels = read_curv(fullfile(fsaverage_folder,'surf','rh.curv'));

for hh = 1%:2

    if hh==1
        hemi = 'l';
        g = gL_infl;
        sulcal_labels = sL_labels;
    elseif hh==2
        hemi = 'r';
        g = gR_infl;
        sulcal_labels = sR_labels;
    end

    views_plot = {[270,0],[90,0],[180,-90],[-40,-30],[40,-30],[-60,30],[60,30]};

    % make a plot with electrode dots
    for vv = 1%:length(views_plot)
        v_d = [views_plot{vv}(1),views_plot{vv}(2)];

        % check which hemispheres have electrodes
        xyzmni_plot = NaN(height(loc_info),3);
        for kk = 1:height(loc_info)
            % if ~isnan and correct hemi
            if ~isnan(loc_info.vertex_fsaverage(kk)) && ismember(loc_info.hemisphere{kk},upper(hemi)) 
                xyzmni_plot(kk,:) = g.vertices(loc_info.vertex_fsaverage(kk),:);
            end
        end

        % get the inflated coordinates
        els = xyzmni_plot;

        % els popout into viewing direction
        % !!! when using popout, do not rotate the rendering !!!
        a_offset = .1*max(abs(els(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els_pop = els+repmat(a_offset,size(els,1),1);
            
        figure
        tH = ieeg_RenderGifti(g,sulcal_labels,[]); hold on
        ieeg_label(els_pop,10,10,loc_info.name)
        ieeg_viewLight(v_d(1),v_d(2)) % correct viewing angle
        tH.FaceAlpha = 0.5;

    end
end


%% inflated brain

% load inflated surfaces and make a gifti
[vertex_coords, faces] = read_surf(fullfile(fsaverage_folder,'surf','lh.inflated'));
g.vertices = vertex_coords; 
g.faces = faces+1;
gL_infl = gifti(g); clear g faces vertex_coords
[vertex_coords, faces] = read_surf(fullfile(fsaverage_folder,'surf','rh.inflated'));
g.vertices = vertex_coords; 
g.faces = faces+1;
gR_infl = gifti(g); clear g faces vertex_coords

% curvature labels
sL_labels = read_curv(fullfile(fsaverage_folder,'surf','lh.curv'));
sR_labels = read_curv(fullfile(fsaverage_folder,'surf','rh.curv'));

for hh = 1%:2

    if hh==1
        hemi = 'l';
        g = gL_infl;
        sulcal_labels = sL_labels;
    elseif hh==2
        hemi = 'r';
        g = gR_infl;
        sulcal_labels = sR_labels;
    end

    views_plot = {[270,0],[90,0],[180,-90],[-40,-30],[40,-30],[-60,30],[60,30]};

    % make a plot with electrode dots
    for vv = 1%:length(views_plot)
        v_d = [views_plot{vv}(1),views_plot{vv}(2)];

        % check which hemispheres have electrodes
        xyzmni_plot = NaN(height(loc_info),3);
        for kk = 1:height(loc_info)
            % if ~isnan and correct hemi
            if ~isnan(loc_info.vertex_fsaverage(kk)) && ismember(loc_info.hemisphere{kk},upper(hemi)) 
                xyzmni_plot(kk,:) = g.vertices(loc_info.vertex_fsaverage(kk),:);
            end
        end

        % get the inflated coordinates
        els = xyzmni_plot;

        % els popout into viewing direction
        % !!! when using popout, do not rotate the rendering !!!
        a_offset = .1*max(abs(els(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els_pop = els+repmat(a_offset,size(els,1),1);
            
        figure
        tH = ieeg_RenderGifti(g,sulcal_labels,[]); hold on
        ieeg_label(els_pop,10,10,loc_info.name)
        ieeg_viewLight(v_d(1),v_d(2)) % correct viewing angle

    end
end

%% show electrode positions on flatmap
% only works if you created a flatmap and saved it in your fs average folder

% curvature labels
sL_labels = read_curv(fullfile(fsaverage_folder,'surf','lh.curv'));
sR_labels = read_curv(fullfile(fsaverage_folder,'surf','rh.curv'));

for hh = 1%:2

    if hh==1
        hemi = 'l';
        sulcal_labels = sL_labels;
        surf = cvnreadsurface([],'lh','full.flat.patch.3d','','surfdir',fullfile(fsaverage_folder,'surf'));
    elseif hh==2
        hemi = 'r';
        sulcal_labels = sR_labels;
        surf = cvnreadsurface([],'rh','full.flat.patch.3d','','surfdir',fullfile(fsaverage_folder,'surf'));
    end

    % check which hemispheres have electrodes
    xyzmni_plot = NaN(height(loc_info),3);
    for kk = 1:height(loc_info)
        % if ~isnan and correct hemi
        if ~isnan(loc_info.vertex_fsaverage(kk)) && ismember(loc_info.hemisphere{kk},upper(hemi)) 
            xyzmni_plot(kk,:) = surf.vertices(loc_info.vertex_fsaverage(kk),:);
        end
    end

    % get the flatmap coordinates
    els = xyzmni_plot;

    % render flatmap
    figure
    tH = ieeg_RenderGifti(surf,sulcal_labels,[]); hold on
    ieeg_viewLight(0,90)

    % add electrode positions 
    for kk = 1:height(loc_info)
        if ~isnan(els(kk,1))
            xyz = els(kk,:);
            el_size = 10;
            el_color = [.1 .1 .1];
            s = scatter3(xyz(1),xyz(2),xyz(3),'o','Filled');
            s.SizeData = el_size.^2;
            s.MarkerFaceColor = el_color;
            s.MarkerEdgeColor = [0 0 0];
            s.MarkerFaceAlpha = 0.5;
        end
    end


end



