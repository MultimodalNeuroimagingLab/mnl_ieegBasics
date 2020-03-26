% This script renders a brian surface with Destrieux maps
% Can potentially add electrodes on top
% dhermes & jvanderaar 2019, UMC Utrecht

% Make sure that this toolbox is in the path:
% can be cloned from: https://github.com/dorahermes/ecogBasicCode.git
addpath('/Fridge/users/jaap/github/ecogBasicCode/render/')

% adding VistaSoft path
addpath('/home/jaap/vistasoft/')


%% Render plain with used electrodes
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';
subjects = {'RESP0621'};
hemi_cap = {'R'};

% pick a viewing angle:
v_dirs = [90 0];%;90 0;90 -60;270 -60;0 0];

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);
    
%     % electrode locations name:
%     dataLocName = fullfile(dataRootPath,'soc_bids',['/sub-' subj],'/ses-01/','ieeg',...
%         ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);

%     % load electrode locations
%     loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
%     elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        figure
        ecog_RenderGifti(g) % render
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
%         % make sure electrodes pop out
%         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
%         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
%         ecog_Label(els,30,12) % add electrode positions
%         el_add(els(els_NBF{s},:),'k',30) % add electrode positions
        
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(dataRootPath,'derivatives','render',...
%             ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end

%% Render Destrieux with electrodes

% comments/questions

% - For sub-chaam, right hemisphere elektrodes, but the tsv-file states
% left hemisphere (other way around?) - is an issue in Gio's github umcu_elec
% - electrode positions TSV file is named ses-1 instead of ses-01

clear all
dataRootPath = '/Fridge/users/jaap/ccep/dataBIDS/';
% add vistasoft for read_annotation
addpath('/home/jaap/vistasoft/external/freesurfer');

subjects = {'RESP0621','chaam'};
sessions = {'1','01'};
hemi_cap = {'R','R'}; 
hemi_small = {'r','r'};

v_dirs = [90 0];%;90 0;90 -60;270 -60;0 0];


for s = 1%:2%1:length(subjects)
    % subject code 
    subj = subjects{s};
    ses_label = sessions{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi_cap{s} '.surf.gii']);
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'label',...
    [hemi_small{s} 'h.aparc.a2009s.annot']);
    % surface_labels = MRIread(surface_labels_name);
    [vertices, label, colortable] = read_annotation(surface_labels_name);
    vert_label = label; % these labels are strange and do not go from 1:76, but need to be mapped to the colortable
    % mapping labels to colortable
    for kk = 1:size(colortable.table,1) % 76 are labels
        vert_label(label==colortable.table(kk,5)) = kk;
    end
    
    % make a colormap for the labels
    cmap = colortable.table(:,1:3)./256;
    
    % electrode locations name:
    dataLocName = [dataRootPath 'sub-' subj '/ses-' ses_label '/ieeg/sub-' subj '_ses-' ses_label '_electrodes.tsv'];
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
  
        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,colortable.struct_names)
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);      
        ecog_Label(els,30,12) % add electrode positions

        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/Wang_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end




