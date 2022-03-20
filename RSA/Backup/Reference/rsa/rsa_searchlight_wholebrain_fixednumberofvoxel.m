%addpath 'D:\matlab_tools\fieldtrip-20160122'
project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));
%%
% get beamformer filter to flexibly beam data

pathin = fullfile(project_folder,'data');
pathout= fullfile(project_folder,'data,','source','lcmv10_rsa','filter');
mkdir(pathout);

all_subs={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

load (fullfile(project_folder,'source','sourcemodel_10mm.mat'))
load(fullfile(project_folder,'scripts','additional_functions','MR_cap_standardElecpos.mat'))
load(fullfile(toolbox_folder, 'fieldtrip-20160122','template','headmodel','standard_bem.mat'))



for n=1:numel(all_subs) 
    sel_sub=all_subs{n};
    load(fullfile(path_in,strcat(sel_sub,'_data')));  
    cd (pathout);   

%covariance matrix

cfg=[];
cfg.covariance         = 'yes';
cfg.covariancewindow   = 'all';
cfg.keeptrials   = 'yes';
cfg.covariancewindow   = [0 4];
erp=ft_timelockanalysis(cfg,data);

% LCMV 
% 
cfg              = []; 
cfg.method       = 'lcmv';
cfg.elec         = elec_mnivol;
cfg.grid = sourcemodel;
cfg.vol          = vol;
cfg.lcmv.lambda       = '5%';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepfilter   = 'yes';
cfg.lcmv.realfilter   = 'yes';
cfg.lcmv.fixedori   = 'yes';
%cfg.normalize   = 'yes'; 
cfg.sensetype='eeg';
sourceAll = ft_sourceanalysis(cfg, erp);

% get only brain voxel (using aal atlas)
atlas=ft_read_atlas(fullfile(toolbox_folder,'fieldtrip-20160122','template','atlas','aal','ROI_MNI_V4.nii'));

% get only brain voxels in aal atlas
cfg=[];
   cfg.atlas        = atlas;
   cfg.roi   =atlas.tissuelabel(1:90);
   cfg.inputcoord   = 'mni'; 
 mask = ft_volumelookup(cfg,atlas);
%
mask_vec=mask(:);
mask_volume=atlas;
mask_volume=rmfield(mask_volume,'tissue');
mask_volume=rmfield(mask_volume,'tissuelabel');
mask_volume.mask=mask_vec;

% interpolate atlas mask brain to chosen source resolution
cfg=[];
cfg.parameter = 'mask';
cfg.interpmethod  = 'nearest';
brain_mask=ft_sourceinterpolate(cfg,mask_volume,sourceAll);
inside_vec=brain_mask.mask(:);

sourceAll.aal_brain=inside_vec;
save (fullfile(pathout,strcat(all_subs{n},'_sourceAll')),'sourceAll')
clear erp sourceAll 
end
%% rsa on source rois

% with the filter we can flexibly beam part of the brain (and only this
% part!)

path_data=fullfile(project_folder,'data');
path_filter=fullfile(project_folder,'data,','source','lcmv10_rsa','filter');

all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

sel_filter='_sourceAll'
%%% definition roi
sphere_radius=64; % in this case not radius, but number of next neighbours
num_roi_vox=sphere_radius;

%%% definition rsa
% define sliding window
cfg_rsa.win=0.2; % in sec
cfg_rsa.slide=0.01;
% define new Sampling rate
cfg_rsa.sr=100; % in Hz
cfg_rsa.step=cfg_rsa.win/(1/cfg_rsa.sr)+1;
cfg_rsa.trial_combi= 'in_trial';
cfg_rsa.corr_type='Spearman';


sel_effect='late_cluster';%or 'early_cluster'
%enco: define start and end of item window
t_start_e=0.2;
t_end_e=0.6;

%cue define start and end of item window
t_start_c=3.6;
t_end_c=4;


path_out =fullfile(project_folder,'data,','source','lcmv10_rsa', strcat(sel_effect,'sphere',num2str(sphere_radius),'_all_trials_item_cue_sr',num2str(cfg_rsa.sr),'window',num2str(cfg_rsa.win.*1000),'slide',num2str(cfg_rsa.slide.*1000),sel_filter,'_zall'));
mkdir(path_out)

% define roi for each source voxel (all sources same voxel definition so this can be done based on one subject)
 load(strcat(path_filter,vp{1},sel_filter))
    inside_vox=sourceAll.inside& sourceAll.aal_brain;
    % define for each voxel the roi
        % sphere definition
        dist2vox=squareform(pdist(sourceAll.pos(inside_vox,:)));
        num_vox=sum(inside_vox);
        sel_percentile=prctile(dist2vox,(num_roi_vox/num_vox).*100);
        
        sel_roi=(dist2vox<=repmat(sel_percentile,num_vox,1));
        
        

for n=1:numel(all_subs)
     sel_sub=all_subs{n};
    load(fullfile(path_in,strcat(sel_sub,'_data')));  
     display(strcat('starting:',sel_sub))

% beam data
 load(fullfile(path_filter,strcat(all_subs{n},sel_filter)))

  sub_vox=sourceAll.inside& sourceAll.aal_brain;

 %combine all trials in a matrix (chan*(time*trials))
trials= [data.trial{1,:}];
%
%combine all filters in one matrix(vox*chan)
 % select only filters of vox in roi sphere
 filters=vertcat(sourceAll.avg.filter{sub_vox,:});
 virtualsensors=filters*trials;
 beameddata=data;
 trialarray=reshape(virtualsensors,[sum(sub_vox),size(data.time{1,1},2),numel(data.trial)]);
 trial=squeeze(num2cell(trialarray,[1,2]))';

 beameddata.trial=trial;
 beameddata.label=cellstr(num2str(find(sub_vox)));
beameddata.time=data.time;
beameddata.trialinfo=data.trialinfo;
beameddata=rmfield(beameddata, 'sampleinfo');

data=beameddata;
clear beameddata vitrualsensors trialarray sig_vox sourceAll stat maxt filters enco_trials virtualsensors

% run rsa
%
    cfg=[];
    cfg.lpfilter='yes';
    cfg.lpfreq=40;
     data=ft_preprocessing(cfg,data);

 cfg=[];
    cfg.resamplefs=cfg_rsa.sr;
    cfg.detrend='no';
    data=ft_resampledata(cfg,data);


     % find trial in enco & cue/reco
    enco_trials=find(data.trialinfo(:,5)==13 |data.trialinfo(:,5)==11);
    cue_trials= find(data.trialinfo(:,5)==13 |data.trialinfo(:,5)==11);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % prepare enco data
    cfg=[];
    cfg.latency=[t_start_e t_end_e+1/cfg_rsa.sr];
    dataenc=ft_selectdata(cfg, data);

    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=enco_trials;
    dataenc=ft_timelockanalysis(cfg, dataenc);

    n_trials=numel(enco_trials);

    %%%%%%%%%%%%%%%%%%%%%%%%
   % z trans across trials
    mean_trials=repmat(nanmean(dataenc.trial,1),n_trials,1,1);
    std_trials=repmat(nanstd(dataenc.trial,1),n_trials,1,1);
    dataenc.trial=(dataenc.trial-mean_trials)./std_trials;
   clear n_trials mean_trials std_trials enco_trials

    %%%%%
   %%%%%%%%%%%%%
   %prepare cue data

       cfg=[];
    cfg.latency=[t_start_c t_end_c+1/cfg_rsa.sr];
    datacue=ft_selectdata(cfg, data);

    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=cue_trials;
    datacue=ft_timelockanalysis(cfg, datacue);
    n_trials=numel(cue_trials);

    %%%%%%%%%%%%%%%%%%%%%%%%
%     % z trans across trials
    mean_trials=repmat(nanmean(datacue.trial,1),n_trials,1,1);
    std_trials=repmat(nanstd(datacue.trial,1),n_trials,1,1);
    datacue.trial=(datacue.trial-mean_trials)./std_trials;
   clear n_trials mean_trials std_trials cue_trials

   [items, use_enco, use_cue]=intersect(dataenc.trialinfo(:,7),datacue.trialinfo(:,7));

   % resort trials/trialsinfo for matching trial order
   dataenc.trial=dataenc.trial(use_enco,:,:);
   dataenc.trialinfo=dataenc.trialinfo(use_enco,:);
   datacue.trial=datacue.trial(use_cue,:,:);
   datacue.trialinfo=datacue.trialinfo(use_cue,:);

    trialinfo=dataenc.trialinfo; %same as for datacue

    dataenc_all=dataenc;
    datacue_all=datacue;
tic
    parfor chan=1:numel(dataenc_all.label)
    cfg=[];
    cfg.channel=dataenc_all.label(sel_roi(:,chan));
    dataenc=ft_preprocessing(cfg,dataenc_all);
    datacue=ft_preprocessing(cfg,datacue_all);
    corr_trials{chan}=mf_trialwise_rsa(cfg_rsa,dataenc,datacue);
    end
  toc

  file_out=strcat(sel_sub,'_all_trials_item_cue_sr',num2str(cfg_rsa.sr),'window',num2str(cfg_rsa.win.*1000),'slide',num2str(cfg_rsa.slide.*1000),'_zall');
  save(fullfile(path_out, strcat(file_out,'_item_cue_alltrials')),'corr_trials');
    clear corr_trials corr_cue_enc_trial trialinfo corr_cue_enc items
end

% 
%% combine searchlight rsa results

path_filter=fullfile(project_folder,'data,','source','lcmv10_rsa','filter');

condition={'TBF_f','TBR_f','TBF_r','TBR_r'};
all_subs={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

sel_filter='_sourceAll'
%%% definition roi
sphere_radius=64; % radius in cm
sel_effect='late_cluster';
c_alpha=0.05;
% define sliding window
cfg_rsa.win=0.2; % in sec
cfg_rsa.slide=0.01;
% define new Sampling rate
cfg_rsa.sr=100; % in Hz
cfg_rsa.step=cfg_rsa.win/(1/cfg_rsa.sr)+1;
cfg_rsa.trial_combi= 'in_trial';
cfg_rsa.corr_type='Spearman';
path_out =fullfile(project_folder,'data,','source','lcmv10_rsa', strcat(sel_effect,'sphere',num2str(sphere_radius),'_all_trials_item_cue_sr',num2str(cfg_rsa.sr),'window',num2str(cfg_rsa.win.*1000),'slide',num2str(cfg_rsa.slide.*1000),sel_filter,'_zall'));

effects={'rehearsal','inhibition'};

for n=1:numel(all_subs)
    sel_sub=all_subs{n};
    cd(path_data)
    % load data of subjects
    file_def=strcat(sel_sub,'_all_trials_item_cue_sr',num2str(cfg_rsa.sr),'window',num2str(cfg_rsa.win.*1000),'slide',num2str(cfg_rsa.slide.*1000),'_zall_item_cue_alltrials')
    load(file_def)
    
    % load sourceAll
    load(fullfile(path_filter,strcat(sel_sub,sel_filter)))
    
    % average for conditions
    trialinfo=corr_trials{1}.trialinfo; %same as for datacue
    tbr_r_ind=trialinfo(:,5)==11&trialinfo(:,10)==1;
    tbr_f_ind=trialinfo(:,5)==11&trialinfo(:,10)==0;
    tbf_r_ind=trialinfo(:,5)==13&trialinfo(:,10)==1;
    tbf_f_ind=trialinfo(:,5)==13&trialinfo(:,10)==0;
    trial_def_vec=tbr_r_ind+(tbr_f_ind.*2)+(tbf_r_ind.*3)+(tbf_f_ind.*4);
    clear tbr_r_ind tbf_r_ind tbr_f_ind tbf_f_ind
    
    source_tmp=sourceAll;
    source_tmp=rmfield(source_tmp, 'avg');
    inside_def=sourceAll.inside & sourceAll.aal_brain;
    source_tmp.avg.pow=nan(size(sourceAll.inside));
    for c=1:numel(condition)
        sel_trials=trial_def_vec==c;
        avg_corr=zeros(numel(corr_trials),1);
        for chan=1:numel(corr_trials)
            avg_corr(chan)=squeeze(mean(mean(mean(corr_trials{chan}.corr_mat(sel_trials,:,:)))));
        end
        source_tmp.avg.pow(inside_def)=avg_corr;
        source_all{c,n}=source_tmp;
    end
    
end
mri=ft_read_mri(fullfile(toolbox_folder, 'fieldtrip-20160122','template','anatomy','single_subj_T1_1mm.nii'));

path_stat=fullfile(path_data,'stats');
mkdir(path_stat)
cd (path_stat)
for e=1:numel(effects)
    
    effect=effects{e};
    
    switch effect
        case 'inhibition'
            sourcega1=source_all(1,:);
            sourcega2=source_all(2,:);
        case 'rehearsal'
            sourcega1=source_all(3,:);
            sourcega2=source_all(4,:);
        otherwise
    end
    
    % run sourcestats
    design(1, :) = repmat(1:length(vp), 1, 2);
    design(2, :) = [ones(1, length(vp)) ones(1, length(vp)) * 2];
    
    cfg = [];
    cfg.dim         = sourcega1{1}.dim;
    %cfg.method      = 'analytic';
    cfg.statistic   = 'depsamplesT';
    cfg.parameter   = 'avg.pow';
    cfg.correctm    = 'cluster';
    cfg.method='montecarlo';
    cfg.numrandomization = 10000;
    cfg.alpha       = 0.05;
    cfg.clusteralpha = c_alpha;
    cfg.tail        = 0;
    cfg.design = design;
    cfg.ivar = 2;% indepva1
    cfg.uvar = 1;% units of observation  should not exist in indepsamp test
    
    stat = ft_sourcestatistics(cfg,sourcega1{:},sourcega2{:});
      
    save(strcat('stat_',effect,'_clusteralpha',num2str(c_alpha),'.mat'),'stat');
        
end

%% surface plots of results

path_stats=fullfile(project_folder,'source','lcmv10_rsa','early_clustersphere64_all_trials_item_cue_sr100window200slide10_sourceAll_zall,','stats');
dir_stat=dir(path_stats)

files_ind=find(strncmp('stat',{dir_stat(:).name},4));
stat_list={dir_stat(files_ind).name};
nsub=18; % for t value threshold
%
path_out=fullfile(path_stats,'figures');
mkdir(path_out)
mri=ft_read_mri(fullfile(toolbox_folder, 'fieldtrip-20160122','template','anatomy','single_subj_T1_1mm.nii'));

hemisphere={'left','right'};

views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90];

% mask with half maximum

load(fullfile(project_folder,'scripts','additional_functions','jet_grey_halfmax.mat'))

for s=1:numel(stat_list)
    load(strcat(path_stats,stat_list{s}))
    statmask=stat;
   cfg.parameter = 'stat';
    cfg.interpmethod  = 'nearest';
    
    statint=ft_sourceinterpolate(cfg, stat,mri);
    statmask.stat=(stat.stat<=tinv(0.1,nsub-1) |stat.stat>=tinv(0.9,nsub-1));
    statmask=ft_sourceinterpolate(cfg, statmask,mri);
    
    statint.mask=statmask.stat~=0;
    for h=1:numel(hemisphere)
        
        if sum(statint.mask)==0
            figure
        else
            
            cfg = [];
            cfg.method         = 'surface';
            cfg.funparameter   = 'stat';
            cfg.maskparameter  = 'mask';
            cfg.funcolorlim    = [-3 3];
            cfg.funcolormap    =jet_grey2;
            cfg.projmethod     = 'nearest';
            cfg.surffile       = strcat('D:\matlab_tools\carret',hemisphere{h},'.mat');
            cfg.surfdownsample = 5;
            %cfg.camlight = 'yes';
            cfg.camlight       = 'no';
            ft_sourceplot(cfg, statint);
 
            material DULL
            lighting gouraud;
            
            view(squeeze(views(h,1,:))');
            c1=camlight(0,0)
            set(c1, 'style', 'infinite');
            
            view(squeeze(views(h,2,:))');
            c2=camlight(0, 0)
            set(c2, 'style', 'infinite');
        end
        
        view(squeeze(views(h,3,:))');
        print('-f1','-dtiff','-r600',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_lat_3_01.tiff'))
        
        view(squeeze(views(h,4,:))');
        print('-f1','-r600','-dtiff',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_med_3_01.tiff'))
        
        view(squeeze(views(h,2,:))');
        print('-f1','-r600','-dtiff',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_vent2_3_01.tiff'))
        clear c1 c2
        
       close all
    end
    
end