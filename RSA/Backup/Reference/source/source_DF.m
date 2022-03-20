%standard bem with standard positions

project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));
%% get standard positions by normalizing my electrode positions with my mr
path_out=fullfile(project_folder,'source')
mkdir(path_out)

%load standard bem
load(fullfile(toolbox_folder, 'fieldtrip-20160122','template','headmodel','standard_bem.mat'))

%load elec
load(fullfile(project_folder,'scripts','additional_functions','MR_cap_standardElecpos.mat'))

cfg = [];
cfg.grid.resolution=10;
cfg.grid.unit  ='mm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -10;
cfg.vol        = vol;
sourcemodel1 = ft_prepare_sourcemodel(cfg);
mri=ft_read_mri(fullfile(toolbox_folder, 'fieldtrip-20160122','template','anatomy','single_subj_T1_1mm.nii'));
mri.coordsys='mni';
cfg=[];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = sourcemodel1;
cfg.grid.nonlinear = 'no'; % use non-linear normalization
cfg.mri            = mri;
sourcemodel        = ft_prepare_sourcemodel(cfg);

%% check if all fits
figure
ft_plot_mesh(vol.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
hold on;
% electrodes
ft_plot_sens(elec_mnivol,'label','on');

figure;
ft_plot_mesh(vol.bnd(1,3),'edgecolor','none', 'facealpha', 0.8,'facecolor', 'red')
hold on
%ft_plot_sens(elec_new);
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:).*10);
plot3(elec_mnivol.elecpos(:,1),elec_mnivol.elecpos(:,2),elec_mnivol.elecpos(:,3),'g*');% insert elecR
%plot3(elec.pnt(:,1),elec.pnt(:,2),elec.pnt(:,3),'g*');% insert elecR

%save 10mm grid
save(fullfile(path_out,'sourcemodel_10mm.mat'),'sourcemodel');%


%% beam data
path_in=fullfile(project_folder,'data');
path_out=fullfile(project_folder,'data','source','lcmv10');
mkdir(path_out)

all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

load (fullfile(project_folder,'source','sourcemodel_10mm.mat'))
load(fullfile(project_folder,'scripts','additional_functions','MR_cap_standardElecpos.mat'))
load(fullfile(toolbox_folder, 'fieldtrip-20160122','template','headmodel','standard_bem.mat'))

for n=1:numel(all_subs)
    sel_sub=all_subs{n};
    load(fullfile(path_in,strcat(sel_sub,'_data')));
    
    cd (pathout);
    
    %covariance matrix
    cfg=[];
    cfg.lpfilter='yes';
    cfg.lpfreq=40;
    data=ft_preprocessing(cfg,data);
    
    cfg=[];
    cfg.covariance         = 'yes';
    cfg.covariancewindow   = 'all';
    cfg.keeptrials   = 'yes';
    cfg.covariancewindow   = [0 4];
    erp=ft_timelockanalysis(cfg,data);
    
    
    % LCMV
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
    cfg.sensetype='eeg';
    sourceAll = ft_sourceanalysis(cfg, erp);
    
    clear erp
    
    %combine all trials in a matrix (chan*(time*trials))
    trials= [data.trial{1,:}];
    
    %combine all filters in one matrix(insidepos*chan)
    insidepos=find(sourceAll.inside);
    filters=vertcat(sourceAll.avg.filter{insidepos,:});
    
    virtualsensors=filters*trials;
    beameddata=data;
    trialarray=reshape(virtualsensors,[numel(insidepos),size(data.time{1,1},2),numel(data.trial)]);
    trial=squeeze(num2cell(trialarray,[1,2]))';
    
    beameddata.trial=trial;
    beameddata.label=cellstr(num2str(insidepos));
    beameddata.time=data.time;
    beameddata.trialinfo=data.trialinfo;
    
    beameddata=rmfield(beameddata, 'sampleinfo');
    
    clear virtualsensors trials trialarray filters trial
    
    mkdir(fullfile(pathout,all_subs{n}));
    save (fullfile(pathout,all_subs{n},'sourceAll'),'sourceAll')
    
    clear sourceAll
    
    %
    cfg=[];
    cfg.method = 'wavelet';
    cfg.width = 5;
    cfg.output     = 'pow';
    cfg.foi = 2:1:30;
    cfg.toi =-0.5:0.05:4.5;
    cfg.keeptrials = 'yes';
    freq = ft_freqanalysis(cfg, beameddata);
    clear  beameddata
    
    cfg=[];
    cfg.time=[0 4];
    freq=z_trans_TF_seltime(cfg,freq);
    %split data in conditions
    
    % Define the conditions based on triggerinformation from trialinfo
    %(11 = TBR, 13 = TBF,  1 = remembered, 0 = forgotten)
    conds.enc_TBR_r=find(indata.trialinfo(:,5)==11 & indata.trialinfo(:,10)==1);
    conds.enc_TBR_f=find(indata.trialinfo(:,5)==11 & indata.trialinfo(:,10)==0);
    conds.enc_TBF_r=find(indata.trialinfo(:,5)==13 & indata.trialinfo(:,10)==1);
    conds.enc_TBF_f=find(indata.trialinfo(:,5)==13 & indata.trialinfo(:,10)==0);
    
    subjinfo.condnames=fieldnames(conds);
    for i=1:length(subjinfo.condnames)
        
        % stores the number of trials in each condition
        trials=getfield(conds,subjinfo.condnames{i});
        
        freq_cond=freq;
        freq_cond.dimord='chan_freq_time';
        freq_cond.powspctrm=squeeze(nanmean(freq.powspctrm(trials,:,:,:),1));
        
        save (fullfile(path_out,all_subs{n},strcat(subjinfo.condnames{i}, 'beamed_TF.mat')), 'freq_cond')
        clear freq_cond trials
    end
end

