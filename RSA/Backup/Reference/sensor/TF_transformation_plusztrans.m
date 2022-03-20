%% TIME-FREQUENCY TRANSFORMATION
%% first step:
project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%% this loads the data of individual subjects, calculates the actual wavelet transformation, baseline corrects and builds grand averages according to conditions
path_in=fullfile(project_folder,'data');
path_out=fullfile(project_folder,'data','TF');

mkdir(path_out);

all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

for n=1:numel(all_subs)
    sel_sub=all_subs{n};
    load(fullfile(path_in,strcat(sel_sub,'_data')));
    indata=data;
    clear data
    
    for i=1:length(indata.time)
        indata.trial{1,i} = horzcat(fliplr(indata.trial{1,i}(:,2:251)),indata.trial{1,i}(:,:),fliplr(indata.trial{1,i}(:,1250:1499)));
        indata.time{1,i}=-2:0.0040:5.9960;
    end
    
    cfg=[];
    cfg.method='wavelet';
    cfg.foi=1:0.5:30;
    cfg.toi=-1:0.05:5;
    cfg.width=5;
    cfg.output='pow';
    cfg.keeptrials='yes';
    data_all=ft_freqanalysis(cfg, indata);
    
    % ztrans across all trials
    cfg=[];
    cfg.time=[-0.5 4.5];
    data_z_all=z_trans_TF_seltime(cfg,data_all);
    
    % Define the conditions based on triggerinformation from trialinfo
    %(11 = TBR, 13 = TBF, 1 = remembered, 0 = forgotten)
    conds.enc_TBR_r=find(indata.trialinfo(:,5)==11 & indata.trialinfo(:,10)==1);
    conds.enc_TBR_f=find(indata.trialinfo(:,5)==11 & indata.trialinfo(:,10)==0);
    conds.enc_TBF_r=find(indata.trialinfo(:,5)==13 & indata.trialinfo(:,10)==1);
    conds.enc_TBF_f=find(indata.trialinfo(:,5)==13 & indata.trialinfo(:,10)==0);
    
    
    % creates a structure that stores the names of the conditions
    subjinfo.condnames=fieldnames(conds);
    subjinfo.subjects=all_subs;
    % initiate wavelet transformation per condition
    for i=1:length(subjinfo.condnames)
        % stores the number of trials in each condition
        subjinfo.trialnum(n,i)=length(getfield(conds,subjinfo.condnames{i}));
        trials=getfield(conds,subjinfo.condnames{i});
        
        %average data
        data_z=data_all;
        data_z.powspctrm=(squeeze(nanmean(data_z_all.powspctrm(trials,:,:,:))));
        data_z.dimord='chan_freq_time';
        
        mkdir(fullfile(path_out,subjinfo.condnames{i}));
        save(fullfile(path_out,subjinfo.condnames{i},strcat(all_subs{n},'_pow_z')),'data_z');
        clear data_z
        
    end
    clear data_all data_z_all
end

save(strcat(path_out,'subjinfo'),'subjinfo');

%% Build grandaverages over subjects
path_in=fullfile(project_folder,'data','TF');
path_out=fullfile(project_folder,'data','TF','GA');
mkdir(path_out)
%
load(strcat(path_in,'subjinfo.mat'));

%select the conditions you want to average
sample.condnames=subjinfo.condnames(selcond);

all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

% extract and load the individual files from the condition directories
for i=1:length(sample.condnames)
    
    for n=1:length(all_subs)
        load(strcat(path_in,sample.condnames{i},'\',all_subs{n},'_pow_z.mat'));
        pow_z{n}=data_z;
    end
    
    %grandaveraging
    cfg=[];
    cfg.keepindividual = 'yes'; %you want to keep this in order to be able to callculate stats
    data = ft_freqgrandaverage(cfg,pow_z{1,:});
    save(strcat(path_out,'GA_z_',sample.condnames{i}),'data');
end
