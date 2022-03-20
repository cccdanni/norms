%% rsa
% Project: Social Norms Learning
% Author: Danni Chen 
% Update Date: Apr-5-2021

clear;clc;

project_folder='/home/chendanni/Documents/Norms/rawData/';
data_folder = fullfile(project_folder, 'EEGAnalysisResults');
toolbox_folder = '/home/chendanni/MATLAB/toolbox';
cd(project_folder);

%% add toolboxes
addpath (fullfile(project_folder,'additional_functions'));
addpath (fullfile(toolbox_folder,'fieldtrip-20210330'));

%% RSA for Prelearning
path_in=fullfile(project_folder,'EEGAnalysisResults');

% define the participants - 
subs = [1204:1242, 1244:1253];
all_subs = {};
for i = 1:length(subs)
    all_subs (i) =  { num2str(subs(i)) } ;
end

cd(path_in);

% define sliding windo
win = 0.02; % in sec 
slide = 0.01;

% define new Sampling rate
sr=100; % in Hz
step=win/(1/sr)+1;

% path_out = fullfile(project_folder,'RSA','data',strcat('all_trials_pre_sr',num2str(sr),'window',num2str(win.*1000),'slide',num2str(slide.*1000)));
path_out = fullfile(project_folder,'RSA','data',strcat('all_trials_pre_sr',num2str(sr),'window','full','slide','none'));

mkdir(path_out)

channels='all';

%define start and end of item window
t_start= 0;
t_end=0.996;
tois=t_start:1/sr:t_end;
t1=t_start:slide:(t_end-win);
t2=t1+win;
ind_t1=1:slide/(1/sr):((numel(tois)-win/(1/sr)));
ind_t2=ind_t1+win/(1/sr);
n_bins=numel(t1);

%% first step: correlate for each subject in all trials

for n=1:numel(all_subs)
    
    sel_sub = all_subs{n};
    sel_sub_folder = fullfile( data_folder, sel_sub );
    
    % load in the preprocessed data of each participant
    EEG = pop_loadset('filename',[sel_sub '_preLearningImplicit_EpochArtRej_FalseRespRej.set'],'filepath',sel_sub_folder);
    data = eeglab2fieldtrip ( EEG, 'preprocessing', 'none' );
    events = EEG.event;
    
    % resample to sr Hz
    cfg=[];
    cfg.resamplefs=sr;
    cfg.detrend='no';
    data=ft_resampledata(cfg,data);
    
%     % select channels
%     cfg=[];
%     cfg.channel=channels;
%     data=ft_preprocessing(cfg,data);
%     
%     % find trial in enco & cue/reco
%     enco_trials=find(data.trialinfo(:,5)==13 |data.trialinfo(:,5)==11);
%     cue_trials= find(data.trialinfo(:,5)==block(1) |data.trialinfo(:,5)==block(2));

    % prepare enco data
    cfg=[];
    cfg.latency=[t_start t_end+1/sr];
    data=ft_selectdata( cfg, data );
    
    % keep trials (face trials only)
    % deal with repeat epoch
    epochid = vertcat(events.epoch)';
    [uniqueid i j] = unique(epochid, 'first');
    indexdups = find( not (ismember(1:numel(epochid), i)));
    idex = 1:numel(epochid);
    idex = idex(~ismember (idex, indexdups) );
    events = events(idex);
    
    object_trials = find(arrayfun(@(events) ismember(9, events.bini), events));
    face_trials = find(arrayfun(@(events) ~ismember(9, events.bini), events));
    
    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=face_trials;
    data=ft_timelockanalysis(cfg, data);
    n_trials=numel(face_trials);
    
    % z trans across trials
    mean_trials = repmat(nanmean(data.trial,1),n_trials,1,1);
    std_trials  = repmat(nanstd(data.trial,1),n_trials,1,1);
    data.trial  = (data.trial-mean_trials)./std_trials;
    clear n_trials mean_trials std_trials
    
%     %%%%%%%%%%%%%
%     %prepare cue data
%     
%     cfg=[];
%     cfg.latency=[t_start_c t_end_c+1/sr];
%     datacue=ft_selectdata(cfg, data);
%     
%     cfg=[];
%     cfg.keeptrials='yes';
%     cfg.trials=cue_trials;
%     datacue=ft_timelockanalysis(cfg, datacue);
%     
%     n_trials=numel(cue_trials);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%
%     %     % z trans across trials
%     mean_trials=repmat(nanmean(datacue.trial,1),n_trials,1,1);
%     std_trials=repmat(nanstd(datacue.trial,1),n_trials,1,1);
%     datacue.trial=(datacue.trial-mean_trials)./std_trials;
%     clear n_trials mean_trials std_trials cue_trials
%     %
    
%     [items, use_enco, use_cue]=intersect(data.trialinfo(:,7),datacue.trialinfo(:,7));
    

    % average trials by items
    codelabel = unique( [data.trialinfo.codelabel].', 'rows', 'stable'); 
    codeseq = [];
    
    for code_nr = 1: size(codelabel, 2)
        this_code = codelabel{code_nr};
        this_code = str2double ( this_code(2:end) );
        codeseq(code_nr,1) = this_code;
        codeseq(code_nr,2) = code_nr;
    end
    
    [~, idx] = sort(codeseq(:,1)); 
    codeseq  = codeseq (idx, :); % label sequence arranging 
    
    trialinfo = table2struct( data.trialinfo );

    erpdata = struct('time',{},'label', {}, 'elec', {}, 'avg', {}, 'var', {}, 'dof', {}, 'dimord', {}, 'cfg', {}, 'codelabel', {});
    for code_nr = 1:size(codelabel, 2)
        clear this_data this_code this_trials
        
        this_code = ['s', num2str(codeseq(code_nr, 1))]; % this trial label
        this_trials = find(arrayfun(@(trialinfo) strcmp(this_code, trialinfo.codelabel), trialinfo));
        
        cfg=[];
        cfg.trials=this_trials;
        this_data = ft_timelockanalysis(cfg, data);
        this_data.codelabel = this_code;
       
        erpdata(code_nr)= this_data;
    end
    
    erpdata2 = zeros(size(erpdata(1).label, 1),  size(erpdata(1).time, 2), size(erpdata,2));
    for code_nr = 1:size(erpdata, 2)
        this_erp = erpdata(code_nr).avg;
        erpdata2(:, :, code_nr) = squeeze(this_erp);
    end
    

%     % resort trials/trialsinfo for matching trial order
%     data.trial=data.trial(use_enco,:,:);
%     data.trialinfo=data.trialinfo(use_enco,:);
%     datacue.trial=datacue.trial(use_cue,:,:);
%     datacue.trialinfo=datacue.trialinfo(use_cue,:);
    
%     trialinfo=data.trialinfo; %same as for datacue
    
%     data_enc_vec=zeros(size(trialinfo,1),n_bins_e,numel(data.label)*step);
%     for bin=1:n_bins_e
%         % vectorize sel_bins: data_vec(trials, nbins, features)
%         data_vec_tmp=data.trial(:,:,ind_t1_e(bin):ind_t2_e(bin));
%         data_enc_vec(:,bin,:)=reshape(data_vec_tmp,size(trialinfo,1),[]);
%     end
%     clear data_vec_tmp dataenc
%     
%     data_cue_vec=zeros(size(trialinfo,1),n_bins_c,numel(datacue.label)*step);
%     for bin=1:n_bins_c
%         % vectorize sel_bins: data_vec(trials, nbins, features)
%         data_vec_tmp=datacue.trial(:,:,ind_t1_c(bin):ind_t2_c(bin));
%         data_cue_vec(:,bin,:)=reshape(data_vec_tmp,size(trialinfo,1),[]);
%     end
%     clear data_vec_tmp datacue
%     
%     data_enc_vec2=reshape(data_enc_vec, size(data_enc_vec,1)*size(data_enc_vec,2),size(data_enc_vec,3));
%     data_cue_vec2=reshape(data_cue_vec, size(data_cue_vec,1)*size(data_cue_vec,2),size(data_cue_vec,3));

    erpdata3 = reshape (erpdata2, size(erpdata2, 1) * size(erpdata2, 2), size (erpdata2, 3));
    corr_res = corr(erpdata3, 'type', 'Spearman');
    
%     corr_tmp=corr(data_enc_vec2', data_cue_vec2', 'type','Spearman');
    
%     % size corr_cue_enc= trials enco x bins enco x trials cue x bins cue
%     corr_cue_enc=  reshape(corr_tmp,size(data_enc_vec,1),size(data_enc_vec,2),size(data_cue_vec,1),size(data_cue_vec,2));
    
%     clear corr_tmp data data_enc_vec data_cue_vec data_enc_vec2 data_cue_vec2
    
%     % fisher z transform correlations
%     corr_cue_enc=  0.5.*log(((ones(size(corr_cue_enc))+corr_cue_enc)./(ones(size(corr_cue_enc))-corr_cue_enc)));
    corr_res=  0.5.*log(((ones(size(corr_res))+corr_res)./(ones(size(corr_res))-corr_res)));

%     
%     % prepare data: get on diagonal i.e only in trial correlations
%     for tr=1:size(trialinfo,1)
%         corr_cue_enc_trial(tr,:,:)= squeeze(corr_res(tr,:,tr,:));
%     end
%     
%     corr_trials.corr_cue_enc_trial=corr_cue_enc_trial;
%     corr_trials.trialinfo=trialinfo;
%     corr_trials.time_item=[t1_e;t2_e];
%     corr_trials.time_cue=[t1_c;t2_c];
    
    save(fullfile(path_out, strcat(all_subs{n},'pre_alltrials_full.mat')),'corr_res','erpdata','erpdata3');
    
    clear corr_trials corr_cue_enc_trial trialinfo corr_cue_enc items corr_res
    clear cfg code_nr codelabel codeseq data EEG erpdata erpdata2 erpdata3 events face_trials object_trials idx i this_code this_data this erp this_trials
end