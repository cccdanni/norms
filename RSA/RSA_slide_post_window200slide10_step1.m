% Project: Social Norms Learning (Time resolved)
% Function: [Postlearning] RSA analysis (Time-resolved) 
% Author: Danni Chen 
% Update Date: Apr-29-2021

clear;clc;

project_folder='/home/chendanni/Documents/Norms/rawData/';
data_folder = fullfile(project_folder, 'EEGAnalysisResults');
toolbox_folder = '/home/chendanni/MATLAB/toolbox';
cd(project_folder);

%% add toolboxes
addpath (fullfile(project_folder,'additional_functions')); 
addpath (fullfile(toolbox_folder,'fieldtrip-20210330'));

%% RSA for Postlearning
path_in=fullfile(project_folder,'EEGAnalysisResults');

% define the participants - 
subs = [1204:1242, 1244:1253];
badsubs = [1228, 1237, 1239, 1229];
subs = setdiff(subs, badsubs);
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

path_out = fullfile(project_folder,'RSA','data',strcat('all_trials_post_sr',num2str(sr),'window','200','slide','10'));

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
    EEG = pop_loadset('filename',[sel_sub '_postLearningImplicit_EpochArtRej_FalseRespRej.set'],'filepath',sel_sub_folder);
    data = eeglab2fieldtrip ( EEG, 'preprocessing', 'none' );
    events = EEG.event;
    
    % resample to sr Hz
    cfg=[];
    cfg.resamplefs=sr;
    cfg.detrend='no';
    data=ft_resampledata(cfg,data);
    
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
    % ERPData: averaged for each stim
    
    erpdata2 = zeros(size(erpdata(1).label, 1),  size(erpdata(1).time, 2), size(erpdata,2));
    erpdata_tmp = zeros(size(erpdata,2), size(erpdata(1).label, 1),  size(erpdata(1).time, 2));
    for code_nr = 1:size(erpdata, 2)
        this_erp = erpdata(code_nr).avg;
        erpdata2(:, :, code_nr) = squeeze(this_erp);
        erpdata_tmp(code_nr,:,:) = squeeze(this_erp);
    end
    % ERPData2: ERPData (Channels X Timepoints X stim).
    % ERPDataTmp: reshape ERPData (stim X channels X timepoints)

    
    erpdata_slidingwin=zeros(size(erpdata,2),n_bins,numel(data.label)*step); 
    % ERPData_slidingwin: stim_nr x bins_nr x (channels_nr x step)
    for bin=1:n_bins
        % vectorize sel_bins: data_vec(trials, nbins, features)
        data_tmp=erpdata_tmp(:,:,ind_t1(bin):ind_t2(bin));
        erpdata_slidingwin(:,bin,:)=reshape(data_tmp,size(erpdata,2),[]);
    end
    clear data_vec_tmp dataenc

%     erpdata3 = reshape (erpdata2, size(erpdata2, 1) * size(erpdata2, 2), size (erpdata2, 3));
%     % erpdata3: (timepoint X channels) X stimulus
%     corr_res = corr(erpdata3, 'type', 'Spearman'); % General correlation 
    
    erpdata4=reshape(erpdata_slidingwin, size(erpdata_slidingwin,1)*size(erpdata_slidingwin,2),size(erpdata_slidingwin,3));
    %erpdata4: (stim X bins) X  (channels X step)
    corr_res = corr(erpdata4', 'type', 'Spearman');
    % corr_res:  (stim X bins) X (stim X bins)
    corr_res=  reshape(corr_res,size(erpdata_slidingwin,1),size(erpdata_slidingwin,2),size(erpdata_slidingwin,1),size(erpdata_slidingwin,2));
    % corr_res: stim X bins X stim X bins
        
    % fisher z transform correlations
    corr_res=  0.5.*log(((ones(size(corr_res))+corr_res)./(ones(size(corr_res))-corr_res)));
    
%     % prepare data: get on diagonal i.e only in trial correlations
%     for tr=1:size(erpdata,2)
%         corr_res_trial(tr,:,:)= squeeze(corr_res(tr,:,tr,:));
%     end

    % prepare data: get on diagonal i.e only in trial correlations (not
    % across bins)
    for tr=1:size(corr_res,2)
        corr_res_trial(tr,:,:)= squeeze(corr_res(:,tr,:,tr));
    end

    corr_trials.corr_res_trial=corr_res_trial;
    corr_trials.trialinfo={erpdata(:).codelabel};
    corr_trials.time=[t1;t2];
    
    save(fullfile(path_out, strcat(all_subs{n},'post_alltrials_timeresolved.mat')),'corr_res','erpdata','corr_res_trial','erpdata4','corr_trials');

    clear corr_trials corr_cue_enc_trial trialinfo corr_cue_enc items corr_res
    clear cfg code_nr codelabel codeseq data EEG erpdata erpdata2 erpdata3 events face_trials object_trials idx i this_code this_data this erp this_trials
end


