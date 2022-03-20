% This is a script for conducting RSA analysis for Norms Learning (2-D)
% 1) RSA between morphed faces & experimental faces
% 2) statistical analysis - a) higher, lower & consistent, b) ttest over
% zero
% Author: Danni Chen, Ziqing Yao & Xiaoqing Hu (The University of Hong Kong)
% Update date: 2022-03-01


%% Basic setting up

subs = [1205:1242,1244:1253];
badsubs   = [1201:1204, 1214, 1215, 1238]; % 1201 - 1204: HK Local; 1214, 1215 & 1238 - visually detected bad subjects
subs = setdiff (subs, badsubs);

work_dir = '/home/chendanni/Documents/Norms/analysis/';
cd (work_dir);

task = "post-minus-pre:";
erp  = "full_range";

stat_test = "yes"; 
perm_test = "yes";
ztrans = "yes"; 
average_erp = "no";
exclude = "exclude_object_control_faces";
n_all = 70;
n_cond = 7; % this is when we are excluding control faces
bin_list = [1,2,3,4,5,6,8];

%% Subject Filename(s)

% Load Subject Datafiles 
data_folder = fullfile (work_dir, 'EEGAnalysisResults', 'ERP'); % set directory of data set (default: pwd)
toolbox_folder = '/home/chendanni/MATLAB/toolbox';

% Save RSA results
saveLocation = fullfile (pwd, 'EEGAnalysisResults','RSA');
if ~exist (saveLocation)
    mkdir (saveLocation);
end

% define the participants 
all_subs = {};
for i = 1:length(subs)
    all_subs (i) =  { num2str(subs(i)) } ;
end

%% define RSA parameters
% define sliding window
win = 0.20;
slide = 0.01;

% define new Sampling rate
sr=100; % in Hz
step=win/(1/sr)+1;

%define start and end of item window
t_start= 0;
t_end=0.996;
tois=t_start:1/sr:t_end;
t1=t_start:slide:(t_end-win);
t2=t1+win;
ind_t1=1:slide/(1/sr):((numel(tois)-win/(1/sr)));
ind_t2=ind_t1+win/(1/sr);
n_bins=numel(t1);

channels='all';

% permutation test setup
rand=5000;
p_cluster=0.05;
p_first=0.01;


%% define path_in & path_out

path_in_pre   = fullfile(saveLocation, 'Results',strcat('RSA2D_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
path_in_post  = fullfile(saveLocation, 'Results',strcat('RSA2D_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
savName_pre = strcat ('/', 'Results_RSA2D_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 
savName_post = strcat ('/', 'Results_RSA2D_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 
path_out = fullfile(saveLocation, 'Stats',strcat('RSA2D_Diff_', task, '_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
mkdir(path_out)


%% add toolboxes
addpath (genpath(fullfile(work_dir,'MyScripts')));
addpath (fullfile(toolbox_folder,'fieldtrip-20210330'));

% %% First step: correlate for each subject in all trials
% 
% cd (path_in);
% 
% for n = 1:numel(all_subs)
%     
%     TStart = tic;
%     
%     sel_sub = all_subs{n};
%     sel_sub_folder = fullfile( data_folder, sel_sub );
%     
%     fprintf (strcat ("Running: Subject ", sel_sub, "\n")); 
%     
%     
%     %% EEG data loading & processing
%     % load in the preprocessed data of each participant
%     EEG = pop_loadset('filename',char(strcat( sel_sub, '_', task, '_EpochArtRej_FalseRespRej.set')),'filepath',fullfile(sel_sub_folder,'/'));
%     data = eeglab2fieldtrip ( EEG, 'preprocessing', 'none' );
%     events = EEG.event;
%    
%     % resample to sr Hz
%     cfg=[];
%     cfg.resamplefs=sr;
%     cfg.detrend='no';
%     data=ft_resampledata(cfg,data);
%     
%     % select channels if needed
%     if ~strcmp(channels, 'all')
%         cfg=[];
%         cfg.channel=channels;
%         data=ft_preprocessing(cfg,data);
%     end
% 
%     % prepare task data (after stimuli onset)
%     cfg=[];
%     cfg.latency=[t_start t_end+1/sr];
%     data=ft_selectdata( cfg, data );
%     
%     % remove unused trials if needed, here we removed object trials and
%     % duplicated epochs    
%     epochid = vertcat(events.epoch)';
%     [uniqueid i j] = unique(epochid, 'first');
%     indexdups = find( not (ismember(1:numel(epochid), i)));
%     idex = 1:numel(epochid);
%     idex = idex(~ismember (idex, indexdups) );
%     events = events(idex); % clean up repeated epochs
%     
%     tinfo = table2struct (data.trialinfo);
%     face_trials = ismember ([tinfo.bini], bin_list);
%     cfg=[];
%     cfg.keeptrials='yes';
%     cfg.trials=face_trials;
%     data=ft_timelockanalysis(cfg, data); % keep only experimental face trials, exclude both object trials and control faces trials
%     n_trials=size(data.trial,1); % number of trials
%     n_stim=length (unique ([data.trialinfo.codelabel])); % number of stimuli
%     
%     % z trans across trials
%     if (strcmp (ztrans, 'yes'))
%         mean_trials = repmat(nanmean(data.trial,1),n_trials,1,1);
%         std_trials  = repmat(nanstd(data.trial,1),n_trials,1,1);
%         data.trial  = (data.trial-mean_trials)./std_trials;
%         clear n_trials mean_trials std_trials 
%     end
%     % As we applied a z-transformation to the EEG activity at every time point 
%     % across all trials of all conditions in order to remove ERP componentsprior 
%     % to calculating correlations (Fellner et al., 2020, Current Biology)    
%     
%     
%     % find trial in experimental & morphed faces
%     tinfo = table2struct (data.trialinfo);
%     exp_trials = ismember ([tinfo.bini], 1:6); 
%     mor_trials = ismember ([tinfo.bini], 8);
%     
%     %%%%%%%%%%%%%%
%     % prepare experimental data
%     cfg = [];
%     cfg.keeptrials = 'yes';
%     cfg.trials = exp_trials;
%     dataexp = ft_timelockanalysis (cfg, data);
%     n_trials_exp = numel (exp_trials);
%     n_stim_exp   = numel (unique([dataexp.trialinfo.codelabel]));
%     
% 
%     %%%%%%%%%%%%%
%     % prepare morphed data
%     cfg = [];
%     cfg.keeptrials = 'yes';
%     cfg.trials = mor_trials;
%     datamor = ft_timelockanalysis (cfg, data);
%     n_trials_mor = numel (mor_trials);
%     n_stim_mor = numel (unique([datamor.trialinfo.codelabel]));
%     
% %     % z trans across trials
% %     if (strcmp (ztrans, 'yes'))
% %         mean_trials = repmat(nanmean(datamor.trial,1),n_trials_mor,1,1);
% %         std_trials  = repmat(nanstd(datamor.trial,1),n_trials_mor,1,1);
% %         datamor.trial  = (datamor.trial-mean_trials)./std_trials;
% %         clear n_trials mean_trials std_trials 
% %     end
%     
%     %%%%%%%%%%%%
%     
%     trialinfo_exp = dataexp.trialinfo; 
%     data_exp_vec=zeros(size(trialinfo_exp,1),n_bins,numel(dataexp.label)*step); % trial x bins x features (step x )
%     for bin=1:n_bins
%         % vectorize sel_bins: data_exp(trials, nbins, features)
%         data_exp_tmp=dataexp.trial(:,:,ind_t1(bin):ind_t2(bin));
%         data_exp_vec(:,bin,:)=reshape(data_exp_tmp,size(trialinfo_exp,1),[]);
%         clear data_exp_tmp
%     end
%     clear dataexp
%     
%     trialinfo_mor = datamor.trialinfo;
%     data_mor_vec=zeros(size(trialinfo_mor,1),n_bins,numel(datamor.label)*step);
%     for bin=1:n_bins
%         % vectorize sel_bins: data_mor(trials, nbins, features)
%         data_mor_tmp=datamor.trial(:,:,ind_t1(bin):ind_t2(bin));
%         data_mor_vec(:,bin,:)=reshape(data_mor_tmp,size(trialinfo_mor,1),[]);
%         clear data_mor_tmp
%     end
%     clear datamor
%     
%     data_exp_vec2=reshape(data_exp_vec, size(data_exp_vec,1)*size(data_exp_vec,2),size(data_exp_vec,3));
%     data_mor_vec2=reshape(data_mor_vec, size(data_mor_vec,1)*size(data_mor_vec,2),size(data_mor_vec,3));
%     corr_tmp=corr(data_exp_vec2', data_mor_vec2', 'type','Spearman');
% 
%     % size corr_exp_mor= trials exp x bins exp x trials morph x bins morph
%     corr_exp_mor =  reshape(corr_tmp,size(data_exp_vec,1),size(data_exp_vec,2),size(data_mor_vec,1),size(data_mor_vec,2));
%     
%     clear corr_tmp data data_exp_vec data_mor_vec data_exp_vec2 data_mor_vec2
% 
%     % fisher z transform correlations
%     corr_exp_mor =  0.5.*log(((ones(size(corr_exp_mor))+corr_exp_mor)./(ones(size(corr_exp_mor))-corr_exp_mor)));
% 
%     % prepare data: get on diagonal i.e only in trial correlations
%     % we have multiple morphed faces trials, they will be averaged
%     % size corr_exp_mor_trial trial_exp * bins_exp * bins_morphs
%     for tr=1:size(trialinfo_exp,1)
%         corr_tmp = mean(corr_exp_mor(tr,:,:,:), 3);
%         corr_exp_mor_trial(tr,:,:)= squeeze(corr_tmp);
%         clear corr_tmp
%     end
%         
%     % reshape the trial_based correlation matrix to stimuli based matrix
%     corr_exp_mor_stim = nan (n_stim_exp, n_bins, n_bins);
%     trialinfo = trialinfo_exp;
%     triallabel = [string(trialinfo.codelabel)];
%     
%     stim_list = [];
%     for tmp_i = 1:n_stim_exp
%         stim_list = [stim_list; sscanf(triallabel(tmp_i),"s%d")];
%     end
%     stim_list = sort (stim_list)';
%     
%     stim_info = struct;
%     for i_stim = 1:n_stim_exp 
%             i_trials = find (strcmp (strcat('s', num2str(stim_list(i_stim))) , triallabel))';
%             corr_tmp = corr_exp_mor_trial (i_trials, :, :);
%             corr_exp_mor_stim (i_stim, :, :) = squeeze (mean (corr_tmp, 1));    
%             stim_info(i_stim).stim_id = i_stim;
%             stim_info(i_stim).stim_label = strcat('s', num2str(stim_list(i_stim)));
%             stim_info(i_stim).trials_i =  table2array(trialinfo(i_trials, 'bepoch'))';
%             stim_info(i_stim).bin_list = table2array(trialinfo(i_trials, 'bini'));
%             stim_info(i_stim).bini = stim_info(i_stim).bin_list(1);
%             clear i_trials corr_tmp;
%     end
%     
%     corr_trials.corr_exp_mor_trial=corr_exp_mor_trial;
%     corr_trials.trialinfo=trialinfo_exp;
%     corr_trials.time_exp=[t1;t2];
%     corr_trials.time_mor=[t1;t2];
%     corr_trials.corr_exp_mor_stim=corr_exp_mor_stim;
%     corr_trials.stiminfo=stim_info;
%    
%     save(fullfile(path_out, strcat(savName, all_subs{n},'.mat')), 'corr_trials');
%     
%     clear corr_trials corr_exp_mor_trial corr_exp_mor_stim trialinfo_exp trialinfo_cue corr_exp_mor
% 
%     toc(TStart);
%     
% end


%% load trial correlations and draw random in each subject



% combine all correlation in one
for n=1:numel(all_subs)
    
    pre_data = load(fullfile(path_in_pre, strcat(savName_pre, all_subs{n},'.mat')));
    post_data = load(fullfile(path_in_post, strcat(savName_post, all_subs{n},'.mat')));
    
%     % option 1 using all trials:
%     trial_corr_all{n}=pre_data.corr_trials.corr_exp_mor_trial - pos;
%     trialinfo_all{n}=corr_trials.trialinfo;
%     
    % option 2 using averaged trials:
    trial_corr_all{n}=post_data.corr_trials.corr_exp_mor_stim - pre_data.corr_trials.corr_exp_mor_stim;
    trial_corr_pre{n}=pre_data.corr_trials.corr_exp_mor_stim;
    trial_corr_post{n}=post_data.corr_trials.corr_exp_mor_stim;
    trialinfo_all{n} =struct2table(pre_data.corr_trials.stiminfo); % should be the same between pre-learning and post-learning
    
    % clear corr_trials;
end
selvps = numel(all_subs);

% prepare data for statistical analysis

for n=1:numel(all_subs)
    
    trial_corr=trial_corr_all{n};
    trialinfo=trialinfo_all{n};
    
    in_high_ind   = find (table2array(trialinfo(:,'bini'))==1);
    in_low_ind    = find (table2array(trialinfo(:,'bini'))==2);
    in_con_ind    = find (table2array(trialinfo(:,'bini'))==3);
    out_high_ind  = find (table2array(trialinfo(:,'bini'))==4);
    out_low_ind   = find (table2array(trialinfo(:,'bini'))==5);
    out_con_ind   = find (table2array(trialinfo(:,'bini'))==6);
    
    % ttest zero
    data_in_high(n,:,:) = nanmean (trial_corr (in_high_ind,:,:));
    data_in_low(n,:,:)  = nanmean (trial_corr (in_low_ind,:,:));
    data_in_con(n,:,:)  = nanmean (trial_corr (in_con_ind,:,:));
    data_out_high(n,:,:)= nanmean (trial_corr (out_high_ind,:,:));
    data_out_low(n,:,:) = nanmean (trial_corr (out_low_ind,:,:));
    data_out_con(n,:,:) = nanmean (trial_corr (out_con_ind,:,:));
  
    % index all
    trial_ind_all{n, 1} = in_high_ind;
    trial_ind_all{n, 2} = in_low_ind;
    trial_ind_all{n, 3} = in_con_ind;
    trial_ind_all{n, 4} = out_high_ind;
    trial_ind_all{n, 5} = out_low_ind;
    trial_ind_all{n, 6} = out_con_ind;
    
end

data_in_hl = data_in_high - data_in_low;
data_in_hc = data_in_high - data_in_con;
data_in_cl = data_in_con - data_in_low;
data_out_hl = data_out_high - data_out_low;
data_out_hc = data_out_high - data_out_con;
data_out_cl = data_out_con - data_out_low;

data_comp_all = {data_in_high,data_in_low,data_in_con,data_out_high,data_out_low,data_out_con,...
    data_in_hl,data_in_hc,data_in_cl,data_out_hl,data_out_hc,data_out_cl};
contrast_type = {'ingroup-higher','ingroup-lower','ingroup-consistent',...
    'outgroup-higher','outgroup-lower','outgroup-consistent',...
    'ingroup-higher-vs-lower','ingroup-higher-vs-consistent','ingroup-consistent-vs-lower',...
    'outgroup-higher-vs-lower','outgroup-higher-vs-consistent','outgroup-consistent-vs-lower'};
contrast_ind = [0,0; 0,0; 0,0; 0,0; 0,0; 0,0;...
    1,2;...
    1,3;...
    3,2;...
    4,5;...
    4,6;...
    6,5];
for i_cont = 1:length(contrast_type)
    contrast_type{i_cont} = char(strcat(task, '-', contrast_type{i_cont}));
end

% pre data 

for n=1:numel(all_subs)
    
    trial_corr=trial_corr_pre{n};
    trialinfo =trialinfo_all{n};
    
    in_high_ind   = find (table2array(trialinfo(:,'bini'))==1);
    in_low_ind    = find (table2array(trialinfo(:,'bini'))==2);
    in_con_ind    = find (table2array(trialinfo(:,'bini'))==3);
    out_high_ind  = find (table2array(trialinfo(:,'bini'))==4);
    out_low_ind   = find (table2array(trialinfo(:,'bini'))==5);
    out_con_ind   = find (table2array(trialinfo(:,'bini'))==6);
    
    % ttest zero
    pre_data_in_high(n,:,:) = nanmean (trial_corr (in_high_ind,:,:));
    pre_data_in_low(n,:,:)  = nanmean (trial_corr (in_low_ind,:,:));
    pre_data_in_con(n,:,:)  = nanmean (trial_corr (in_con_ind,:,:));
    pre_data_out_high(n,:,:)= nanmean (trial_corr (out_high_ind,:,:));
    pre_data_out_low(n,:,:) = nanmean (trial_corr (out_low_ind,:,:));
    pre_data_out_con(n,:,:) = nanmean (trial_corr (out_con_ind,:,:));
  
end
data_comp_pre = {pre_data_in_high,pre_data_in_low,pre_data_in_con,...
    pre_data_out_high,pre_data_out_low,pre_data_out_con};

% post data

for n=1:numel(all_subs)
    
    trial_corr=trial_corr_post{n};
    trialinfo =trialinfo_all{n};
    
    in_high_ind   = find (table2array(trialinfo(:,'bini'))==1);
    in_low_ind    = find (table2array(trialinfo(:,'bini'))==2);
    in_con_ind    = find (table2array(trialinfo(:,'bini'))==3);
    out_high_ind  = find (table2array(trialinfo(:,'bini'))==4);
    out_low_ind   = find (table2array(trialinfo(:,'bini'))==5);
    out_con_ind   = find (table2array(trialinfo(:,'bini'))==6);
    
    % ttest zero
    post_data_in_high(n,:,:) = nanmean (trial_corr (in_high_ind,:,:));
    post_data_in_low(n,:,:)  = nanmean (trial_corr (in_low_ind,:,:));
    post_data_in_con(n,:,:)  = nanmean (trial_corr (in_con_ind,:,:));
    post_data_out_high(n,:,:)= nanmean (trial_corr (out_high_ind,:,:));
    post_data_out_low(n,:,:) = nanmean (trial_corr (out_low_ind,:,:));
    post_data_out_con(n,:,:) = nanmean (trial_corr (out_con_ind,:,:));
  
end
data_comp_post = {post_data_in_high,post_data_in_low,post_data_in_con,...
    post_data_out_high,post_data_out_low,post_data_out_con};

% statistical analysis

Stats_Results = struct;

for i_cont = 1:length(contrast_type)
    
    clear data h stat data_mask_neg data_L_neg data_num_neg data_negt ind_negt data;
    clear data_L_pos data_num_pos data_post ind_post;
    
    % find cluster in data
    data=data_comp_all{i_cont};

    [h,~,~,stat]=ttest(data,0,'Alpha',p_first);

    data_mask_neg=squeeze(h.*(stat.tstat<0));
    data_mask_neg(isnan(data_mask_neg))=0;
    [data_L_neg,data_num_neg] = bwlabel(data_mask_neg);
    for neg=1:data_num_neg
        m=find(data_L_neg==neg);
        data_negt(neg)=sum(stat.tstat(m));
    end
    if isempty(neg)
        data_negt=0;
    end
    [data_negt,ind_negt]=sort(data_negt,'ascend');

    data_mask_pos=squeeze(h.*(stat.tstat>0));
    data_mask_pos(isnan(data_mask_pos))=0;
    [data_L_pos,data_num_pos] = bwlabel(data_mask_pos);
    for pos=1:data_num_pos
        m=find(data_L_pos==pos);
        data_post(pos)=sum(stat.tstat(m));
    end
    if isempty(pos)
        data_post=0;
    end
    [data_post,ind_post]=sort(data_post,'descend');
    % keep tstat for later plotting
    data_stat=stat;

    % save data
    Stats_Results(i_cont).data = data;
    Stats_Results(i_cont).contrast = contrast_type{i_cont};
    Stats_Results(i_cont).stat = stat;
    Stats_Results(i_cont).data_negt = data_negt;
    Stats_Results(i_cont).ind_negt  = ind_negt;
    Stats_Results(i_cont).data_L_neg = data_L_neg;
    Stats_Results(i_cont).data_num_neg = data_num_neg;
    Stats_Results(i_cont).data_pos  = data_post;
    Stats_Results(i_cont).ind_post  = ind_post;
    Stats_Results(i_cont).data_L_pos = data_L_pos;
    Stats_Results(i_cont).data_num_pos = data_num_pos;
    
    clear data_num_neg data_num_pos ci ans min_neg min_pos neg_tsum pos_tsum
    
    % rand data
    null_mat=zeros(size(data));
    rand_data=null_mat;

    for r=1:rand
        
        for n=1:selvps
            
            trial_corr=trial_corr_all{n};
            
            if i_cont <= 6 % compare current matrix with zero
                
%                 trial_ind=trial_ind_all{n, i_cont};
%                 % randomly select all trials
%                 rand_ind =randperm(size(trial_corr,1), numel(trial_ind));
%                 rand_data(n,:,:) = nanmean(trial_corr (rand_ind,:,:));             
                dat_1 = data_comp_pre{i_cont};
                dat_2 = data_comp_post{i_cont};
                dat_1 = dat_1(n,:,:);
                dat_2 = dat_2(n,:,:);
                rand_ind = randi([0 1], size(dat_1));
                dat_1_rand = rand_ind.*dat_1 + (~rand_ind).*dat_2;
                dat_2_rand = rand_ind.*dat_2 + (~rand_ind).*dat_1;
                rand_data(n,:,:) = dat_1_rand - dat_2_rand;
                
            elseif i_cont > 6 % switch two conditions
                
                dat_1 = data_comp_all{contrast_ind(i_cont,1)};
                dat_2 = data_comp_all{contrast_ind(i_cont,2)};
                dat_1 = dat_1(n,:,:);
                dat_2 = dat_2(n,:,:);
                rand_ind = randi([0 1], size(dat_1));
                dat_1_rand = rand_ind.*dat_1 + (~rand_ind).*dat_2;
                dat_2_rand = rand_ind.*dat_2 + (~rand_ind).*dat_1;
                rand_data(n,:,:) = dat_1_rand - dat_2_rand;
             
            end 
            
        end

        [h,~,~,stat]=ttest(rand_data,0,'Alpha',p_first);

        rand_mask_neg=squeeze(h.*(stat.tstat<0));

        rand_mask_neg(isnan(rand_mask_neg))=0;
        [rand_L_neg,rand_num_neg] = bwlabel(rand_mask_neg);

        for neg=1:rand_num_neg
            m=find(rand_L_neg==neg);
            rand_negt(neg)=sum(stat.tstat(m));
        end
        if isempty(neg)
            rand_negt=0;
        end

        rand_mask_pos=squeeze(h.*(stat.tstat>0));
        rand_mask_pos(isnan(rand_mask_pos))=0;
        [rand_L_pos,rand_num_pos] = bwlabel(rand_mask_pos);
        for pos=1:rand_num_pos
            m=find(rand_L_pos==pos);
            rand_post(pos)=sum(stat.tstat(m));
        end
        if isempty(pos)
            rand_post=0;
        end
        pos_tsum{r}=sort(rand_post,'descend');
        neg_tsum{r}=sort(rand_negt,'ascend');

    end
    
    max_pos=max(cellfun(@numel,pos_tsum));
    max_neg=max(cellfun(@numel,neg_tsum));

    for x=1:rand
        pos_tsum{x}= [pos_tsum{x},zeros(1,max_pos-numel(pos_tsum{x}))];
        neg_tsum{x}= [neg_tsum{x},zeros(1,max_neg-numel(neg_tsum{x}))];
    end

    pos_dist=reshape([pos_tsum{:}],max_pos,rand);
    pos_dist=sort(pos_dist,2,'descend');

    neg_dist=reshape([neg_tsum{:}],max_neg,rand);
    neg_dist=sort(neg_dist,2,'ascend');
    
    Stats_Results(i_cont).pos_tsum = pos_tsum;
    Stats_Results(i_cont).neg_tsum = neg_tsum;
    
    clear rand_L_neg rand_L_pos rand_num_neg rand_num_pos ci ans  neg_tsum pos_tsum
    clear rand_data rand_def_vec rand_inhibtion rand_mask_neg rand_mask_pos
    clear rand_negt rand_post rand_rehearsal rep rand_vec pos neg

    % significance?
    % check pos clusters
    ind=1;
    p_threshold=p_cluster;
    sig_mask_pos=zeros(size(data_mask_pos));
    pos_cluster_p=[];
    pos_cluster(ind)=nearest(pos_dist(ind,:),data_post(ind))./rand;
    while nearest(pos_dist(ind,:),data_post(1))<=round(p_threshold.*rand)
        p_threshold=p_threshold/2;
        sig_mask_pos=sig_mask_pos+(data_L_pos==ind_post(ind));
        pos_cluster(ind)=nearest(pos_dist(ind,:),data_post(ind))./rand;
        ind=ind+1;
        if ind>numel(data_post)|ind> max_pos
            break
        end
    end

    % check neg clusters
    ind=1;
    p_threshold=p_cluster;
    sig_mask_neg=zeros(size(data_mask_neg));
    neg_cluster_p=[];
    neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand;
    while nearest(neg_dist(ind,:),data_negt(ind))<=round(p_threshold.*rand)
        p_threshold=p_threshold/2;
        sig_mask_neg=sig_mask_neg+(data_L_neg==ind_negt(ind));
        neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand;
        ind=ind+1;
        if ind>numel(data_negt)|ind> max_neg
            break
        end

    end
    
    clear ind p_threshold
    
    % plot results
    mask_alpha=sig_mask_pos+sig_mask_neg;

    contrast = contrast_type{i_cont};
    cmp = customcolormap_preset('red-white-blue');
    
    corr_trials = pre_data.corr_trials;
    
    fig =  figure
    subplot(2,2,1:2)
    H = imagesc(mean(corr_trials.time_exp),mean(corr_trials.time_mor),squeeze(data_stat.tstat),[-5 5])
    title(strcat(contrast))
    colormap(cmp)
    c = colorbar;
    c.Label.String = 't-value'
    set(gca,'YDir','normal')
    %set(H,'AlphaData',mask_alpha)
    hold on
    contour(mean(corr_trials.time_mor),mean(corr_trials.time_exp),mask_alpha,1,'LineColor','k','LineWidth',1)

    subplot(2,2,3)
    hist(pos_dist(1,:))
    hold on
    plot([data_post(1) data_post(1)],[0 rand],'r')
    title('distribution positive cluster')
    text(data_post(1),rand-100,strcat('pcorr=',num2str(pos_cluster)))

    subplot(2,2,4)
    hist(neg_dist(1,:))
    hold on
    plot([data_negt(1) data_negt(1)],[0 rand],'r')
    title('distribution negative cluster')
    text(data_negt(1), rand-100,strcat('pcorr=',num2str(neg_cluster)))
    savefig(fig,fullfile(path_out,strcat('clusterstat_alltrialscorr','_',contrast)))

    Stats_Results(i_cont).pos_dist_rand = pos_dist;
    Stats_Results(i_cont).neg_dist_rand = neg_dist;
    Stats_Results(i_cont).sig_mask_pos_rand = sig_mask_pos;
    Stats_Results(i_cont).sig_mask_neg_rand = sig_mask_neg;
    Stats_Results(i_cont).timex = corr_trials.time_exp;
    Stats_Results(i_cont).timey = corr_trials.time_mor;
    Stats_Results(i_cont).pneg = neg_cluster;
    Stats_Results(i_cont).ppos = pos_cluster;

    
    close all
    
end

cd (path_out);
save(fullfile(path_out,strcat('clusterstat_alltrialscorr')),'Stats_Results');
