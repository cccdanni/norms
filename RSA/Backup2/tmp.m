
project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));


%% RSA for DF: correlate enco with enco
path_in=fullfile(project_folder,'data');

% define the participants - 
all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

cd(path_in);
% define sliding window
win=0.2; % in sec
slide=0.01;
% define new Sampling rate
sr=100; % in Hz
step=win/(1/sr)+1;

path_out = fullfile(project_folder,'RSA','data',strcat('all_trials_item_cue_sr',num2str(sr),'window',num2str(win.*1000),'slide',num2str(slide.*1000)));
mkdir(path_out)

channels='all';

%enco: define start and end of item window
t_start_e=-0.1;
t_end_e=0.7;
tois_e=t_start_e:1/sr:t_end_e;
t1_e=t_start_e:slide:(t_end_e-win);
t2_e=t1_e+win;
ind_t1_e=1:slide/(1/sr):((numel(tois_e)-win/(1/sr)));
ind_t2_e=ind_t1_e+win/(1/sr);
n_bins_e=numel(t1_e);

%cue define start and end of item window
t_start_c=1.9;
t_end_c=4.1;
tois_c=t_start_c:1/sr:t_end_c;
t1_c=t_start_c:slide:(t_end_c-win);
t2_c=t1_c+win;
ind_t1_c=1:slide/(1/sr):((numel(tois_c)-win/(1/sr)));
ind_t2_c=ind_t1_c+win/(1/sr);
n_bins_c=numel(t1_c);

% first step: correlate for each subject in all trials
block=[11,13]; % for enco

for n=1:numel(all_subs)
    sel_sub=all_subs{n};
    %load in the preprocessed data of each participant
    load(fullfile(path_in,strcat(sel_sub,'_data')));
    
    cfg=[];
    cfg.lpfilter='yes';
    cfg.lpfreq=40;
    % for bs filter
    %     cfg.bsfilter='yes';
    %     cfg.bsfilter=[7 20];
    data=ft_preprocessing(cfg,data);
    
    cfg=[];
    cfg.resamplefs=sr;
    cfg.detrend='no';
    data=ft_resampledata(cfg,data);
    cfg=[];
    cfg.channel=channels;
    data=ft_preprocessing(cfg,data);
    
    % find trial in enco & cue/reco
    enco_trials=find(data.trialinfo(:,5)==13 |data.trialinfo(:,5)==11);
    cue_trials= find(data.trialinfo(:,5)==block(1) |data.trialinfo(:,5)==block(2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % prepare enco data
    cfg=[];
    cfg.latency=[t_start_e t_end_e+1/sr];
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
    cfg.latency=[t_start_c t_end_c+1/sr];
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
    %
    [items, use_enco, use_cue]=intersect(dataenc.trialinfo(:,7),datacue.trialinfo(:,7));
    
    % resort trials/trialsinfo for matching trial order
    dataenc.trial=dataenc.trial(use_enco,:,:);
    dataenc.trialinfo=dataenc.trialinfo(use_enco,:);
    datacue.trial=datacue.trial(use_cue,:,:);
    datacue.trialinfo=datacue.trialinfo(use_cue,:);
    
    trialinfo=dataenc.trialinfo; %same as for datacue
    
    data_enc_vec=zeros(size(trialinfo,1),n_bins_e,numel(dataenc.label)*step);
    for bin=1:n_bins_e
        % vectorize sel_bins: data_vec(trials, nbins, features)
        data_vec_tmp=dataenc.trial(:,:,ind_t1_e(bin):ind_t2_e(bin));
        data_enc_vec(:,bin,:)=reshape(data_vec_tmp,size(trialinfo,1),[]);
    end
    clear data_vec_tmp dataenc
    
    data_cue_vec=zeros(size(trialinfo,1),n_bins_c,numel(datacue.label)*step);
    for bin=1:n_bins_c
        % vectorize sel_bins: data_vec(trials, nbins, features)
        data_vec_tmp=datacue.trial(:,:,ind_t1_c(bin):ind_t2_c(bin));
        data_cue_vec(:,bin,:)=reshape(data_vec_tmp,size(trialinfo,1),[]);
    end
    clear data_vec_tmp datacue
    
    data_enc_vec2=reshape(data_enc_vec, size(data_enc_vec,1)*size(data_enc_vec,2),size(data_enc_vec,3));
    data_cue_vec2=reshape(data_cue_vec, size(data_cue_vec,1)*size(data_cue_vec,2),size(data_cue_vec,3));
    corr_tmp=corr(data_enc_vec2', data_cue_vec2', 'type','Spearman');
    
    % size corr_cue_enc= trials enco x bins enco x trials cue x bins cue
    corr_cue_enc=  reshape(corr_tmp,size(data_enc_vec,1),size(data_enc_vec,2),size(data_cue_vec,1),size(data_cue_vec,2));
    
    clear corr_tmp data data_enc_vec data_cue_vec data_enc_vec2 data_cue_vec2
    
    % fisher z transform correlations
    corr_cue_enc=  0.5.*log(((ones(size(corr_cue_enc))+corr_cue_enc)./(ones(size(corr_cue_enc))-corr_cue_enc)));
    
    % prepare data: get on diagonal i.e only in trial correlations
    for tr=1:size(trialinfo,1)
        corr_cue_enc_trial(tr,:,:)= squeeze(corr_cue_enc(tr,:,tr,:));
    end
    
    corr_trials.corr_cue_enc_trial=corr_cue_enc_trial;
    corr_trials.trialinfo=trialinfo;
    corr_trials.time_item=[t1_e;t2_e];
    corr_trials.time_cue=[t1_c;t2_c];
    
    save(fullfile(path_out, strcat(all_subs{n},'item_cue_alltrials.mat')),'corr_trials');
    
    clear corr_trials corr_cue_enc_trial trialinfo corr_cue_enc items
end


%% load trial correlations and draw random in each subject
path_out = fullfile(project_folder,'RSA','data','all_trials_item_cue_sr100window200slide10');

all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

contrast='interaction';
load(fullfile(project_folder,'scripts','additional_functions','jet_grey2.mat'))

rand=1000;
p_cluster=0.05;
p_first=0.01;
% combine all vps in o
for n=1:numel(all_subs)
    load(fullfile(path_in, strcat(all_subs{n},'item_cue_alltrials.mat')))
    trial_corr_all{n}=corr_trials.corr_cue_enc_trial;
    trialinfo_all{n}=corr_trials.trialinfo;
end

% data stats
for n=1:numel(all_subs)
    trial_corr=trial_corr_all{n};
    trialinfo=trialinfo_all{n};
    tbr_r_ind=trialinfo(:,5)==11&trialinfo(:,10)==1;
    tbr_f_ind=trialinfo(:,5)==11&trialinfo(:,10)==0;
    
    tbf_r_ind=trialinfo(:,5)==13&trialinfo(:,10)==1;
    tbf_f_ind=trialinfo(:,5)==13&trialinfo(:,10)==0;
    
    trial_def_vec=tbr_r_ind+(tbr_f_ind.*2)+(tbf_r_ind.*3)+(tbf_f_ind.*4);
    
    data_inhibition(n,:,:)=nanmean(trial_corr(trial_def_vec==4,:,:))-nanmean(trial_corr(trial_def_vec==2,:,:));
    data_rehearsal(n,:,:)=nanmean(trial_corr(trial_def_vec==3,:,:))-nanmean(trial_corr(trial_def_vec==1,:,:));
    
    trial_def_vec_all{n}=trial_def_vec;
end


% find cluster in data
data=data_inhibition-data_rehearsal;

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
clear  data_num_neg data_num_pos ci ans min_neg min_pos neg_tsum pos_tsum
clear tbr_r_ind tbf_r_ind tbr_f_ind tbf_f_ind
% rand data


null_mat=zeros(size(data_inhibition));
rand_inhibition=null_mat;
rand_rehearsal=null_mat;

for r=1:rand
    for n=1:numel(selvps)
        trial_corr=trial_corr_all{n};
        trial_def_vec=trial_def_vec_all{n};
        
        % two different randoptions: all trials random, or in conditions random
        % switch all trials
        rand_vec=randperm(numel(trial_def_vec));
        rand_def_vec=trial_def_vec(rand_vec);
        
        rand_inhibition(n,:,:)=nanmean(trial_corr(rand_def_vec==4,:,:))-nanmean(trial_corr(rand_def_vec==2,:,:));
        rand_rehearsal(n,:,:)=nanmean(trial_corr(rand_def_vec==3,:,:))-nanmean(trial_corr(rand_def_vec==1,:,:));
    end
    
    rand_data=rand_inhibition-rand_rehearsal;
    
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

fig= figure
subplot(2,2,1:2)
H= imagesc(mean(corr_trials.time_cue),mean(corr_trials.time_item),squeeze(data_stat.tstat),[-5 5])
title(strcat(contrast))
colormap(jet_grey)
set(gca,'YDir','normal')
%set(H,'AlphaData',mask_alpha)
hold on
contour(mean(corr_trials.time_cue),mean(corr_trials.time_item),mask_alpha,1,'LineColor','k','LineWidth',1)


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
savefig(fig,fullfile(path_in,strcat('clusterstat_alltrialscorr','_',contrast)))

clusterstat.stat=data_stat.tstat;
clusterstat.timex=corr_trials.time_cue;
clusterstat.timey=corr_trials.time_item;
clusterstat.maskpos=sig_mask_pos;
clusterstat.maskneg=sig_mask_neg;
clusterstat.pneg=neg_cluster;
clusterstat.ppos=pos_cluster;
clusterstat.negdist=neg_dist;
clusterstat.posdist=pos_dist;
save(fullfile(path_in,strcat('clusterstat_alltrialscorr','_',contrast)),'clusterstat')
close all
clear


%% plot data out of clusters

% load clusterstats
path_in = fullfile(project_folder,'RSA','data','all_trials_item_cue_sr100window200slide10');
all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

contrast='interaction';
load(fullfile(project_folder,'scripts','additional_functions','jet_grey2.mat'))

conditions={'TBF-F','TBR-F','TBF-R','TBR-R'};
cond_def=[13,0;11,0;13,1;11,1];
load(fullfile(path_in,strcat('clusterstat_alltrialscorr','_',contrast,'.mat')))
neg_check=0;
pos_check=0;

% check for sig cluster
if any(clusterstat.pneg<0.05)
    neg_cluster_def=find(clusterstat.maskneg);
    neg_check=1;
end

if any(clusterstat.ppos<0.05)
    pos_cluster_def=find(clusterstat.maskpos);
    pos_check=1;
end

for n=1:numel(all_subs)
    load(fullfile(path_in,strcat(all_subs{n},'item_cue_alltrials.mat')))
        for c=1:numel(conditions)
        
        if~isfield(corr_trials,'corr_cue_enc_trial')
            cond_order=[4,2,3,1];
            tmp =squeeze(corr_trials.corr_cue_enc_cond(cond_order(c),:,:));
        else
            sel_ind=corr_trials.trialinfo(:,5)==cond_def(c,1)&corr_trials.trialinfo(:,10)==cond_def(c,2);
            tmp =squeeze(nanmean(corr_trials.corr_cue_enc_trial(sel_ind,:,:)));
        end
        if pos_check
            pos_cluster_data(n,c)=nanmean(tmp(pos_cluster_def));        
        end
        
        if neg_check
            neg_cluster_data(n,c)=nanmean(tmp(neg_cluster_def));
        end
    end
end


%% plots
fig1=figure


% plot all bars mean centered for f and r
neg_cluster_zdata(:,1:2)=neg_cluster_data(:,1:2)-repmat(mean(neg_cluster_data(:,1:2),2),numel(selvps,1));
neg_cluster_zdata(:,3:4)=neg_cluster_data(:,3:4)-repmat(mean(neg_cluster_data(:,3:4),2),numel(selvps,1));
pos_cluster_zdata(:,1:2)=pos_cluster_data(:,1:2)-repmat(mean(pos_cluster_data(:,1:2),2),numel(selvps,1));
pos_cluster_zdata(:,3:4)=pos_cluster_data(:,3:4)-repmat(mean(pos_cluster_data(:,3:4),2),numel(selvps,1));

% plot early and late cluster separately
fig=figure
neg_inhibition=neg_cluster_zdata(:,1)-neg_cluster_zdata(:,2);
neg_rehearsal=neg_cluster_zdata(:,3)-neg_cluster_zdata(:,4);

figure
subplot(2,2,1)
bar([nanmean(neg_cluster_zdata(:,1:2)),nanmean(pos_cluster_zdata(:,1:2))])
hold on
scatter([ones(1,numel(selvps)),ones(1,numel(selvps)).*2,ones(1,numel(selvps)).*3,ones(1,numel(selvps)).*4]...
    ,reshape([neg_cluster_zdata(:,1:2),pos_cluster_zdata(:,1:2)],[],1),'ok')
xticklabels([conditions(1:2),conditions(1:2)])
ylim([-0.06,0.06])
title(strcat('inhibition:',contrast))
subplot(2,2,2)
axis off
[h,p,~,stat]=ttest(neg_inhibition);
text(0,1,strcat('inhibition early: t(',num2str(stat.df),')=',num2str(stat.tstat),', p=',num2str(p) ))
[h,p,~,stat]=ttest(pos_inhibition);
text(0,0.8,strcat('inhibition late: t(',num2str(stat.df),')=',num2str(stat.tstat),', p=',num2str(p) ))


subplot(2,2,3)
bar([nanmean(neg_cluster_zdata(:,3:4)),nanmean(pos_cluster_zdata(:,3:4))])
hold on
scatter([ones(1,numel(selvps)),ones(1,numel(selvps)).*2,ones(1,numel(selvps)).*3,ones(1,numel(selvps)).*4]...
    ,reshape([neg_cluster_zdata(:,3:4),pos_cluster_zdata(:,3:4)],[],1),'ok')
xticklabels([conditions(3:4),conditions(3:4)])
ylim([-0.06,0.06])
title(strcat('rehearsal:',contrast))
subplot(2,2,4)
axis off
[h,p,~,stat]=ttest(neg_rehearsal);
text(0,1,strcat('rehearsal early: t(',num2str(stat.df),')=',num2str(stat.tstat),', p=',num2str(p) ))
[h,p,~,stat]=ttest(pos_rehearsal);
text(0,0.8,strcat('rehearsal late: t(',num2str(stat.df),')=',num2str(stat.tstat),', p=',num2str(p) ))
savefig(fig,fullfile(path_in,'bargraphs_meanconditions_effectorder'))

%% correlations inhibition and rehearsal

rsa_inhibition=neg_cluster_data(:,1)-neg_cluster_data(:,2);
rsa_rehearsal=pos_cluster_data(:,3)-pos_cluster_data(:,4);

% check for correlation inhibition/rehearsal effect: rsa
[rho,p]=corr(rsa_inhibition,rsa_rehearsal,'Type','Spearman')
