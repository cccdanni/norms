
project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%% conjunction analysis of inhibition & rehearsal contrast
path_out = fullfile(project_folder,'RSA','data','all_trials_item_cue_sr100window200slide10');
all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

load(fullfile(project_folder,'scripts','additional_functions','jet_grey_halfmax.mat'))

rand=1000;
p_first=0.05;
% combine all vps in o
for n=1:numel(all_subs)
    load(strcat(path_in, all_subs{n},'item_cue_alltrials'))
    
    trialinfo=corr_trials.trialinfo;
    tbr_r_ind=trialinfo(:,5)==11&trialinfo(:,10)==1;
    tbr_f_ind=trialinfo(:,5)==11&trialinfo(:,10)==0;
    
    tbf_r_ind=trialinfo(:,5)==13&trialinfo(:,10)==1;
    tbf_f_ind=trialinfo(:,5)==13&trialinfo(:,10)==0;
    
    % 'inhibition'
    data1(n,:,:)=nanmean(corr_trials.corr_cue_enc_trial(tbf_f_ind,:,:))-nanmean(corr_trials.corr_cue_enc_trial(tbr_f_ind,:,:));
    % 'rehearsal'
    data2(n,:,:)=nanmean(corr_trials.corr_cue_enc_trial(tbf_r_ind,:,:))-nanmean(corr_trials.corr_cue_enc_trial(tbr_r_ind,:,:));
    
end



% conjunction on t value basis
[h1,p1,~,stat1]=ttest(data1,0,'Alpha',p_first);
[h2,p2,~,stat2]=ttest(data2,0,'Alpha',p_first);
crit005 = tinv(1 - (p_first/2), numel(selvps)-1);
crit01= tinv(1 - (0.05), numel(selvps)-1);
crit02= tinv(1 - (0.1), numel(selvps)-1);

% conjunction t_map:
% sign(stat1.tstat)==sign(stat2.tstat) & min(abs(stat1.stat),
sign_mask=sign(stat1.tstat)==sign(stat2.tstat);
min_t = min([abs(stat1.tstat);abs(stat2.tstat)],[],1).*(sign_mask.*sign(stat1.tstat));
max_p= max([p1;p2],[],1).*(sign_mask);
max_p(max_p==0)=1;

tcdf(min_t,numel(selvps)-1)
figure
subplot(3,4,[1 2 3 5 6 7 9 10 11])
imagesc(mean(corr_trials.time_cue),mean(corr_trials.time_item),squeeze(min_t),[-5 5])
title(['inhibition & rehearsal conjunction: absolute min t-value & same direction of t-vales'])
colormap(gca,jet_grey)
set(gca,'YDir','normal')
colorbar
subplot(3,4,4)
imagesc(mean(corr_trials.time_cue),mean(corr_trials.time_item),squeeze(max_p),[0 1])
title('threshold p<0.05')
colormap(gca,'hot')
set(gca,'YDir','normal')
colorbar
