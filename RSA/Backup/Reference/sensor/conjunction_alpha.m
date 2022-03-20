project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%% conjunction alpha effects

load(fullfile(project_folder,'scripts','additional_functions','jet_grey_halfmax.mat'))
% load stat maps


% load stat rehearsal
load(fullfile(project_folder,'data','TF','statistics','3d_GA_z_enc_TBF_r_GA_z_enc_TBR_r_2_30_2000_4000_clusterp0.01.mat'))
stat_rehearsal=stat;
% load stat inhibition
load(fullfile(project_folder,'data','TF','statistics','3d_GA_z_enc_TBF_f_GA_z_enc_TBR_f_2_30_2000_4000_clusterp0.01.mat'))
stat_inhibition=stat;


stat1=stat_inhibition.stat;
stat2=stat_rehearsal.stat;

% conjunction t_map:
% sign(stat1.tstat)==sign(stat2.tstat) & min(abs(stat1.stat),
sign_mask=sign(stat1)==sign(stat2);
stat_all(1,:,:,:)=stat1;
stat_all(2,:,:,:)=stat2;
min_t = squeeze(min(abs(stat_all),[],1)).*(sign_mask.*sign(stat1));

figure
imagesc(squeeze(mean(min_t)));

imagesc(squeeze(sum(abs(min_t)>2.1098)))


foi=[8 13];
f1=nearest(stat.freq,foi(1));
f2=nearest(stat.freq,foi(2));

% sort electrodes for plots: first right/left (dim1) then anterior
% posterior (dim2)

[sorted,ind]= sortrows(elec_mnivol.chanpos,[2,1],'descend');
[sorted_label,indelec,indstat]=intersect(elec_mnivol.label(ind),stat.label,'stable');


figure
subplot(1,3,1)
imagesc(stat.time, 1:64, squeeze(nanmean(stat1(indstat,f1:f2,:),2)),[-4 4])
colormap(jet_grey2)
title('inhibition')
colorbar
subplot(1,3,2)
imagesc(stat.time, 1:64, squeeze(nanmean(stat2(indstat,f1:f2,:),2)),[-4 4])
colormap(jet_grey2)
title('rehearsal')
colorbar
subplot(1,3,3)
imagesc(stat.time, 1:64, squeeze(nanmean(min_t(indstat,f1:f2,:),2)),[-4 4])
colormap(jet_grey2)
colorbar

title('conjunction (minimum t): inhibition & rehearsal')
