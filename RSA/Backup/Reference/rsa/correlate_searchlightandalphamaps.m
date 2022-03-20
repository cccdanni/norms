project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%%
% correlate searchlight rsa with alpha power effects


effect='rehearsal'; % or 'inhibition'
timing='late'; % or 'early'

path_alpha=fullfile(project_folder,'source','source_statplots','lcmv10');

if strcmp(effect, 'inhibition')
path_rsa=fullfile(project_folder,'data,','source','lcmv10','early_clustersphere64_all_trials_item_cue_sr100window200slide10_sourceAll_zall','stats');
file_rsa='stat_inhibition_clusteralpha0.05.mat';
    if strcmp(timing, 'early')
    file_alpha='stat_inhibition_clusteralpha0.018_13_hz2000_3500.mat';
    elseif strcmp(timing, 'late')
    file_alpha='stat_inhibition_clusteralpha0.018_13_hz2500_3000.mat';
    end
elseif strcmp(effect, 'rehearsal')
path_rsa=fullfile(project_folder,'data,','source','lcmv10','late_clustersphere64_all_trials_item_cue_sr100window200slide10_sourceAll_zall','stats');
 file_rsa='stat_rehearsal_clusteralpha0.05.mat';
    if strcmp(timing, 'early')
    file_alpha='stat_rehearsal_clusteralpha0.018_13_hz3000_3500.mat';
    elseif strcmp(timing, 'late')
    file_alpha='stat_rehearsal_clusteralpha0.018_13_hz3500_4000.mat';
    end
end

load(fullfile(path_rsa,file_rsa))
stat_rsa=stat;

load(fullfile(path_alpha,file_alpha))
stat_alpha=stat;

% correlate values
    %remove Nans
    sel_ind=~isnan(stat_alpha.stat);
    sel_ind_rsa=~isnan(stat_rsa.stat);
    % check, both sel ind should be same
    check=sum(sel_ind-sel_ind_rsa)==0    

[rho,p]=corr(stat_alpha.stat(sel_ind),stat_rsa.stat(sel_ind), 'Type','Spearman'	)

sel_ind=find(sel_ind);
% random rho distribution
for r=1:1000
    rand_ind=sel_ind(randperm(numel(sel_ind)));
[rho_rand(r)]=corr(stat_alpha.stat(sel_ind),stat_rsa.stat(rand_ind), 'Type','Spearman'	);

end
if rho<0
    pcorr=sum(rho_rand<rho)/1000;
elseif rho>0
    pcorr=sum(rho_rand>rho)/1000;
end

figure
scatterhistogram(stat_alpha.stat,stat_rsa.stat,'HistogramDisplayStyle','smooth')
title(strcat('alpha-searchlight map correlation, r=,',num2str(rho),', pcorr=',num2str(pcorr)))
xlabel('alpha')
ylabel('rsa')

% median
sel_ind=~isnan(stat_alpha.stat);
mean_rsa=mean(stat_rsa.stat(sel_ind));
mean_alpha=mean(stat_alpha.stat(sel_ind));

vec_up=0:1/90:(1-1/90);
vec_down=1-vec_up;
circmap=[ones(90,1),vec_down',vec_down';...
        vec_down',zeros(90,1),zeros(90,1);...
        zeros(90,1),zeros(90,1),vec_up';...
        vec_up',vec_up',ones(90,1)];
colormap(circmap)
figure
imagesc(1:1:360)
colormap(circmap)

% get radians of each voxel pair 
vec_rsa_alpha=complex(stat_alpha.stat-mean_alpha, stat_rsa.stat-mean_rsa);
angle_rsa_alpha=rad2deg(angle(vec_rsa_alpha));%-deg2rad(135);
stat_angle=stat_rsa;
stat_angle.stat=((angle_rsa_alpha))+(((angle_rsa_alpha)<0).*360);
save(fullfile(path_rsa,strcat('angle_',timing,file_rsa)),'stat_angle')

figure
scatterhist(stat_alpha.stat,stat_rsa.stat,'Location','NorthEast','Kernel','on')
title(strcat('alpha-searchlight map correlation, r=,',num2str(rho),', pcorr=',num2str(pcorr)))
xlabel('alpha')
ylabel('rsa')

if strcmp(effect, 'rehearsal')
circmap=[circmap(end-134:end,:);circmap(1:225,:)]
elseif strcmp(effect, 'inhibition')
    circmap=[circmap(end-44:end,:);circmap(1:315,:)]

end
    figure
imagesc(1:360)
colormap(circmap)

figure
angle2color=(ceil(angle_rsa_alpha(sel_ind)))+((ceil(angle_rsa_alpha(sel_ind))<0).*360)+1;
sel_color=circmap(angle2color,:);
scatter(stat_alpha.stat(sel_ind)-mean_alpha,stat_rsa.stat(sel_ind)-mean_rsa,[],sel_color,'filled','MarkerEdgeColor',[0 0 0])
title(strcat('alpha-searchlight map correlation, r=,',num2str(rho),', pcorr=',num2str(pcorr)))
xlabel('alpha')
ylabel('rsa')


figure
angle2color=(ceil(angle_rsa_alpha(sel_ind)))+((ceil(angle_rsa_alpha(sel_ind))<0).*360)+1;
sel_color=circmap(angle2color,:);
if strcmp(effect, 'rehearsal')
mask=find(angle2color<180|angle2color>270);
elseif strcmp(effect, 'inhibition')
mask=find(angle2color<270);

end
sel_color(mask,:)=repmat([1 1 1],numel(mask),1);
scatter(stat_alpha.stat(sel_ind),stat_rsa.stat(sel_ind),[],sel_color,'filled','MarkerEdgeColor',[0 0 0])
title(strcat('alpha-searchlight map correlation, r=,',num2str(rho),', pcorr=',num2str(pcorr)))
xlabel('alpha')
ylabel('rsa')


 
%% correlate different alphas with each other for statistical comparision of correlations (test using cocorr)

path_alpha=fullfile(project_folder,'source','source_statplots','lcmv10');
file_alpha='stat_inhibition_clusteralpha0.018_13_hz2000_2500.mat';
% path_alpha=fullfile(project_folder,'source','source_statplots','lcmv10');
% file_alpha='stat_rehearsal_clusteralpha0.018_13_hz3000_3500.mat';
load(strcat(path_alpha,file_alpha))
stat_alpha1=stat;

path_alpha=fullfile(project_folder,'source','source_statplots','lcmv10');
file_alpha='stat_inhibition_clusteralpha0.018_13_hz2500_3000.mat';
% path_alpha=fullfile(project_folder,'source','source_statplots','lcmv10');
% file_alpha='stat_rehearsal_clusteralpha0.018_13_hz3500_4000.mat';
load(strcat(path_alpha,file_alpha))
stat_alpha2=stat;

path_rsa=fullfile(project_folder,'source','lcmv10_rsa','early_clustersphere64_all_trials_item_cue_sr100window200slide10_sourceAll_zall','stats');
file_rsa='stat_inhibition_clusteralpha0.05.mat';
% path_rsa=fullfile(project_folder,'source','lcmv10_rsa','late_clustersphere64_all_trials_item_cue_sr100window200slide10_sourceAll_zall','stats');
% file_rsa='stat_rehearsal_clusteralpha0.05.mat';
load(strcat(path_rsa,file_rsa))
stat_rsa=stat;

sel_ind=~isnan(stat_alpha1.stat);
[rho_alphas,p]=corr(stat_alpha1.stat(sel_ind),stat_alpha2.stat(sel_ind))
[rho_alpha1,p]=corr(stat_alpha1.stat(sel_ind),stat_rsa.stat(sel_ind))
[rho_alpha2,p]=corr(stat_alpha2.stat(sel_ind),stat_rsa.stat(sel_ind))


%% surface plots of results
effect='inhibition'
load(fullfile(project_folder,'scripts','additional_functions','circmap_black2blue.mat'))

if strcmp(effect,'inhibition')
path_stats=fullfile(project_folder,'source','lcmv10_rsa','early_clustersphere64_all_trials_item_cue_sr100window200slide10_sourceAll_zall','stats');
circmap=[circmap(end-44:end,:);circmap(1:315,:)];
elseif strcmp(effect,'rehearsal')
path_stats=fullfile(project_folder,'source','lcmv10_rsa','late_clustersphere64_all_trials_item_cue_sr100window200slide10_sourceAll_zall','stats');
circmap=[circmap(end-134:end,:);circmap(1:225,:)]
end


dir_stat=dir(path_stats)

files_ind=find(strncmp('angle',{dir_stat(:).name},4));
stat_list={dir_stat(files_ind).name};
nsub=18; % for t value threshold
%
path_out=fullfile(path_stats,'figures');
mkdir(path_out)
mri=ft_read_mri(fullfile(toolbox_folder, 'fieldtrip-20160122','template','anatomy','single_subj_T1_1mm.nii'));

hemisphere={'left','right'};

views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90];


for s=1:numel(stat_list)
        load(strcat(path_stats,stat_list{s}))
        stat=stat_angle;
        statmask=stat;
         cfg = [];
       cfg.parameter = 'stat';
       cfg.interpmethod  = 'nearest';

      statint=ft_sourceinterpolate(cfg, stat,mri);
     if strcmp(effect,'inhibition')
      statmask.stat=stat.stat>270;
     elseif strcmp(effect, 'rehearsal')
        statmask.stat= stat.stat>180&stat.stat<270;
     end
      statmask=ft_sourceinterpolate(cfg, statmask,mri); 
      statint.mask=statmask.stat~=0;
      
for h=1:numel(hemisphere)
    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'stat';
    cfg.maskparameter  = 'mask';
    cfg.funcolorlim    = [0 360];
    cfg.funcolormap    =circmap;
    cfg.projmethod     = 'nearest';
    cfg.surffile       = fullfile(project_folder,'scripts','additional_functions',strcat('carret',hemisphere{h},'.mat'));
    cfg.surfdownsample = 5;
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

view(squeeze(views(h,3,:))');
print('-f1','-dtiff','-r600',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_lat_mask.tiff')) 

view(squeeze(views(h,4,:))');
print('-f1','-r600','-dtiff',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_med_mask.tiff')) 

view(squeeze(views(h,2,:))');
print('-f1','-r600','-dtiff',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_vent2_mask.tiff')) 
clear c1 c2 

close all
    end

end