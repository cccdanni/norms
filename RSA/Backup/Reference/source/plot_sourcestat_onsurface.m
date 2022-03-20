project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%% gathers all source stats
path_stats=(fullfile(project_folder,'source','source_statplots',beam_type));
dir_stat=dir(path_stats)

files_ind=find(strncmp('stat',{dir_stat(:).name},4));
stat_list={dir_stat(files_ind).name};
nsub=18; % for t value threshold
%
path_out=fullfile(path_stats,'figures');
mkdir(path_out)
mri=ft_read_mri(fullfile(toolbox_folder,'SPM8','spm8','templates','T1.nii'))

hemisphere={'left','right'};

views(1,:,:)=[-90,30;90 -30;-90,0;90,0;0,-90];
views(2,:,:)=[90,30;-90 -30;90,0;-90,0;0,-90];

% mask with half maximum
load(fullfile(project_folder,'scripts','additional_functions','jet_grey_halfmax.mat'))

for s=1:numel(stat_list)
        load(strcat(path_stats,stat_list{s}))
        statmask=stat;
       cfg.parameter = 'stat';
       cfg.interpmethod  = 'linear';
       statint=ft_sourceinterpolate(cfg, stat,mri);
       statmask.stat=(stat.stat<=tinv(0.01,nsub-1) |stat.stat>=tinv(0.99,nsub-1));
       statmask=ft_sourceinterpolate(cfg, statmask,mri);
 
       statint.mask=statmask.stat~=0;
for h=1:numel(hemisphere)

if sum(statint.mask)==0
    figure
else

    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'stat';
    %cfg.maskparameter = cfg.funparameter;
    cfg.maskparameter  = 'mask';
    cfg.funcolorlim    = [-5 5];
    %cfg.funcolorlim='maxabs';
    cfg.funcolormap    =jet_grey2;
    %cfg.projmethod     = 'nearest';
    cfg.projmethod='sphere_avg';
    cfg.sphereradius=3.5;
    cfg.surffile       = fullfile(project_folder,'scripts','additional_functions',strcat('carret',hemisphere{h},'.mat'));
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
print('-f1','-dtiff','-r600',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_lat_clustercorr_crit01_5.tiff')) 

view(squeeze(views(h,4,:))');
print('-f1','-r600','-dtiff',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_med_clustercorr_crit01_5.tiff')) 

view(squeeze(views(h,2,:))');
print('-f1','-r600','-dtiff',strcat(path_out,stat_list{s}(1:end-4),hemisphere{h},'_vent2_clustercorr_crit01_5.tiff')) 
clear c1 c2 

close all
    end

end