project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%% build grandaverages for sourcestats
toi=[2 3];
foi=[8 13];

beam_type='lcmv10';
path_in= fullfile(project_folder,'source',beam_type);
path_out=(fullfile(project_folder,'source','source_statplots',beam_type));
mkdir(path_out)

load (fullfile(project_folder,'source','sourcemodel_10mm.mat'))

mri=ft_read_mri(fullfile(toolbox_folder, 'fieldtrip-20160122','template','anatomy','single_subj_T1_1mm.nii'));

condition={'TBF_f','TBR_f','TBF_r','TBR_r'};

effects={'rehearsal','inhibition'};

all_subs={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};
c_alpha=0.01;
atlas=ft_read_atlas(fullfile(toolbox_folder,'fieldtrip-20160122','template','atlas','aal','ROI_MNI_V4.nii'));
cfg=[];
   cfg.atlas        = atlas;
   cfg.roi   =atlas.tissuelabel(1:90);
   cfg.inputcoord   = 'mni'; 
 mask = ft_volumelookup(cfg,atlas);
%
mask_vec=mask(:);
mask_volume=atlas;
mask_volume=rmfield(mask_volume,'tissue');
mask_volume=rmfield(mask_volume,'tissuelabel');
mask_volume.mask=mask_vec;


for e=1:numel(effects)
    effect=effects{e};
for n=1:numel(all_subs)
close all
        cd (path_in);
      
         load (fullfile(path_in,all_subs{n},'sourceAll.mat'))

switch effect
     case 'inhibition'
         load (fullfile(path_in,all_subs{n},strcat('enc_',condition{1},'beamed_TF')));
        freq1=freq_cond;
        clear freq_cond
         load (fullfile(path_in,all_subs{n},strcat('enc_',condition{2},'beamed_TF')));
        freq2=freq_cond;          
    case 'rehearsal'
         load (fullfile(path_in,all_subs{n},strcat('enc_',condition{3},'beamed_TF')));
        freq1=freq_cond;
        clear freq_cond
         load (fullfile(path_in,all_subs{n},strcat('enc_',condition{4},'beamed_TF')));
        freq2=freq_cond;
    otherwise
end
               
t1=nearest(freq1.time,toi(1));
t2=nearest(freq1.time,toi(2));

f1=nearest(freq1.freq,foi(1));
f2=nearest(freq1.freq,foi(2));

pow1=squeeze(nanmean(nanmean((freq1.powspctrm(:,f1:f2,t1:t2)),2),3));
pow2=squeeze(nanmean(nanmean((freq2.powspctrm(:,f1:f2,t1:t2)),2),3));

sourceAll.pos = template.sourcemodel.pos;
sourceAll.dim = template.sourcemodel.dim;

insidepos=find(sourceAll.inside);

source1=sourceAll;
source1.avg.pow(insidepos)=pow1;

source2=sourceAll;
source2.avg.pow(insidepos)=pow2;

cfg=[];
cfg.parameter = 'mask';
cfg.interpmethod  = 'nearest';
brain_mask=ft_sourceinterpolate(cfg,mask_volume,source1);
inside_vec=brain_mask.mask(:);

source1.inside=inside_vec;
source2.inside=inside_vec;
source1.avg.pow=source1.avg.pow.*inside_vec;
source2.avg.pow=source2.avg.pow.*inside_vec;

sourcega1{n}=source1;
sourcega2{n}=source2;

clear cfg design source1 source2 sourceAll 
end

design(1, :) = repmat(1:length(all_subs), 1, 2);
design(2, :) = [ones(1, length(all_subs)) ones(1, length(all_subs)) * 2];

cfg = [];
cfg.dim         = sourcega1{1}.dim;
cfg.statistic   = 'depsamplesT';
cfg.parameter   = 'avg.pow';
cfg.correctm    = 'cluster';
cfg.method='montecarlo';
cfg.numrandomization = 10000;
cfg.alpha       = 0.05;
cfg.clusteralpha = c_alpha;
cfg.tail        = 0;
cfg.design = design;
cfg.ivar = 2;% indepva1 
cfg.uvar = 1;% units of observation  should not exist in indepsamp test    

% restrict analysis to in brain areas
stat = ft_sourcestatistics(cfg,sourcega1{:},sourcega2{:});

%  sourcewrite
foi_str=strcat(num2str(foi(1)),'_',num2str(foi(2)));
toi_str=strcat(num2str(toi(1)*1000),'_',num2str(toi(2)*1000));

cd (path_out);
save(strcat('stat_',effect,'_clusteralpha',num2str(c_alpha),foi_str,'_hz',toi_str,'.mat'),'stat');
end


%% load and combine source data
beam_type='lcmv';

path_in= fullfile(project_folder,'source',beam_type);

all_subs={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

load (fullfile(project_folder,'source','sourcemodel_10mm.mat'))
mri=ft_read_mri(fullfile(toolbox_folder, 'fieldtrip-20160122','template','anatomy','single_subj_T1_1mm.nii'));

condition={'TBF_f','TBR_f','TBF_r','TBR_r'};

atlas=ft_read_atlas(fullfile(toolbox_folder,'fieldtrip-20160122','template','atlas','aal','ROI_MNI_V4.nii'));
cfg=[];
   cfg.atlas        = atlas;
   cfg.roi   =atlas.tissuelabel(1:90);
   cfg.inputcoord   = 'mni'; 
 mask = ft_volumelookup(cfg,atlas);
%
mask_vec=mask(:);
mask_volume=atlas;
mask_volume=rmfield(mask_volume,'tissue');
mask_volume=rmfield(mask_volume,'tissuelabel');
mask_volume.mask=mask_vec;

cfg=[];
cfg.parameter = 'mask';
cfg.interpmethod  = 'nearest';
brain_mask=ft_sourceinterpolate(cfg,mask_volume,source);
inside_vec=brain_mask.mask(:);

inside_source=sourceAll.inside;
source_ind=find(inside_source);
inside_pow=logical(inside_vec(source_ind));
% prep data

for c=1:numel(condition)
    for n=1:numel(all_subs)
        
        cd (pathin);
        load (fullfile(path_in,all_subs{n},strcat('enc_',condition{c},'beamed_TF')));
        freq=freq_cond;
        
        % select only inside source channels to match to source stats
        freq.label=freq.label(inside_pow);
        freq.powspctrm=freq.powspctrm(inside_pow,:,:);
        pow_all{n}=freq;        
    end
    cfg=[];
    cfg.keepindividual='yes';
    pow_GA{c}=ft_freqgrandaverage(cfg,pow_all{:});   
end

%% do test in stat defined source roi

toi=[2.5 3];
foi=[8 13];
c_alpha=0.01;
effect='inhibition';
combine_cluster='yes'

% load stat with roi def
path_stat=fullfile(project_folder,'source','source_statplots',beam_type);
cd (path_stat)

foi_str=strcat(num2str(foi(1)),'_',num2str(foi(2)));
toi_str=strcat(num2str(toi(1)*1000),'_',num2str(toi(2)*1000));
stat_file=strcat('stat_',effect,'_clusteralpha',num2str(c_alpha),foi_str,'_hz',toi_str,'.mat');
load(stat_file);

load(fullfile(project_folder,'scripts','additional_functions','jet_grey2.mat'))

atlas=ft_read_atlas(fullfile(toolbox_folder,'fieldtrip-20160122','template','atlas','aal','ROI_MNI_V4.nii'));

neg_check=0;
if isfield(stat, 'negclusters')
    if ~isempty(stat.negclusters)
        neg_check=1;
    end
end

pos_check=0;
if isfield(stat, 'posclusters')
    if ~isempty(stat.posclusters)
        pos_check=1;
    end
end

maxt=[];
max_ind=[];
max_pos=[];
sig_check=[];
max_label=[];

if pos_check
    for clus=1:numel(stat.posclusters)
        % find local maxima of clusters
        stat_tmp=stat.stat.*(stat.posclusterslabelmat==clus);
        inside=~isnan(stat.stat);
        stat_tmp=stat_tmp(inside);
        pos_tmp=stat.pos(inside,:);
        
        clustersize(clus)=sum(stat.posclusterslabelmat==clus);
        [maxt(clus),max_ind(clus)]=max(stat_tmp);
        max_pos(clus,:)=pos_tmp(max_ind(clus),:)
        max_label{clus} = atlas_lookup(atlas, max_pos(clus,:).*10,'inputcoord', 'mni','queryrange',5)
        sig_check(clus)=stat.posclusters(clus).prob<=0.06;
    end
    clusterdef.maxt=maxt;
    clusterdef.max_ind=max_ind;
    clusterdef.max_post=max_pos;
    clusterdef.max_label=max_label;
    clusterdef.sig_check=sig_check;
end

if neg_check
    for clus=1:numel(stat.negclusters)
        
        stat_tmp=stat.stat.*(stat.negclusterslabelmat==clus);
        inside=~isnan(stat.stat);
        stat_tmp=stat_tmp(inside);
        pos_tmp=stat.pos(inside,:);
        clustersize(clus)=sum(stat.negclusterslabelmat==clus);
        [maxt(clus),max_ind(clus)]=min(stat_tmp);
        max_pos(clus,:)=pos_tmp(max_ind(clus),:)
        max_label{clus} = atlas_lookup(atlas, max_pos(clus,:).*10,'inputcoord', 'mni','queryrange',5)
        sig_check(clus)=stat.negclusters(clus).prob<=0.05;
    end
    clusterdef.maxt=maxt;
    clusterdef.max_ind=max_ind;
    clusterdef.max_post=max_pos;
    clusterdef.max_label=max_label;
    clusterdef.sig_check=sig_check;
end

save(fullfile(path_stat,strcat('clusterdef_',stat_file,'.mat')),'clusterdef')

switch combine_cluster
    case 'yes'
        num_cluster=1;
    case 'no'
        num_cluster=sum( clusterdef.sig_check)
end

for clust=1:num_cluster
    switch combine_cluster
        case 'yes'
            clust_label='all';
            % plot tf for different cluster
            %%%%% take care: in stat there are also outside brain voxel!
            if neg_check
                tmp=    stat.negclusterslabelmat(stat.inside);
                sel_cluster=tmp<=max(find(clusterdef.sig_check))&tmp>0;
                
            end
            if pos_check
                tmp=    stat.posclusterslabelmat(stat.inside);
                sel_cluster=tmp<=max(find(clusterdef.sig_check))&tmp>0;
            end
            
        case 'no'
            clust_label=strcat('cluster ',num2str(clust));
            if neg_check
                tmp=    stat.negclusterslabelmat(stat.inside);
                sel_cluster=tmp==clust;
                
            end
            if pos_check
                tmp=    stat.posclusterslabelmat(stat.inside);
                sel_cluster=tmp==clust;
            end
    end
    
    
    % plot bargraph of average power in sig cluster
    t1=nearest(pow_GA{1}.time,toi(1));
    t2=nearest(pow_GA{1}.time,toi(2));
    f1=nearest(pow_GA{1}.freq,foi(1));
    f2=nearest(pow_GA{1}.freq,foi(2));
    
    for c=1:numel(condition)
        submeandata(:,c)= squeeze(nanmean(nanmean(nanmean(pow_GA{c}.powspctrm(:,sel_cluster,f1:f2,t1:t2),2),3),4));
    end
    % barplot for each condition
    fig_title=strcat(path_stat,'avgpow_sigcluster_bars',clust_label,stat_file,'.mat');
    f=figure
    subplot(1,2,1)
    hold on
    bar(mean(submeandata))
    scatter(reshape(repmat(1:4,size(submeandata,1),1),1,[]),reshape(submeandata,1,[]),'k')
    ylabel('mean z');
    set(gca,'XTick',[1 2 3 4]);
    set(gca,'XTickLabel',condition);
    set(gca,'FontSize',16);
    title(fig_title)
    
    subplot(1,2,2)
    % calculate contrasts
    
    [h_inhibition,p_inhibition,x,tstat_inhibition]=ttest(submeandata(:,1),submeandata(:,2));
    [h_rehearsal,p_rehearsal,x,tstat_rehearsal]=ttest(submeandata(:,3),submeandata(:,4));
    [h_interaction,p_interaction,x,tstat_interaction]=ttest(submeandata(:,1)-submeandata(:,2),submeandata(:,3)-submeandata(:,4));
    [h_cond,p_cond,x,tstat_cond]=ttest((submeandata(:,1)+submeandata(:,3)).*0.5,(submeandata(:,2)+submeandata(:,4)).*0.5);
    [h_sme,p_sme,x,tstat_sme]=ttest((submeandata(:,3)+submeandata(:,4)).*0.5,(submeandata(:,1)+submeandata(:,2)).*0.5);
    
    axis off
    hold on
    text(0,0.9,strcat('inhibition: p=',num2str(p_inhibition),' ,t(',num2str(tstat_inhibition.df),')=',num2str(tstat_inhibition.tstat)))
    text(0,0.75,strcat('rehearsal: p=',num2str(p_rehearsal),' ,t(',num2str(tstat_rehearsal.df),')=',num2str(tstat_rehearsal.tstat)))
    text(0,0.6,strcat('interaction: p=',num2str(p_interaction),' ,t(',num2str(tstat_interaction.df),')=',num2str(tstat_interaction.tstat)))
    text(0,0.45,strcat('cond: p=',num2str(p_cond),' ,t(',num2str(tstat_cond.df),')=',num2str(tstat_cond.tstat)))
    text(0,0.3,strcat('sme: p=',num2str(p_sme),' ,t(',num2str(tstat_sme.df),')=',num2str(tstat_sme.tstat)))
       
    savefig(f,strcat(fig_title,'.fig'))
    
   
    figure
    for c=1:numel(condition)
        subplot(2,2,c)
        imagesc(pow_GA{c}.time,pow_GA{c}.freq,squeeze(nanmean(nanmean(pow_GA{c}.powspctrm(:,sel_cluster,:,:),1),2)),[-0.1 0.1])
        set(gca,'YDir','normal')
        title(clust_label)
    end
           
    % run cluster stats    
    GA_inhibition=pow_GA{1};
    GA_inhibition.powspctrm=pow_GA{1}.powspctrm-pow_GA{2}.powspctrm;
    GA_rehearsal=pow_GA{1};
    GA_rehearsal.powspctrm=pow_GA{3}.powspctrm-pow_GA{4}.powspctrm;
    
    cfg = [];    
    % this part defines the time and frequency range you want to investigate
    cfg.latency          = [2 4]; % for the directed forgetting analyses during encoding, the party starts around 2 seconds!
    tois={cfg.latency};
    cfg.frequency        = [3 30];
    fois={cfg.frequency};
    
    cfg.avgoverfreq = 'no';
    cfg.avgovertime = 'no';
    cfg.avgoverchan = 'yes';
    cfg.channel=pow_GA{1}.label(sel_cluster);

    cfg.tail             = 0;
    cfg.statistic        = 'depsamplesT'; % within subjects design = depsamples!
    cfg.alpha            = 0.05;
    
    % first level
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = 10000;
    cfg.computecritval = 'yes';
   
            % for within-subjects (depsamplesT)
            Nsub = size(pow_GA{1}.powspctrm,1);                                       %# of subjects?
            design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
            design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
            
            cfg.uvar     = 2;
            cfg.ivar     = 1;
    cfg.design = design;
    
    % rehearsal
    [stat_source] = ft_freqstatistics(cfg, pow_GA{3},pow_GA{4});
    h=figure
    mask=squeeze(stat_source.mask);
    hold on
    imagesc(stat_source.time,stat_source.freq, squeeze(stat_source.stat),[-4 4]);
    colormap(jet_grey2)
    xlabel('time (s)');
    ylabel('frequency (Hz)');
    set(gca,'YDir','normal')
    cRange = caxis;
    contour(stat_source.time,stat_source.freq,mask,1,'LineColor','k','LineWidth',1)
    caxis(cRange);
    
    [~,~,~,~,sig_text]=mf_clustercheck(stat_source,0.05)
    title(strcat('rehearsal',clust_label,sig_text))
    savefig(h,fullfile(path_stat,'figures',strcat('tf_plot_rehearsal_',clust_label,stat_file,'.fig')))
    
    [stat_source] = ft_freqstatistics(cfg, pow_GA{1},pow_GA{2})
    g=figure
    mask=squeeze(stat_source.mask);
    hold on
    imagesc(stat_source.time,stat_source.freq, squeeze(stat_source.stat),[-4 4]);
    colormap(jet_grey2)
    xlabel('time (s)');
    ylabel('frequency (Hz)');
    set(gca,'YDir','normal')
    cRange = caxis;
    contour(stat_source.time,stat_source.freq,mask,1,'LineColor','k','LineWidth',1)
    caxis(cRange);
    [~,~,~,~,sig_text]=mf_clustercheck(stat_source,0.05)
    title(strcat('inhibition',clust_label,sig_text))
    savefig(g,fullfile(path_stat,'figures',strcat('tf_plot_inhibition_',clust_label,stat_file,'.fig')))
    
end

 %%
 % save average power in cluster for correlations
 t1=nearest(pow_GA{1}.time,toi(1));
 t2=nearest(pow_GA{1}.time,toi(2));
 f1=nearest(pow_GA{1}.freq,foi(1));
 f2=nearest(pow_GA{1}.freq,foi(2));
 for i=1:4
 avg_alpha(:,i)=nanmean(nanmean(nanmean(pow_GA{i}.powspctrm(:,sel_cluster,f1:f2,t1:t2),2),3),4);
 end
 alpha_cond=condition;
save(fullfile(path_stat,strcat('extractedavg_alphaavg_pow_source_',stat_file,'.mat')),'avg_alpha','alpha_cond')

 %% correlations inhibition and rehearsal
cd(path_stat)
load('extractedmax_alphaavg_pow_source_stat_inhibition_clusteralpha0.018_13_hz2000_2500.mat.mat')
inhibition_alpha1=avg_alpha;
load('extractedmax_alphaavg_pow_source_stat_inhibition_clusteralpha0.018_13_hz2500_3000.mat.mat')
inhibition_alpha2=avg_alpha;

%load avg rehearsal alpha in rehearsal cluster
load('extractedmax_alphaavg_pow_source_stat_rehearsal_clusteralpha0.018_13_hz3000_3500.mat.mat')
rehearsal_alpha1=avg_alpha;
load('extractedmax_alphaavg_pow_source_stat_rehearsal_clusteralpha0.018_13_hz3500_4000.mat.mat')
rehearsal_alpha2=avg_alpha;

% check for correlation inhibition/rehearsal effect: power
[rho,p]=corr(alpha1_inhibition,alpha1_rehearsal,'Type','Spearman')
[rho,p]=corr(alpha2_inhibition,alpha2_rehearsal,'Type','Spearman')