%% PLOTTING AND CLUSTERSTATISTICS

%% first step:
% download fieldtrip,data,layout etc
% save everything in a folder (or usb)

project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes
% add path with your fieldtrip toolbox

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
ft_default.spmversion = 'spm12';
ft_defaults
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%
%% 3D CLUSTERSTATISTICS
path_in=fullfile(project_folder,'data','TF','GA');
path_out=fullfile(project_folder,'data','TF','statistics');
mkdir(path_out);
cd(path_in);

contrasts{1}={'GA_z_enc_TBF_f','GA_z_enc_TBR_f'};
contrasts{2}={'GA_z_enc_TBF_r','GA_z_enc_TBR_r'};

subconds={'GA_z_enc_TBR_r','GA_z_enc_TBR_f','GA_z_enc_TBF_r','GA_z_enc_TBF_f'};

for c=1: numel(contrasts)
    
    conds=contrasts{c};
    ananame=strcat(conds{1},'_',conds{2});
    
    % This defines a filter to smooth the data in time and frequency to make
    % the cluster statistics work more reliable
    cfg=[];
    cfg.fwhm_t=0.2;
    cfg.fwhm_f=2;
    
    % load your data
    load(strcat(conds{1},'.mat'));
    data1=smooth_TF_GAgroup(cfg,data);
    load(strcat(conds{2},'.mat'));
    data2=smooth_TF_GAgroup(cfg,data);
    
    % define neighbours: as you know, clusters are also defined spatially, so
    % fieldtrip needs to know the electrode setup
    cfg                 = [];
    cfg.layout          = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
    cfg.feedback        = 'no';
    cfg.method          = 'triangulation';    ighbours möglich
    cfg.neighbours      = ft_prepare_neighbours(cfg, data2);
    neighbours = cfg.neighbours;
    
    
    cfg = [];
    % this part defines the time and frequency range you want to investigate
    cfg.latency          = [2 4]; 
    cfg.frequency        = [2 30];
    
    % open 3d cluster
    cfg.avgoverfreq = 'no';
    cfg.avgovertime = 'no';
    cfg.avgoverchan = 'no';
    
    % this part defines the statistical methods
    % second-level (= cluster correction)
    cfg.tail             = 0;
    cfg.statistic        = 'depsamplesT'; 
    cfg.alpha            = 0.05;
    cfg.correcttail      = 'prob';
    cfg.neighbours       = neighbours;
    cfg.minnbchan        = 2; 
    
    % first level
    cfg.method           = 'montecarlo';
    cfg.correctm         = 'cluster';
    cfg.clustertail      = 0;
    cfg.clusteralpha     = 0.01;
    cfg.clusterstatistic = 'maxsum';
    cfg.numrandomization = 1000;
    cfg.computecritval = 'yes';
    
    % for between-subjects (depsamplesT)
    design = zeros(1,size(data1.powspctrm,1) + size(data2.powspctrm,1));
    design(1,1:size(data1.powspctrm,1)) = 1;
    design(1,(size(data1.powspctrm,1)+1):(size(data1.powspctrm,1) + size(data2.powspctrm,1)))= 2;
    cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
    
    cfg.design = design;
    [stat] = ft_freqstatistics(cfg, data1, data2);
    
    save(fullfile(path_out,strcat('3d_',ananame,'_',num2str(cfg.frequency(1)),'_',num2str(cfg.frequency(2)),'_',num2str(cfg.latency(1)*1000),'_',num2str(cfg.latency(2)*1000),'_clusterp',num2str(cfg.clusteralpha),'.mat')),'stat');
end

%% 3D Plotting: When & where do your clusters occur? 
path_stats=fullfile(project_folder,'data','TF','statistics');
path_out=fullfile(project_folder,'data','TF','statistics','figures');
mkdir(path_out);

dir_stat=dir(path_stats);

files_ind=find(strncmp('3d_GA',{dir_stat(:).name},4));
stat_list={dir_stat(files_ind).name};
%  load a GA file as plotting template
load(fullfile(project_folder,'data','TF','GA','GA_AVG_z_enc_TBF.mat'));
data1=data;
load(fullfile(project_folder,'scripts','additional_functions','jet_grey2.mat'))


for s=1:numel(stat_list)
    load(fullfile(path_stats,stat_list{s}));
    
    % check for sig clusters
    % positive clusters
    
    if isfield(stat, 'posclusters')
        if ~isempty(stat.posclusters)
            pos_ind=find([stat.posclusters(1,:).prob]<=0.05);
            
            for p=1:numel(pos_ind)
                
                f=figure;
                imagesc(stat.time,stat.freq, squeeze(sum((stat.posclusterslabelmat==pos_ind(p)).*stat.stat)))
                title(['Positive Cluster ',num2str(p),'pvalue=',num2str(stat.posclusters(p).prob)]);
                xlabel('time (s)');
                ylabel('frequency (Hz)');
                set(gca,'YDir','normal')
                colormap('hot')
                
                posmat=zeros(size(stat.posclusterslabelmat));
                posmat(stat.posclusterslabelmat==1)=p;
                % the data structure possig contains the time, frequencies, and labels
                % that make up the significant positive clusters
                possig{p}.time=stat.time(sum(sum(posmat,1),2)>0);
                possig{p}.freq=stat.freq(sum(sum(posmat,1),3)>0);
                possig{p}.labels=stat.label(sum(sum(posmat,2),3)>0);
                
                fig_file=fullfile(path_out,strcat(stat_list{s},'poscluster',num2str(p),'.fig'));
                
                savefig(f,fig_file)
                close all
                
                % plot average t TF in cluster
                f=figure;
                mask=squeeze(sum(posmat,1))~=0;
                hold on
                imagesc(stat.time,stat.freq, squeeze(nanmean(stat.stat(sum(sum(posmat,2),3)>0,:,:))),[-4 4]);
                title(['Positive Cluster ',num2str(p),'pvalue=',num2str(stat.posclusters(p).prob)]);
                colormap(jet_grey2)
                xlabel('time (s)');
                ylabel('frequency (Hz)');
                set(gca,'YDir','normal')
                cRange = caxis;
                contour(stat.time,stat.freq,mask,1,'LineColor','k','LineWidth',1)
                caxis(cRange);
                
                fig_file=fullfile(path_out,strcat(stat_list{s},'poscluster_t',num2str(p),'.fig'));
                
                savefig(f,fig_file)
                close all
                % get dummy data structure to fill and this is achieved by taking the
                % data you still have in your workspace from the statistics or by the
                % file loaded in above
                
                data2plot=data1;                
                data2plot.avg=repmat(squeeze(nansum(nansum(((stat.posclusterslabelmat==p).*stat.stat),2),3)),1,2);
                data2plot.dimord='chan_time';
                data2plot.time=[1, 2];
                data2plot=rmfield (data2plot, 'freq');
                % data2plot.freq=nanmean(stat.freq);
                figure
                cfg=[];
                cfg.parameter          = 'avg';
                cfg.layout = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
                cfg.interactive='no';
                cfg.gridscale          =300;
                cfg.style  ='straight';
                cfg.markersize    =6;
                ft_topoplotER(cfg, data2plot)
                title(['Positive Cluster ',num2str(p),'pvalue=',num2str(stat.posclusters(p).prob)]);
                colormap('hot')
                fig_file=fullfille(path_out,strcat(stat_list{s},'topo_poscluster',num2str(p),'.fig'));
                
                savefig(fig_file)
                close all
                
                
                % plot average t-topo
                f1=min(find(sum(sum(posmat,1),3)>0));
                f2=max(find(sum(sum(posmat,1),3)>0));
                t1=min(find(sum(sum(posmat,1),2)>0));
                t2=max(find(sum(sum(posmat,1),2)>0));
                data2plot.avg=repmat(squeeze(nanmean(nanmean(stat.stat(:,f1:f2,t1:t2),2),3)),1,2);
                
                figure
                cfg=[];
                cfg.parameter          = 'avg';
                cfg.layout = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
                cfg.interactive='no';
                cfg.gridscale          =300;
                cfg.zlim=[-4 4];%
                cfg.style  ='straight_imsat';
                cfg.colormap=jet_grey2;
                cfg.markersize    =2;
                cfg.highlightchannel   =possig{p}.labels;
                cfg.highlightsize      = 6;
                cfg.highlightsymbol      = 'o';
                cfg.markersymbol             =  '.';
                cfg.highlight          = 'on';
                ft_topoplotER(cfg, data2plot)
                title(['Positive Cluster ',num2str(p),'pvalue=',num2str(stat.posclusters(p).prob),...
                    ' freq=',num2str(stat.freq(f1)),'to',num2str(stat.freq(f2)),...
                    ' time=',num2str(stat.time(t1)),'to',num2str(stat.time(t2))]);
                
                fig_file=fullfile(path_out,strcat(stat_list{s},'topo_poscluster_t',num2str(p),'.fig'));
                
                savefig(fig_file)
                close all
            end
        end
    end
    % negative clusters
    
    if isfield(stat, 'negclusters')
        if ~isempty(stat.negclusters)
            neg_ind=find([stat.negclusters(1,:).prob]<=0.05);
            
            for n=1:numel(neg_ind)
                
                f= figure
                imagesc(stat.time,stat.freq, -1*squeeze(sum((stat.negclusterslabelmat==neg_ind(n)).*stat.stat)))
                title(['Negative Cluster ',num2str(n),'pvalue=',num2str(stat.negclusters(n).prob)]);
                xlabel('time (s)');
                ylabel('frequency (Hz)');
                set(gca,'YDir','normal')
                colormap('hot')
                
                negmat=zeros(size(stat.negclusterslabelmat));
                negmat(find(stat.negclusterslabelmat==1))=n;
                % the data structure negsig contains the time, frequencies, and labels
                % that make up the significant negative clusters
                negsig{n}.time=stat.time(sum(sum(negmat,1),2)>0);
                negsig{n}.freq=stat.freq(sum(sum(negmat,1),3)>0);
                negsig{n}.labels=stat.label(sum(sum(negmat,2),3)>0);
                fig_file=fullfile(path_out,strcat(stat_list{s},'negcluster',num2str(n),'.fig'));
                
                savefig(f,fig_file)
                close all
                
                f=figure
                mask=squeeze(sum(negmat,1))~=0;
                hold on
                imagesc(stat.time,stat.freq, squeeze(nanmean(stat.stat(sum(sum(negmat,2),3)>0,:,:))),[-4 4]);
                title(['Negative Cluster ',num2str(n),'pvalue=',num2str(stat.negclusters(n).prob)]);
                xlabel('time (s)');
                ylabel('frequency (Hz)');
                set(gca,'YDir','normal')
                colormap(jet_grey2)
                cRange = caxis;
                contour(stat.time,stat.freq,mask,1,'LineColor','k','LineWidth',1)
                caxis(cRange);
                fig_file=fullfile(path_out,strcat(stat_list{s},'negcluster_t',num2str(n),'.fig'));
                
                savefig(f,fig_file)
                close all
                
                data2plot=data1;
                data2plot.avg=repmat(squeeze(nansum(nansum(((stat.negclusterslabelmat==n).*stat.stat),2),3)),1,2)*-1;
                data2plot.dimord='chan_time';
                data2plot.time=[1, 2];
                data2plot=rmfield (data2plot, 'freq');
                figure
                cfg=[];
                cfg.parameter          = 'avg';
                cfg.layout = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
                cfg.interactive='no';
                cfg.gridscale          =300;
                cfg.style  ='straight';
                cfg.zlim='zeromax';
                cfg.markersize    =6;
                ft_topoplotER(cfg, data2plot)
                title(['Negative Cluster ',num2str(n),'pvalue=',num2str(stat.negclusters(n).prob)]);
                colormap('hot')
                fig_file=fullfile(path_out,strcat(stat_list{s},'topo_negcluster',num2str(n),'.fig'));
                
                savefig(fig_file)
                close all
                
                % plot average t-topo
                f1=min(find(sum(sum(negmat,1),3)>0));
                f2=max(find(sum(sum(negmat,1),3)>0));
                t1=min(find(sum(sum(negmat,1),2)>0));
                t2=max(find(sum(sum(negmat,1),2)>0));
                data2plot.avg=repmat(squeeze(nanmean(nanmean(stat.stat(:,f1:f2,t1:t2),2),3)),1,2);
                
                figure
                cfg=[];
                cfg.parameter          = 'avg';
                cfg.layout = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
                cfg.interactive='no';
                cfg.gridscale          =300;
                cfg.zlim=[-4 4];%
                cfg.style  ='straight_imsat';
                cfg.colormap=jet_grey2;
                cfg.markersize    =2;
                cfg.highlightsize      = 6;
                cfg.highlightsymbol      = 'o';
                cfg.markersymbol             =  '.';
                cfg.highlight          = 'on';
                cfg.highlightchannel   =negsig{n}.labels;
                ft_topoplotER(cfg, data2plot)
                title(['Negative Cluster ',num2str(n),'pvalue=',num2str(stat.negclusters(n).prob),...
                    ' freq=',num2str(stat.freq(f1)),'to',num2str(stat.freq(f2)),...
                    ' time=',num2str(stat.time(t1)),'to',num2str(stat.time(t2))]);
                
                fig_file=fullfile(path_out,strcat(stat_list{s},'topo_negcluster_t',num2str(n),'.fig'));
                savefig(fig_file)
                close all
                
            end
        end
    end
end
