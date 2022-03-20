
project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%% erp analysis for supplements

path_in=fullfile(project_folder,'data');
path_out=fullfile(project_folder,'data','erp')
mkdir(path_out)
all_subs={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

condition={'TBF-f','TBR-f','TBF-r','TBR-r'};
%
for n=1:numel(all_subs)
    sel_sub=all_subs{n};
    load(fullfile(path_in,strcat(sel_sub,'_data')));
    cfg=[];
    cfg.lpfilter='yes';
    cfg.lpfreq=10;
    data=ft_preprocessing(cfg,data);
    for c=1:numel(condition)
        sel_condition=condition{c}
        switch sel_condition
            case 'TBF-f'
                trial_def=find(data.trialinfo(:,5)==13 & data.trialinfo(:,10)==0);
            case 'TBF-r'
                trial_def=find(data.trialinfo(:,5)==13 & data.trialinfo(:,10)==1);
            case 'TBR-f'
                trial_def=find(data.trialinfo(:,5)==11 & data.trialinfo(:,10)==0);
            case 'TBR-r'
                trial_def=find(data.trialinfo(:,5)==11 & data.trialinfo(:,10)==1);
        end
        cfg=[];
        cfg.trials=trial_def;
        cfg.keeptrials='no';
        erp=ft_timelockanalysis(cfg,data)
        
        cfg=[];
        cfg.baseline=[1.8 2];
        erp=ft_timelockbaseline(cfg,erp);
        all_erp{n,c}=erp;
        clear erp
    end
end

for c=1:numel(condition)
    cfg=[];
    cfg.keepindividual='yes';
    ga{c}=ft_timelockgrandaverage(cfg,all_erp{:,c})
end
clear all_erp

%% do stats

% define contrasts

sel_contrast='rehearsal';

switch sel_contrast
    case 'interaction'
        data1=ga{1};
        data2=ga{3};
        data1.individual=ga{1}.individual-ga{2}.individual;
        data2.individual=ga{3}.individual-ga{4}.individual;
        
    case 'rehearsal'
        data1=ga{3};
        data2=ga{4};
    case 'inhibition'
        data1=ga{1};
        data2=ga{2};
    otherwise
        error('no contrast defined')
end


% define neighbours
cfg                 = [];
cfg.layout          = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
cfg.method          = 'triangulation';      
cfg.neighbours      = ft_prepare_neighbours(cfg, data2);
ft_neighbourplot(cfg, data1) 
neighbours = cfg.neighbours;
%%
% Now the actual statistical calculation is defined
load(fullfile(project_folder,'scripts','additional_functions','jet_grey_halfmax.mat'))

cfg = [];
% this part defines the time  range you want to investigate
cfg.latency          = [2.5 3]; %
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'no';

% this part defines the statistical methods
% first-level (= t-tests)
cfg.tail             = 0;
cfg.statistic        = 'depsamplesT'; % within subjects design = depsamples!
cfg.alpha            = 0.05;

cfg.neighbours       = neighbours;
cfg.minnbchan        = 2; % how many neighbouring channels need to exceed the cfg.alpha together in order to make a cluster

% second level (= cluster correction)
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clustertail      = 0;
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.numrandomization = 1000;
cfg.computecritval = 'yes';

% for within-subjects (depsamplesT)
Nsub = size(data1.individual,1);                                       %# of subjects?
design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

cfg.uvar     = 2;
cfg.ivar     = 1;

cfg.design = design;

[stat] = ft_timelockstatistics(cfg, data1,data2);


data2plot=data1;
data2plot.individual=nanmean(stat.stat,2);
data2plot.time=nanmean(stat.time);
data2plot.dimord='chan_time'
f=figure
cfg=[];
cfg.parameter          = 'individual';
cfg.layout          = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
cfg.interactive='no';
cfg.gridscale          =300;
cfg.style  ='straight_imsat';
cfg.colormap=jet_grey2;
cfg.highlightchannel   =(find(sum(stat.mask,2)));
cfg.highlightsize      = 6;
cfg.highlightsymbol      = 'o';
cfg.highlight='on';
cfg.markersymbol             =  '.';
cfg.zlim=[-5 5];%
cfg.comment            = 'no';
ft_topoplotER(cfg, data2plot)
colorbar

savefig(f, fullfile(path_out,strcat(sel_contrast,'_topo_time',num2str(tois(1).*1000),'to',num2str(tois(2).*1000))),'compact')
close all


% erps for each condition
plot_toi=[-0.2 4];
t1=nearest(data1.time, plot_toi(1));
t2=nearest(data1.time, plot_toi(2));
% t-values for inhibition and rehearsal
inhibition=ga{1}.individual-ga{2}.individual;
rehearsal=ga{3}.individual-ga{4}.individual;
[p,h,~,stats]=ttest(inhibition);
t_inhibition=squeeze(stats.tstat);
[p,h,~,stats]=ttest(rehearsal);
t_rehearsal=squeeze(stats.tstat);

if isfield(stat, 'posclusters')
    if ~isempty(stat.posclusters)
        if any([stat.posclusters(:).prob]<=0.05)
            sig_elecs=find(sum(stat.posclusterslabelmat<=sum([stat.posclusters(:).prob]<=0.05)&...
                stat.posclusterslabelmat~=0,2));
            
            %plot erps for all conditions for these electrodese
            for c=1:numel(condition)
                erp(c,:)=squeeze(mean(mean(ga{c}.individual(:,sig_elecs,t1:t2))));
            end
            
            f=figure
            subplot(1,2,1)
            plot(ga{1}.time(t1:t2),erp')
            legend(condition)
            xlabel(strcat('poscluster, p=',num2str(stat.posclusters(1).prob)))
            
            subplot(1,2,2)
            plot(ga{1}.time(t1:t2),squeeze(mean(mean(inhibition(:,sig_elecs,t1:t2)))),'r')
            hold on
            plot(ga{1}.time(t1:t2),squeeze(mean(mean(rehearsal(:,sig_elecs,t1:t2)))),'g')
            legend('inhibition','rehearsal')
            savefig(f, fullfile(path_out,strcat(sel_contrast,'_poscluster_time',num2str(tois(1).*1000),'to',num2str(tois(2).*1000))),'compact')
        end
    end
end

if isfield(stat, 'negclusters')
    if ~isempty(stat.negclusters)
        if any([stat.negclusters(:).prob]<=0.05)
            sig_elecs=find(sum(stat.negclusterslabelmat<=sum([stat.negclusters(:).prob]<=0.05)&...
                stat.negclusterslabelmat~=0,2));
            
            %plot erps for all conditions for these electrodese
            for c=1:numel(condition)
                erp(c,:)=squeeze(mean(mean(ga{c}.individual(:,sig_elecs,t1:t2))));
            end
            f=figure
            subplot(1,2,1)
            plot(ga{1}.time(t1:t2),erp')
            legend(condition)
            xlabel(strcat('negcluster, p=',num2str(stat.negclusters(1).prob)))
            subplot(1,2,2)
            plot(ga{1}.time(t1:t2),squeeze(mean(mean(inhibition(:,sig_elecs,t1:t2)))),'r')
            hold on
            plot(ga{1}.time(t1:t2),squeeze(mean(mean(rehearsal(:,sig_elecs,t1:t2)))),'g')
            legend('inhibition','rehearsal')
            savefig(f, fullfile(path_out,strcat(sel_contrast,'_negcluster_time',num2str(tois(1).*1000),'to',num2str(tois(2).*1000))),'compact')
        end
    end
end
