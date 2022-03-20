project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes

addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));

%% trial based permutation for alpha power effects
path_in=fullfile(project_folder,'data');
path_out=fullfile(project_folder,'data','TF','single_trial_alpha');
foi=[8 13];

mkdir(path_out);
cd(path_in);

all_subs ={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

for n=1:numel(all_subs)
    sel_sub=all_subs{n};
    load(fullfile(path_data,strcat(sel_sub,'_data')));
    
    % this piece of script is actually there to correct that we had too
    % short trial lengths
    for i=1:length(data.time);
        data.trial{1,i} = horzcat(fliplr(data.trial{1,i}(:,2:251)),data.trial{1,i}(:,:),fliplr(data.trial{1,i}(:,1250:1499)));
        data.time{1,i}=-2:0.0040:5.9960;
    end
    % run tf
    
    cfg=[];
    cfg.method='wavelet';
    cfg.foi=foi(1):0.5:foi(2);
    cfg.toi=-1:0.05:5;
    cfg.width=5;
    cfg.output='pow';
    cfg.keeptrials='yes';
    cfg.trials=enco_trials;
    data=ft_freqanalysis(cfg, data);
    
    % ztrans across all trials
    cfg=[];
    cfg.time=[-0.5 4.5];
    data=z_trans_TF_seltime(cfg,data);
    
    % average for alpha 
    data.powspctrm=nanmean(data.powspctrm,3);
    data.freq=[10];
    data.cumtapcnt=data.cumtapcnt(:,1)
    all_data{n}=data;
end

% get neighbours
cfg                 = [];
cfg.layout          = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
cfg.feedback        = 'no';
cfg.method          = 'triangulation';      % 3 Methoden zur Bestimmung der Neighbours möglich
neighbours      = ft_prepare_neighbours(cfg, data);

% define test settings (time and condition to contrast)

nrand=1000;
all_contrasts={'inhibition','rehearsal'};

all_toi=[2,2.5;2.5,3;3,3.5;3.5,4]
for t=1:size(all_toi,1)
    toi=all_toi(t,:);
    for c=1:numel(all_contrasts)
        contrast=all_contrasts{c};
        
        
        for n=1:numel(selvps)
            % select trials based on condition
            switch contrast
                case 'inhibition'
                    sel_trials1{n}=find(all_data{n}.trialinfo(:,5)==13 &all_data{n}.trialinfo(:,10)==0);
                    sel_trials2{n}=find(all_data{n}.trialinfo(:,5)==11 &all_data{n}.trialinfo(:,10)==0);
                case 'rehearsal'
                    sel_trials1{n}=find(all_data{n}.trialinfo(:,5)==13 &all_data{n}.trialinfo(:,10)==1);
                    sel_trials2{n}=find(all_data{n}.trialinfo(:,5)==11 &all_data{n}.trialinfo(:,10)==1);
            end
        end
        
        data_tmp=rmfield(data,'powspctrm')
        data_tmp.powspctrm=zeros(numel(selvps),numel(data.label),1,numel(data.time));
        data_tmp.dimord='subj_chan_freq_time';
        for r=1:nrand
            data1=data_tmp;
            data2=data_tmp;
            for n=1:numel(selvps)
                % generate random trials
                num1=numel(sel_trials1{n});
                num2=numel(sel_trials2{n});
                tmp_trials=[sel_trials1{n};sel_trials2{n}];
                rand_ind=randperm(num1+num2);
                tmp_trials1=tmp_trials(rand_ind(1:num1));
                tmp_trials2=tmp_trials(rand_ind(num1+1:end));
                
                data1.powspctrm(n,:,:,:)=nanmean(all_data{n}.powspctrm(tmp_trials1,:,:,:));
                data2.powspctrm(n,:,:,:)=nanmean(all_data{n}.powspctrm(tmp_trials2,:,:,:));
                
            end
            % run freqstats
            cfg = [];           
            cfg.latency          = toi;
            cfg.frequency        = [10 10];
            cfg.avgoverfreq = 'yes';
            cfg.avgovertime = 'yes';
            cfg.avgoverchan = 'no';
            cfg.tail             = 0;
            cfg.statistic        = 'depsamplesT'; % within subjects design = depsamples!
            cfg.alpha            = 0.05;
            cfg.neighbours       = neighbours;
            cfg.minnbchan        = 2; % how many neighbouring channels need to exceed the cfg.alpha together in order to make a cluster
            cfg.method           = 'montecarlo';
            cfg.correctm         = 'cluster';
            cfg.clustertail      = 0;
            cfg.clusteralpha     = 0.05;
            cfg.clusterstatistic = 'maxsum';
            cfg.numrandomization = 1;
            cfg.computecritval = 'yes';
            
            Nsub = size(data1.powspctrm,1);                                       %# of subjects?
            design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
            design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
            
            cfg.uvar     = 2;
            cfg.ivar     = 1;
            cfg.design = design;
            
            stat = ft_freqstatistics(cfg, data1, data2);
            
            if isfield(stat,'posclusters')
                if ~isempty(stat.posclusters)
                    poscluster(r)=stat.posclusters(1).clusterstat;
                else
                    poscluster(r)=0;
                end
            else
                poscluster(r)=0;
            end
            
            if isfield(stat,'negclusters')
                if ~isempty(stat.negclusters)
                    negcluster(r)=stat.negclusters(1).clusterstat;
                else
                    negcluster(r)=0;
                end
            else
                negcluster(r)=0;
            end
            
        end
        
        poscluster=sort(poscluster,'descend');
        negcluster=sort(negcluster,'ascend');
        
        % run freqstat on real data
        for n=1:numel(selvps)
            data1.powspctrm(n,:,:,:)=nanmean(all_data{n}.powspctrm(sel_trials1{n},:,:,:));
            data2.powspctrm(n,:,:,:)=nanmean(all_data{n}.powspctrm(sel_trials2{n},:,:,:));
        end
        
        
        cfg = [];
        cfg.latency          = toi;
        cfg.frequency        = [10 10];
        cfg.avgoverfreq = 'yes';
        cfg.avgovertime = 'yes';
        cfg.avgoverchan = 'no';
        cfg.tail             = 0;
        cfg.statistic        = 'depsamplesT'; % within subjects design = depsamples!
        cfg.alpha            = 0.05;
        cfg.neighbours       = neighbours;
        cfg.minnbchan        = 2; % how many neighbouring channels need to exceed the cfg.alpha together in order to make a cluster
        cfg.method           = 'montecarlo';
        cfg.correctm         = 'cluster';
        cfg.clustertail      = 0;
        cfg.clusteralpha     = 0.05;
        cfg.clusterstatistic = 'maxsum';
        cfg.numrandomization = 1;
        cfg.computecritval = 'yes';
        
        Nsub = size(data1.powspctrm,1);                                       %# of subjects?
        design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
        design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
        
        cfg.uvar     = 2;
        cfg.ivar     = 1;
        cfg.design = design;
        
        stat = ft_freqstatistics(cfg, data1, data2);
        
        
        % get pvalue using random distribution
        if isfield(stat,'posclusters')
            if ~isempty(stat.posclusters)
                p_pos=nearest(poscluster, stat.posclusters(1).clusterstat)/nrand
            else
                p_pos=1
            end
        else
            p_pos=1
        end
        
        if isfield(stat,'negclusters')
            if ~isempty(stat.negclusters)
                p_neg=nearest(negcluster, stat.negclusters(1).clusterstat)/nrand
            else
                p_neg=1
            end
        else
            p_neg=1
        end
        load(fullfile(project_folder,'scripts','additional_functions','jet_grey2.mat'))
        data2plot=data1;
        data2plot.avg=stat.stat;
        data2plot.dimord='chan_time';
        data2plot.time=[1];
        data2plot.freq=nanmean(stat.freq);
        f=figure
        cfg=[];
        cfg.parameter          = 'avg';
        cfg.layout          = fullfile(project_folder,'scripts','additional_functions','BrainCap64_1020_lay.mat');
        cfg.interactive='no';
        cfg.gridscale          =300;
        cfg.zlim=[-4 4];
        cfg.colorbar ='yes';
        cfg.comment='no';
        cfg.style  ='straight_imsat';
        cfg.colormap=jet_grey2;
        cfg.markersize    =2;
        cfg.highlightsize      = 6;
        cfg.highlightsymbol      = 'o';
        cfg.markersymbol             =  '.';
        cfg.highlight          = 'on';
        if p_pos<0.05
            cfg.highlightchannel=stat.label(stat.posclusterslabelmat==1);
        elseif p_neg<0.05
            cfg.highlightchannel=stat.label(stat.negclusterslabelmat==1);
        end
        cfg.colorbartext       =strcat('p pos=',num2str(p_pos),', p neg=',num2str(p_neg));
        ft_topoplotER(cfg, data2plot)
        savefig(f,fullfile(path_out,strcat('fig_',contrast,'toi',num2str(toi(1).*1000),'to',num2str(toi(1).*1000))));
        close all        
    end
end
