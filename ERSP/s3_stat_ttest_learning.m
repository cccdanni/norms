% ERSP analysis
% step 3: statistical analysis
% Project: Social Norms Learning
% Dependency: FieldTrip
% Author: Danni Chen 
% Update Date: 2022-03-02

clear all; clc; close all;

%% basic information setup

% Subject info
subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238, 1210, 1211, 1223, 1239]; 
subNameD = setdiff(subNameD, badsubNameD);

% folder info
workFolder          = '/home/chendanni/Documents/Norms/analysis/';
ERSPParentFolder    = fullfile (workFolder, 'EEGAnalysisResults', 'ERSP');
ScriptsFolder       = fullfile(workFolder, 'MyScripts');
FieldTripFolder     = '/home/chendanni/MATLAB/toolbox/fieldtrip-20210330';
addpath ( genpath(ScriptsFolder) );
addpath ( genpath(FieldTripFolder) );

% parameter setup
bltw = [-0.300 0.000]; % in secs
toiw = [-0.300 3.000]; % in secs
bindescr = {'ingroup-higher', 'ingroup-lower', 'ingroup-consistent',...
        'outgroup-higher', 'outgroup-lower', 'outgroup-consistent'};
task = {'Learning'};

% ROI info
ROI = struct;
ROIs = {'FC','Ce','CP','LP','RP', 'LOT','ROT','OT'}; 
for i = 1:8 ROI(i).ROIName = ROIs{i}; end
ROI(1).ChanList = {'Fz', 'FCz', 'F1', 'F2', 'FC1', 'FC2'}; %FC
ROI(2).ChanList = {'Cz', 'C1', 'C2'}; %Ce
ROI(3).ChanList = {'CP1', 'CP2', 'Pz', 'P1', 'P2'}; %CP
ROI(4).ChanList = {'P3', 'P5', 'CP3', 'CP5'}; %LP
ROI(5).ChanList = {'P4', 'P6', 'CP4', 'CP6'}; %RP
ROI(6).ChanList = {'T7','TP7','P7','PO7'}; %LOT
ROI(7).ChanList = {'T8','TP8','P8','PO8'}; %ROT
ROI(8).ChanList = {'T7','TP7','P7','PO7','T8','TP8','P8','PO8'}; %OT

% two conditions ttest
ConditionName = {'ingroup-higher-vs-lower', 'ingroup-higher-vs-consistent', 'ingroup-consistent-vs-lower',...
    'outgroup-higher-vs-lower', 'outgroup-higher-vs-consistent', 'outgroup-consistent-vs-lower',...
    'higher-ingroup-vs-outgroup', 'lower-ingroup-vs-outgroup','consistent-ingroup-vs-outgroup'};
ConditionBin = [1,2; 1,3; 3,2; 4,5; 4,6; 6,5; 1,4; 2, 5; 3,6];

%% statistical analysis

for iTask = 1:length(task)

    % locate dataset
    curTask = task {iTask};
    outputERSPFolder = fullfile(ERSPParentFolder,'Stats', strcat(curTask, 'Stats'), 'depsampleT');
    if ~exist(outputERSPFolder)
        mkdir(outputERSPFolder);
    end
    inputERSPFolder = fullfile(ERSPParentFolder,'Results', strcat(curTask, 'Avg'));
    
    %% dependent sample t-test ----------
    
    for ibin = 1:length(ConditionName) % loop for each condition
        
        % load data
        avg_all = [];
        curbin = ConditionName{ibin};
        curbincond = ConditionBin(ibin, :);
        curbincond_name = {bindescr{curbincond(1)}, bindescr{curbincond(2)}}
        
        data = {};
        
        for i_cond = 1:length(curbincond)
        
            tmpname = curbincond_name{i_cond};
        
            if ~strcmp(curTask, 'postMinuspre')
                load (fullfile (inputERSPFolder, strcat(curTask,'-',tmpname,'-grandavg.mat')));
            else
                pre = load (fullfile (fullfile(ERSPParentFolder, 'Results',strcat(task{1}, 'Avg')),  strcat(task{1},'-',tmpname,'-grandavg.mat')));
                post = load (fullfile (fullfile(ERSPParentFolder, 'Results', strcat(task{2}, 'Avg')),  strcat(task{2},'-',tmpname,'-grandavg.mat')));
                avg_all = pre.avg_all;
                for k = 1:length(avg_all)
                    avg_all{k}.powspctrm = post.avg_all{k}.powspctrm - pre.avg_all{k}.powspctrm;
                end
            end
            
            data{i_cond} = avg_all; % include cond 1 & cond 2
            clear avg_all;
        
        end
        
        
        % statistical analysis
        for iROI = 1:length(ROIs) % loop for each ROI
            
            curROI = ROI(iROI).ROIName;
            COI = ROI(iROI).ChanList;
            
            avg1 = data{1};
            avg2 = data{2};


            cfg               = [];
            cfg.channel       = COI;
            cfg.latency       = [-0.300 2.996];
            cfg.frequency     = 'all';
            cfg.avgoverchan   = 'yes';
            cfg.avgovertime   = 'no';
            cfg.avgoverfreq   = 'no';
            cfg.parameter     = 'powspctrm';
            cfg.correctm      = 'cluster';
            cfg.clusteralpha  = 0.05;
            cfg.clusterstatistic = 'maxsum';
            cfg.method        = 'montecarlo';
            cfg.statistic     = 'ft_statfun_depsamplesT'; 
            cfg.numrandomization = 1000;
            cfg.computestat   = 'yes';
            cfg.computecritical = 'yes';
            cfg.computeprop   = 'yes';
            cfg.tail = 1;

            subj = length(avg1);
            design = zeros (2, subj);
            for i = 1:subj
                design(1,i) = i;
            end
            for i = 1:subj
                design(1,subj+i) = i;
            end
            design(2, 1:subj) = 1;
            design(2, (subj+1):end) = 2;
            design = design';

            cfg.design        = design;
            cfg.uvar          = 1; % number or list with indices, unit variable(s), corresponds to the subject number (in a within-subject manipulation)
            cfg.ivar          = 2; % number or list with indices, independent variable(s)

            stat = ft_freqstatistics(cfg, avg1{:}, avg2{:});

            cfg = [];
            cfg.channel = COI;
            cfg.latency = [-0.300 2.996];
            freq1_grandavg = ft_freqgrandaverage(cfg, avg1{:});
            freq2_grandavg = ft_freqgrandaverage(cfg, avg2{:});
            freq1_grandavg = ft_freqdescriptives(cfg, freq1_grandavg);
            freq2_grandavg = ft_freqdescriptives(cfg, freq2_grandavg);

            % plotting the results
            stat.raweffect = freq1_grandavg.powspctrm - freq2_grandavg.powspctrm;
            stat.raweffect = mean(stat.raweffect, 1);
            
            cfg = [];
            ztickformat('%.2f');
            cfg.channel       = 'all';
            cfg.avgoverchan   = 'yes';
            cfg.zlim          = [-0.10 0.10];
            cfg.parameter     = 'raweffect';
            cfg.maskparameter = 'mask';
            %cfg.maskstyle    = 'outline';
            cfg.maskalpha     = 0.2;
            figure;
            cfg.colormap      = '*RdBu';
            ft_singleplotTFR(cfg,stat);
            
            
            % save figure and stats
            cd(outputERSPFolder);
            
            saveas(gcf,strcat( curTask,'-', curbin, '-', curROI),'png');
            close all;
            save (strcat( curTask,'-', curbin, '-', curROI, '-depsampleT.mat'), 'stat', 'curROI', 'COI', 'curbin')
            
        end
        
        
    end

end