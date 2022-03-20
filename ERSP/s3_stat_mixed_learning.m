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
toiw = [-0.200 3.000]; % in secs
bindescr = {'ingroup-higher', 'ingroup-lower', 'ingroup-consistent',...
        'outgroup-higher', 'outgroup-lower', 'outgroup-consistent'};
task          = {'Learning'};

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

% F-test conditions
ConditionName = {'mixed'};
ConditionBin = [1:6];

%% statistical analysis

for iTask = 1:length(task)

    % locate dataset
    curTask = task {iTask};
    outputERSPFolder = fullfile(ERSPParentFolder,'Stats', strcat(curTask, 'Stats'), 'onesampleT');
    if ~exist(outputERSPFolder)
        mkdir(outputERSPFolder);
    end
    inputERSPFolder = fullfile(ERSPParentFolder,'Results', strcat(curTask, 'Avg'));
    
    %% f-test 
    
    for icond = 1:length(ConditionName)
        outputERSPFolder = fullfile(ERSPParentFolder,'Stats', strcat(curTask, 'Stats'), 'ftest');
        if ~exist (outputERSPFolder)
            mkdir (outputERSPFolder);
        end
        
        curcondbin = ConditionBin(icond, :);
        curcondname = ConditionName {icond};
        grandavg_all = [];      
        
        % load data
        avg_all = [];
        if ~strcmp(curTask, 'postMinuspre')
            for ibin = 1:length(curcondbin)
                avg_all = [];
                curbin = bindescr{curcondbin(ibin)};
                load (fullfile (inputERSPFolder, strcat(curTask,'-',curbin,'-grandavg.mat')));
                grandavg_all{ibin} = avg_all;
            end
        else
            for ibin = 1:length(curcondbin)
                avg_all = [];
                curbin = bindescr{curcondbin(ibin)};
                pre = load (fullfile (fullfile(ERSPParentFolder, 'Results',strcat(task{1}, 'Avg')),  strcat(task{1},'-',curbin,'-grandavg.mat')));
                post = load (fullfile (fullfile(ERSPParentFolder, 'Results',strcat(task{2}, 'Avg')),  strcat(task{2},'-',curbin,'-grandavg.mat')));
                avg_all = pre.avg_all;
                for k = 1:length(avg_all)
                    avg_all{k}.powspctrm = post.avg_all{k}.powspctrm - pre.avg_all{k}.powspctrm;
                end
                grandavg_all{ibin} = avg_all;
            end            
        end
        
        for iROI = 1:length(ROIs) % loop for each ROI
            
            curROI = ROI(iROI).ROIName;
            COI = ROI(iROI).ChanList;
            
            avg1 = grandavg_all{1};
            avg2 = grandavg_all{2};
            avg3 = grandavg_all{3};
            avg4 = grandavg_all{4};
            avg5 = grandavg_all{5};
            avg6 = grandavg_all{6};

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
            cfg.statistic     = 'depsamplesFmultivariate'; 
            cfg.numrandomization = 1000;
            cfg.computestat   = 'yes';
            cfg.computecritical = 'no';
            cfg.computeprop   = 'yes';
            cfg.tail          = 1; % this has to be 1 for depsamplesF
            cfg.clustertail   = 1; % this has to be 1 for depsamplesF

            subj = length(avg);
            design = zeros (2, subj*3);
            for i = 1:subj
                design(1,i) = i;
            end
            for i = 1:subj
                design(1,subj+i) = i;
            end
            for i = 1:subj
                design(1,subj*2+i) = i;
            end
            design(2, 1:subj) = 1;
            design(2, (subj+1):(subj*2)) = 2;
            design(2, (subj*2+1):end) = 3;
            % design = design';

            cfg.design        = design;
            cfg.uvar          = 1; % number or list with indices, unit variable(s), corresponds to the subject number (in a within-subject manipulation)
            cfg.ivar          = 2; % number or list with indices, independent variable(s)

            stat = ft_freqstatistics(cfg, avg1{:}, avg2{:}, avg3{:});

            % plotting the results
            cfg = [];
            cfg.channel   = 'all';
            cfg.avgoverchan = 'yes';
            cfg.parameter = 'stat';
            cfg.maskparameter = 'mask';
            %cfg.maskstyle = 'outline';
            cfg.maskalpha  = 0.2;
            figure;
            cfg.colormap = '*RdBu';
            ft_singleplotTFR(cfg,stat);
            
            % save figure and stats
            cd(outputERSPFolder);
            saveas(gcf,strcat( curTask,'-', curcondname, '-', curROI),'png');
            close all;
            save (strcat( curTask,'-', curcondname, '-', curROI, '-ftest.mat'), 'stat', 'curROI', 'COI', 'curcondname')
        
        end
    end

end