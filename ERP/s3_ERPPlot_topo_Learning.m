% Plot ERP topography with FieldTrip
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: 8/30/2021

%% basic setting

% Subject info
clear all; clc; close all;
subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238, 1210, 1211, 1223, 1239]; 
% 1201 - 1204: HK Local; 
% 1214, 1215 & 1238 - visually detected bad subjects
% 1210, 1211, 1223, 1239 - removed more than 100 trials
subNameD = subNameD(~ismember(subNameD, badsubNameD));

% folder info
workFolder = '/home/chendanni/Documents/Norms/analysis/';
task = {'Learning'};
ERPParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'ERP'); % set your output path
PreprocessParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing');
ScriptsFolder = '/home/chendanni/Documents/Norms/analysis/MyScripts/ERP';
FieldTripFolder = '/home/chendanni/MATLAB/toolbox/fieldtrip-20210330';
addpath ( genpath(ScriptsFolder) );
addpath ( genpath(FieldTripFolder) );
addpath ( genpath("/home/chendanni/Documents/Norms/analysis/MyScripts/additional_functions/boundedline-pkg-master") );

bl = [-300 0];

%% Plot ERP topo

for iTask = 1:length(task)
    
    curTask = task{iTask};
    
    % component of interest info
    switch curTask
        case 'Learning'
            bltw = [-0.30 0];
            tw   = [-300 2996];
            demoName = 'demo_learning.mat';
            complist = {'lpc','n400','frn','lpp', 'p300'};
            twlist = [0.55 0.55;...
              0.40 0.40;...
              0.35 0.35;...
              1.30 1.30;...
              0.30 0.30];
        case 'preLearningImplicit'
            bltw = [-0.20 0];
            tw   = [-200 996];
            complist = {'n170','lpc','epn'};
            demoName = 'demo_implicit.mat';
            twlist = [0.17 0.17;...
              0.55 0.55;...
              0.25 0.25];
        case 'postLearningImplicit'
            bltw = [-0.20 0];
            tw   = [-200 996];
            demoName = 'demo_implicit.mat';
            complist = {'n170','lpc','epn'};
            twlist = [0.17 0.17;...
              0.55 0.55;...
              0.25 0.25];
    end
    
    % load demo 
    cd(ScriptsFolder);
    data = [];
    load (demoName);
    
    % load ERP data
    outputERPFolder = fullfile (ERPParentFolder, strcat(curTask, 'ERPlabGrandAvg'));
    cd (outputERPFolder);    
    ERP = pop_loaderp( 'filename', [curTask, '_GrandAvg_ManuRej', num2str(bl(1)), '.erp'], 'filepath',pwd);
    
    % topography
    for iBin = 1:6 % loop each bin 
        for iComp = 1:length(complist) % loop each component
            
            % key parameter
            compName = complist{iComp};
            compWin  = twlist(iComp, :);
            binName  = ERP.bindescr{iBin};
            binName(binName=='_')='-';
            idx_start = find (ERP.times == (tw(1)));
            idx_end   = find (ERP.times == (tw(2)));
            
            % replace demo 
            erpdata   = ERP.bindata(1:length(data.label), idx_start:idx_end, iBin);
            data_upt  = data;
            data_upt.avg = erpdata;
            
            % plot topo
            cfg = [];
            cfg.xlim     = compWin;
            cfg.zlim  = 'maxmin';
            cfg.baseline = 'no';
            cfg.baselinetype = 'absolute';
            cfg.layout = 'EEG1005.lay'; 
            cfg.parameter = 'avg';
            cfg.colormap = '*RdBu';
            figure; ft_topoplotER(cfg,data_upt); colorbar
            title(strcat(curTask,'-', compName, '-', binName));
            
            % save image
            saveas(gcf,strcat('topo_', curTask,'-', compName, '-', binName, num2str(bl(1))),'epsc');
            saveas(gcf,strcat('topo_', curTask,'-', compName, '-', binName, num2str(bl(1))),'png'); 
            
            % clear everything
            data_upt = []; erpdata = [];
            close all;
        end
    end
    
    ERP = []; data = [];
    
end
