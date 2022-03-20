% Plot ERP topography with FieldTrip
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: 8/30/2021

%% basic setting

% Subject info
clear all; clc; close all;
subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238]; % 1201 - 1204: HK Local; 1214, 1215 & 1238 - visually detected bad subjects
subNameD = subNameD(~ismember(subNameD, badsubNameD));

% folder info
workFolder = '/home/chendanni/Documents/Norms/analysis/';
task = {'post-minus-pre'};
ERPParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'ERP'); % set your output path
PreprocessParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing');
ScriptsFolder = '/home/chendanni/Documents/Norms/analysis/MyScripts/ERP';
FieldTripFolder = '/home/chendanni/MATLAB/toolbox/fieldtrip-20210330';
addpath ( genpath(ScriptsFolder) );
addpath ( genpath(FieldTripFolder) );
addpath ( genpath("/home/chendanni/Documents/Norms/analysis/MyScripts/additional_functions/boundedline-pkg-master") );

%% Plot ERP topo

for iTask = 1:length(task)
    
    curTask = task{iTask};
    
    % component of interest info
    switch curTask
        case 'Learning'
            bltw = [-0.30 0];
            tw   = [-300 2996];
            complist = {'n170','lpc','n400','frn','epn','lpp'};
            demoName = 'demo_learning.mat';
            twlist = [0.12 0.22;...
              0.30 0.80;...
              0.30 0.50;...
              0.25 0.45;...
              0.23 0.28;...
              1.00 1.60];
            roilist = {{'LOT','ROT','OT'},...
               {'CP','FC','LP','RP'},...
               {'Ce','CP'},...
               {'FC'}...
               {'LOT','ROT','OT'}...
               {'LP', 'RP', 'CP', 'FC'}};
        case 'preLearningImplicit'
            bltw = [-0.20 0];
            tw   = [-200 996];
            complist = {'n170','lpc','n400','frn','epn'};
            demoName = 'demo_implicit.mat';
            twlist = [0.12 0.22;...
              0.30 0.80;...
              0.30 0.50;...
              0.25 0.45;...
              0.23 0.28];
            roilist = {{'LOT','ROT','OT'},...
               {'CP','FC','LP','RP'},...
               {'Ce','CP'},...
               {'FC'}...
               {'LOT','ROT','OT'}};
        case 'postLearningImplicit'
            bltw = [-0.20 0];
            tw   = [-200 996];
            complist = {'n170','lpc','n400','frn','epn'};
            demoName = 'demo_implicit.mat';
            twlist = [0.12 0.22;...
              0.30 0.80;...
              0.30 0.50;...
              0.25 0.45;...
              0.23 0.28];
            roilist = {{'LOT','ROT','OT'},...
               {'CP','FC','LP','RP'},...
               {'Ce','CP'},...
               {'FC'}...
               {'LOT','ROT','OT'}};
        case 'post-minus-pre'
            bltw = [-0.20 0];
            tw   = [-200 996];
            complist = {'n170','lpc','n400','frn','epn'};
            demoName = 'demo_implicit.mat';
            twlist = [0.12 0.22;...
              0.30 0.80;...
              0.30 0.50;...
              0.25 0.45;...
              0.23 0.28];
            roilist = {{'LOT','ROT','OT'},...
               {'CP','FC','LP','RP'},...
               {'Ce','CP'},...
               {'FC'}...
               {'LOT','ROT','OT'}};    
    end
    
    
    % load demo 
    cd(ScriptsFolder);
    data = [];
    load (demoName);
    
    % load ERP data
    PreERPName = strcat ('preLearningImplicit_GrandAvg_ManuRej.erp'); 
    PreERPFolder = [ERPParentFolder '/preLearningImplicitERPlabGrandAvg'];
    PreERP = pop_loaderp('filename', PreERPName, 'filepath', PreERPFolder);
    
    PostERPName = strcat ('postLearningImplicit_GrandAvg_ManuRej.erp'); 
    PostERPFolder = [ERPParentFolder '/postLearningImplicitERPlabGrandAvg'];
    PostERP = pop_loaderp('filename', PostERPName, 'filepath', PostERPFolder);
    
    ERP = PostERP;
    ERP.bindata = PostERP.bindata - PreERP.bindata;
    ERP.filename = 'post-minus-pre';
    ERPFolder = [ERPParentFolder '/post-minus-preERPlabGrandAvg'];
    if ~exist (ERPFolder)
        mkdir (ERPFolder);
    end
    cd (ERPFolder);
    
    
    clear PreERPName PostERPName PreERP PostERP;
    
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
            saveas(gcf,strcat('topo_', curTask,'-', compName, '-', binName),'epsc');
            saveas(gcf,strcat('topo_', curTask,'-', compName, '-', binName),'png'); 
            
            % clear everything
            data_upt = []; erpdata = [];
            close all;
        end
    end
    
    ERP = []; data = [];
    
end
