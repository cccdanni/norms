%% ReadME
% Evaluate ERP Data Quality & Extract Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: 2022-02-27

%% basic setting

% Subject info
clear all; clc;
subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238, 1210, 1211, 1223, 1239]; 
% 1201 - 1204: HK Local; 
% 1214, 1215 & 1238 - visually detected bad subjects
% 1210, 1211, 1223, 1239 - removed more than 100 trials
subNameD = subNameD(~ismember(subNameD, badsubNameD));

% folder info
workFolder = '/home/chendanni/Documents/Norms/analysis/';
% task = {'preLearningImplicit', 'postLearningImplicit', 'Learning'};
task = {'Learning'};
ERPParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'ERP'); % set your output path
PreprocessParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing');
ScriptsFolder = '/home/chendanni/Documents/Norms/analysis/MyScripts/ERP';
addpath ( ScriptsFolder );

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

bl = [-300 0];

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

% other setting
fs = 250;
admtw = 0.100;

% type 
typelist = {'Mean', 'AdaptiveMean', 'Peak'};

%% extract component

ERPDataQuality = struct;
ERPCompData    = struct;
subAll = []; % initial data 

cntDQ = 1; % initial count 

for iSubject = 1: length(subNameD) % loop each subject

    subName = num2str(subNameD(iSubject));
    
    colNames{1} = 'SubjectID';
    colNames{2} = 'Epoch';
    colNames{3} = 'bini';
    colNames{4} = 'task';
    cntNr = 5;
    
    thisSub = [];
    
    for iTask = 1: length(task) % loop each task
        
        curTask = char (task(iTask));
        
        thisSub = [];
        
        %% load data
        EEG = []; ALLEEG = []; CURRENTSET = []; ERP = [];
        EEG = pop_loadset('filename', strcat(subName, '_', curTask, '_rmbl', num2str(bl(1)), '.set'),...
                          'filepath', fullfile(ERPParentFolder, subName));
        ERP = pop_loaderp('filename', strcat(subName, '_', curTask, '_ERPAvg', num2str(bl(1)), '.erp'),...
                          'filepath', fullfile(ERPParentFolder, subName));
        
        event = [EEG.event.bepoch];
        bini = [EEG.event.bini];
        [event, IA, IC] = unique(event); % remove replicate event number
        bini = bini(IA);
        thisSub = [thisSub ones(length(event),1)*subNameD(iSubject)];
        thisSub = [thisSub event'];
        thisSub = [thisSub bini'];
        thisSub = [thisSub ones(length(event),1)*iTask];
                      
        %% data quality
        ERPDataQuality(cntDQ).SubjectID = subNameD(iSubject);
        thisDQ = ERP.dataquality;
        aSMEnr = find(strcmp({thisDQ.type},'aSME'));
        ERPDataQuality(cntDQ).type  = thisDQ(aSMEnr).type;
        ERPDataQuality(cntDQ).times = thisDQ(aSMEnr).times;
        ERPDataQuality(cntDQ).data  = thisDQ(aSMEnr).data;
        cntDQ = cntDQ + 1;
        
        % component of interest info
        switch curTask
            case 'Learning'
                bltw = bl;
                complist = {'n170','lpc','n400','frn','epn','lpp','p300'};
                twlist = [0.12 0.22;...
                  0.30 0.80;...
                  0.30 0.50;...
                  0.25 0.45;...
                  0.23 0.28;...
                  1.00 1.60;...
                  0.20 0.40];
                roilist = {{'LOT','ROT','OT'},...
                   {'CP','FC','LP','RP'},...
                   {'Ce','CP'},...
                   {'FC'}...
                   {'LOT','ROT','OT'}...
                   {'LP', 'RP', 'CP', 'FC'}...
                   {'LP', 'RP', 'CP'}};
            case 'preLearningImplicit'
                bltw = [-0.20 0];
                complist = {'n170','lpc','n400','frn','epn'};
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
                complist = {'n170','lpc','n400','frn','epn'};
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
        
        %% component extraction
        for iComp = 1:length(complist) % loop all components
            
            thisCompName = complist{iComp};
            thisCompROI = roilist{iComp};
            tw = twlist(iComp, :);
            
            for iROI = 1:length(thisCompROI) % loop all ROIs
                
                thisROIName = thisCompROI{iROI};
                ROIChans = ROI(find(strcmp({ROI.ROIName},thisROIName))).ChanList;
                
                for itype = 1:length(typelist) % loop all types
                    
                    %% record type
                    thistypeName = typelist{itype};
                    thisColName = strcat(thisCompName, thisROIName, thistypeName);
                    colNames{cntNr} = thisColName;
                    
                    %% extract comp
                    CompVol = CompExtraction(EEG, tw, bltw, fs, ROIChans, thistypeName, admtw); % NOT do baseline correction
                    thisSub = [thisSub, CompVol'];
                    
                    cntNr = cntNr + 1;
                    
                end
                
            end
            
        end
        
        subAll = [subAll; thisSub];
    
    end

end

cd (ERPParentFolder);
T = array2table(subAll);
colNames = colNames(1:size(subAll,2));
T.Properties.VariableNames = colNames;
tName = convertCharsToStrings(char(task));
writetable(T, strcat(tName, 'ERPComp',num2str(bl(1)),date(),'.xlsx'));