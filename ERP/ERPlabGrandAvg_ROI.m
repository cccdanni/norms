%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: Jan-13-2021
% Mostly revised based on EEG Preprocessing Scripts of Hu Lab


clear all; clc;

subNameD = [1205:1220]; 
badsubNameD   = [1201:1204, 1228, 1237, 1239, 1229]; % 1201 - 1204: Local, 1228, 129, 1237, 1239 (remained less than 50%)
subNameD = subNameD(~ismember(subNameD, badsubNameD));

workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
%task = {'preLearningImplicit', 'Learning', 'postLearningImplicit'};
task = {'preLearningImplicit', 'postLearningImplicit'}

% change outputParentFolder before ERP OR ERSP
outputParentFolder = [workFolder,'EEGAnalysisResults', 'ERP']; % set your output path
% change accordingly to your eeglab version and location

ROIList = struct;
ROIs = {'FC','Ce','CP','LP','RP','BiP','LOT','ROT','OT'}; 
for i = 1:8 ROIList(i).ROIName = ROIs{i}; end
ROIList(1).ChanList = {'Fz', 'FCz', 'F1', 'F2', 'FC1', 'FC2'}; %FC
ROIList(2).ChanList = {'Cz', 'C1', 'C2'}; %Ce
ROIList(3).ChanList = {'CP1', 'CP2', 'Pz', 'P1', 'P2'}; %CP
ROIList(4).ChanList = {'P3', 'P5', 'CP3', 'CP5'}; %LP
ROIList(5).ChanList = {'P4', 'P6', 'CP4', 'CP6'}; %RP
ROIList(6).ChanList = {'T7','TP7','P7','PO7'}; %LOT
ROIList(7).ChanList = {'T8','TP8','P8','PO8'}; %ROT
ROIList(8).ChanList = {'T7','TP7','P7','PO7','T8','TP8','P8','PO8'}; %OT

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iROI = 1:length(ROIList)
    
    thisROI = ROIList(iROI).ROIName;
    thisROIChans = ROIList(iROI).ChanList;

for iTask = 1:length(task)
    curTask = task{iTask};
    
    ERPsetList = [curTask thisROI 'GrandAvgAll_epoch.txt']; % create a txt file to save the subject level ERP sets to be involved in grand averge
    
    outputERPFolder = [outputParentFolder '/' curTask 'ERPlabGrandAvg'];
    if ~exist(outputERPFolder)
        mkdir(outputERPFolder)
    end
    
    cd(outputERPFolder);
    
    fid = fopen(ERPsetList,'w');
    
    for iSubject = 1:length(subNameD)
        
        %number to string
        subName = num2str(subNameD(iSubject));
        
        %locate sub folder
        outputSubjectFolder = [outputParentFolder subName];
        
        % load existing preprocessed dataset
        EEG = pop_loadset('filename',[subName '_' curTask '_EpochArtRej_FalseRespRej.set'],'filepath',outputSubjectFolder);
        % select ROI channels
        EEG = pop_select( EEG, 'channel', thisROIChans );
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        %Average
        ERP = pop_averager( ALLEEG , 'Criterion', 'good', 'ExcludeBoundary', 'on', 'SEM', 'on' );
        ERP = pop_savemyerp(ERP, 'erpname',...
            [subName '_' curTask thisROI '_ERP_Avg'], 'filename', [subName '_' curTask thisROI '_ERPAvg.erp'], 'filepath', outputSubjectFolder, 'Warning','off');
        
        fprintf(fid,[outputSubjectFolder '/' subName '_' curTask thisROI '_ERPAvg.erp\n']);
        
        ALLEEG = [];
    
    end
    
    fclose(fid);
    
    % Grandaverage over ERP sets written in the txt file of ERPsetList.
    
    ERP = pop_gaverager( [outputParentFolder '/' curTask 'ERPlabGrandAvg/' ERPsetList] , 'SEM', 'on',...
        'Weighted', 'on' );
    ERP = pop_savemyerp(ERP, 'erpname', 'GrandAvgAll', 'filename', [curTask thisROI '_GrandAvg_ManuRej.erp'], 'filepath',...
        outputERPFolder, 'Warning', 'off');

end

end

%% Done.
