%% ReadME
% ERPlab Grand Average
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: Aug 29, 2021 


clear all; clc;
subNameD = [1205:1242,1244:1253];
%subNameD = 1205:1206;
badsubNameD   = [1201:1204, 1220, 1229, 1238, 1239, 1249]; % 1201 - 1204: Local, 1229, 1238, 1239,1249 (reject epochs > 100)
subNameD = subNameD(~ismember(subNameD, badsubNameD));

workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
%task = {'preLearningImplicit', 'Learning', 'postLearningImplicit'};
task = {'preLearningImplicit', 'postLearningImplicit'}

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'ERP'); % set your output path
inputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing')
% change accordingly to your eeglab version and location

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task)
    
    curTask = task{iTask};
    
    ERPsetList = [curTask 'GrandAvgAll_epoch.txt']; % create a txt file to save the subject level ERP sets to be involved in grand averge
    
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
        outputSubjectFolder = [outputParentFolder '/' subName];
        
        if ~exist(outputSubjectFolder)
            mkdir(outputSubjectFolder);
        end
        
        % load existing preprocessed dataset
        EEG = pop_loadset('filename',[subName '_' curTask '_EpochArtRej_FalseRespRej.set'],'filepath',fullfile(inputParentFolder, subName));
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        %Average
        ERP = pop_averager( ALLEEG , 'Criterion', 'good', 'ExcludeBoundary', 'on', 'SEM', 'on' );
        ERP = pop_erpchanoperator( ERP, '/home/chendanni/Documents/Norms/analysis/MyScripts/ERP/ROIs.txt',...
                                'ErrorMsg', 'popup', 'KeepLocations',  0, 'Warning', 'off' );
        ERP = pop_savemyerp(ERP, 'erpname',...
            [subName '_' curTask '_ERP_Avg'], 'filename', [subName '_' curTask '_ERPAvg.erp'], 'filepath', outputSubjectFolder, 'Warning','off');
        fprintf(fid,[outputSubjectFolder '/' subName '_' curTask '_ERPAvg.erp\n']);
        
        ALLEEG = [];
    
    end
    
    fclose(fid);
    
    % Grandaverage over ERP sets written in the txt file of ERPsetList.
    
    ERP = pop_gaverager( [outputParentFolder '/' curTask 'ERPlabGrandAvg/' ERPsetList] , 'SEM', 'on',...
        'Weighted', 'on' );
    ERP = pop_savemyerp(ERP, 'erpname', 'GrandAvgAll', 'filename', [curTask '_GrandAvg_ManuRej.erp'], 'filepath',...
        outputERPFolder, 'Warning', 'on');

end

%% Done.
