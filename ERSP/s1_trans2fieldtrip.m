% ERSP analysis
% step 1: transfer preprocessed EEG data to fieldtrip version (note: not
% baseline corrected)
% Project: Social Norms Learning
% Dependency: FieldTrip
% Author: Danni Chen 
% Update Date: 2022-03-01

%% basic setting

% Subject info
clear all; clc;
subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238]; % 1201 - 1204: HK Local; 1214, 1215 & 1238 - visually detected bad subjects
subNameD = subNameD(~ismember(subNameD, badsubNameD));

% folder info
workFolder = '/home/chendanni/Documents/Norms/analysis/';
task = {'preLearningImplicit', 'postLearningImplicit', 'Learning'};
PreprocessParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing');
ERSPParentFolder = fullfile (workFolder, 'EEGAnalysisResults', 'ERSP', 'Results');
ScriptsFolder = fullfile(workFolder, 'MyScripts');
FieldTripFolder = '/home/chendanni/MATLAB/toolbox/fieldtrip-20210330';
addpath ( genpath(ScriptsFolder) );
addpath ( genpath(FieldTripFolder) );


for iTask = 1:length(task)
    
    curTask = task{iTask};
    
    avg_in_high = {};
    avg_in_low  = {};
    avg_in_con  = {};
    avg_out_high = {};
    avg_out_low  = {};
    avg_out_con  = {};
    avg_all = {};
    
    % transfer from EEGlab to FieldTrip
    for iSubject = 1:length(subNameD)
        
        subName = num2str (subNameD(iSubject));
        
        subjectFolder = fullfile(PreprocessParentFolder, subName);
        cd (subjectFolder);
        
        % load dataset (not baseline corrected)
        switch curTask
            case 'Learning'
                fileName = [subName '_' curTask '_EpochArtRej.set'];
            case 'preLearningImplicit'
                fileName = [subName '_' curTask '_EpochArtRej_FalseRespRej.set'];
            case 'postLearningImplicit'
                fileName = [subName '_' curTask '_EpochArtRej_FalseRespRej.set'];
        end
        EEG = pop_loadset('filename',fileName,'filepath',pwd); 
        
        % transfer to fieldtrip  
        data = eeglab2fieldtrip (EEG, 'preprocessing', 'none');
        outputSubjectFolder = fullfile (ERSPParentFolder, subName);
        if ~exist(outputSubjectFolder)
            mkdir(outputSubjectFolder);
        end
        cd(outputSubjectFolder);
        save (strcat(subName,'_',curTask,'_processed'), 'data');
        
        for k = 1:6
            cfg = []; task1 = [];
            cfg.trials = find (data.trialinfo.bini == k);
            task1 = ft_timelockanalysis(cfg, data);
            
            switch k 
                case 1
                    avg_in_high{iSubject} = task1;
                case 2
                    avg_in_low{iSubject}  = task1;
                case 3
                    avg_in_con{iSubject}  = task1;
                case 4
                    avg_out_high{iSubject} = task1;
                case 5
                    avg_out_low{iSubject} = task1;
                case 6
                    avg_out_con{iSubject} = task1;
            end
        end
        avg_all{iSubject} = data;
        
        data = []; EEG = [];
    end
    
    ERSPtaskFolder = fullfile(ERSPParentFolder, strcat(curTask, 'Avg'));
    if ~exist(ERSPtaskFolder)
        mkdir(ERSPtaskFolder);
    end
    cd (ERSPtaskFolder);
    
    save (strcat(curTask,'_processed_all'),...
        'avg_all', 'avg_in_high', 'avg_in_low', 'avg_in_con',...
         'avg_out_high', 'avg_out_low', 'avg_out_con', '-v7.3');
    
end
