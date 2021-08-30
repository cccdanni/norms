%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Input: *_EpochRejICARej.set
% Output: *_EpochArtDet; *_EpochArtRej
% Update Date: July-15-2021
% Last run date: July-15-2021


%% Step 7: Detect artefacts in Epoch in ERPlab
% This step generates two .set files for each subject. One has a suffix of
% "EpochArtDet" which marks out trials with artifacts but maintain all of
% them. The other has a suffix of "EpochArtRej" which remove/reject the
% trials with artifacts. 
% An excel file called "RejectTrialNumber" will be saved under
% outputParentFolder. You are advised to check this sheet to see how many 
% trials are rejected for each subject. If the number is too large, load 
% the "...EpochArtDet.set" file to manually check the trials marked with 
% artifacts and reject manually. Then save your manual-rejected version as 
% "...EpochArtRej.set" for later analysis.
% 
% Authors: Yao Ziqing, Lin Xuanyi, Hu Xiaoqing
% Create Date: 2018.07
% Update Date: 2020.07

clear; clc;

% subNameD = [1205:1220];
subNameD = [1205:1242,1244:1253]; 
workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
%task = {'preLearningImplicit', 'Learning', 'postLearningMemory', 'postLearningImplicit'};
task = {'preLearningImplicit', 'postLearningImplicit'};

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path

%load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task)
    curTask = task{iTask};
    
    for iSubject = 1:length(subNameD)

        %number to string
        subName = num2str(subNameD(iSubject));
        
        %locate sub folder
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);
        
        % load existing preprocessed dataset
        EEG = pop_loadset('filename',[subName '_' curTask '_EpochRejICARej.set'],'filepath',outputSubjectFolder);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
        
        % load file saved in s1 before removing bad channels
        EEG = pop_loadset('filename',[ subName '_' curTask '_250Hz05HP30LPNotch.set'],...
            'filepath',[outputSubjectFolder]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
        
        % interpolate channels on dataset 1, based on dataset 2(dataset from step 1),
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0);
        EEG = eeg_checkset( EEG );
        EEG = pop_interp(EEG, ALLEEG(2).chanlocs, 'spherical');
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',...
            [subName '_' curTask '_EpochRejICARej_inter.set'],'overwrite','on','gui','off');
        
        % Artefact Detection - Moving window peak to peak. Set your Twindow as
        % your epoch length
        % The Threshold used here is 100uv. 
        switch curTask
            case 'preLearningImplicit'
                lwindow = [-200 996];
            case 'Learning'
                lwindow = [-300 2996];
            case 'postLearningImplicit'
                lwindow = [-200 996];
            case 'postLearningMemory'
                lwindow = [-200 996];
        end
        

        EEG = pop_artmwppth( EEG , 'Channel',  1:61, 'Flag',  1, 'Threshold',  100, 'Twindow', lwindow,...
        'Windowsize',  200, 'Windowstep',100 ); % GUI: 15-Jul-2021 23:25:26
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off');
        
        EEG = pop_saveset( EEG, 'filename', [ subName '_' curTask '_EpochArtDet'],...
            'filepath', [outputSubjectFolder '/']);
        
        RejNum(iSubject,1) = subNameD(iSubject);
        RejNum(iSubject,2) = size(find(EEG.reject.rejmanual==1),2);
        
        % Reject detected artefacts
        EEG = eeg_checkset( EEG );
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, find(EEG.reject.rejmanual==1),0);
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'savenew',...
            [outputSubjectFolder filesep subName '_' curTask '_EpochArtRej'],'gui','off');
   
    end
    
    cd(outputParentFolder);
    xlswrite([curTask 'RejectTrialNumber75' date], RejNum);
    
    disp(['Completed: ' curTask]);
    
    ALLEEG = []; EEG = []; CURRENTSET = [];

end

%% DONE!