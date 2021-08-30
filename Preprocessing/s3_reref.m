%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Input: *_250Hz05HP30LPNotch.set
% Output: *_Inter.set (did not remove bad channel); *_InterReRef.set (remove bad channel) 
% Update Date: July-12-2021
% Last run date: July-12-2021

%% Step 3: Before ICA,
% 1) Remove bad channels saved in s2;
% 2) Then implement interpolation before re-referece to common average;
% 3) Remove interpolated channels before epoch & ICA.

clear all; clc;

subNameD = [1201:1242,1244:1253]; 
workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
%task = {'preLearningImplicit', 'Learning', 'postLearningMemory', 'postLearningImplicit'};
%task = {'preLearningImplicit', 'postLearningImplicit'};
task = {'postLearningImplicit'};


% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults','Preprocessing'); % set your output path
channelLocationFile = 'C:\toolbox\eeglab2019_1\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp';% Channel location file(path according to your file location)
% change accordingly to your eeglab version and location

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task)
    curTask = task{iTask};
    
for iSubject = 1:length(subNameD)
    % number to string
    subName = num2str(subNameD(iSubject));
    
    % locate sub folder
    outputSubjectFolder = fullfile(outputParentFolder,subName);
    cd(outputSubjectFolder);
    
    % load badSection (must have this file in the same subject folder
    % load([ subName 'badSection.mat']);
    % badSection = TMPREJ(:,1:2);
    load([curTask '_' subName '_badChan.mat']);
    
    % load file saved in s1 before removing bad sections/channels
    EEG = pop_loadset('filename',[ subName '_' curTask '_250Hz05HP30LPNotch.set'],...
        'filepath',[outputSubjectFolder '/']);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
    
    % %remove bad section
    % EEG = eeg_eegrej( EEG, badSection);
    % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',...
    %     [ subName '250Hz01BP30BadSec']...
    %     ,'overwrite','on','gui','off');
    % EEG = eeg_checkset( EEG );
    
    % remove bad channel
    EEG = pop_select(EEG, 'nochannel', badChan);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',...
        [ subName '_250Hz005HP30LPNotchBadChan']...
        ,'overwrite','on','gui','off');
    
    % interpolate channels on dataset 3, based on dataset 1(dataset from step 1),
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',3,'study',0);
    EEG = eeg_checkset( EEG );
    EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname',...
        [subName '_Inter.set'],'overwrite','on','gui','off');
    
    % Re-referece to common average
    EEG = pop_reref( EEG, []);
    
    % Remove bad channel again, then save set before next step
    EEG = pop_select(EEG, 'nochannel', badChan);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'setname',...
        [subName '_' curTask '_InterReRef.set'],...
        'savenew', [subName '_' curTask '_InterReRef.set'],'gui','off');
    
end
end

eeglab redraw

% Next step (s4): Epoch

