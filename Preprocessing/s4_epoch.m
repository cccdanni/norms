%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Input: *_InterReRef.set (remove bad channel); *_binassigned.txt
% Output: *_Epoch.set
% Update Date: July-12-2021
% Last run date: July-12-2021

%% Step 4: Epoch for ERP in ERPlab
% You need to preprare your binlist before running, and store them under 
% binlistFolder.

clear all; clc;

subNameD = [1201:1242,1244:1253]; 
workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 

%task = {'preLearningImplicit', 'Learning', 'postLearningMemory', 'postLearningImplicit'};
%task = {'preLearningImplicit', 'postLearningImplicit'};
task = {'postLearningImplicit'};


% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path
binlistFolder = '/home/chendanni/Documents/Norms/analysis/Binlist'; % path to your binlist
binassname = '_binassigned.txt';

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
    EEG = pop_loadset('filename',[subName '_' curTask '_InterReRef.set'],'filepath',outputSubjectFolder);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    
    % Trigger/Binlist file name, if you have counterbalancing, set
    % accordingly
    switch curTask
        case 'preLearningImplicit' 
            binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '_implicit.txt'];
            for i = 1:(size(EEG.urevent,2)-624-1)
               EEG.urevent(i).type = 'practice';
               EEG.event(i).codelabel = 'practice';
               EEG.event(i).type = 199;
            end
        case 'Learning'
            binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '.txt'];
        case 'postLearningMemory'
            binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '.txt'];
        case 'postLearningImplicit'
            binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '_implicit.txt'];
            for i = 1:(size(EEG.urevent,2)-624-1)
               EEG.urevent(i).type = 'practice';
               EEG.event(i).codelabel = 'practice';
               EEG.event(i).type = 199;
            end
    end
    
    % create eventlist
    EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, ...
        'BoundaryString', { 'boundary' }, 'Eventlist', [outputSubjectFolder '/' subName  '_' curTask '_eventlist.txt']); % GUI: 05-Jun-2018 11:00:10
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName 'InterReRef.set'],'gui','off');
    
    % Assign Binlister
    % binfile = strcat(binlistFolder,filesep,binlistfile)
    binfile = strcat(binlistFolder,filesep,binlistfile);
    outputeventfile = [outputSubjectFolder filesep subName '_' curTask binassname];
    inputeventfile = [outputSubjectFolder filesep subName '_' curTask '_eventlist.txt'];
    EEG  = pop_binlister( EEG , 'BDF', binfile, 'ExportEL',...
        outputeventfile, 'ImportEL',...
        inputeventfile, 'IndexEL',  1, 'SendEL2', 'EEG&Text', 'Voutput', 'EEG'  );
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
   
    % extract bin-based epochs
    switch curTask
        case 'preLearningImplicit'
            EEG = pop_epochbin( EEG , [-200.0  1000.0],  'pre'); % set the epoch length and baseline correction window
        case 'Learning'
            EEG = pop_epochbin( EEG , [-300.0  3000.0],  'pre'); % set the epoch length and baseline correction window
        case 'postLearningMemory'
            EEG = pop_epochbin( EEG , [-200.0  1000.0],  'pre'); % set the epoch length and baseline correction window
        case 'postLearningImplicit'
            EEG = pop_epochbin( EEG , [-200.0  1000.0],  'pre'); % set the epoch length and baseline correction window
    end
    
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',...
        [outputSubjectFolder '/' subName '_' curTask '_Epoch.set'],'gui','off');
   
end

end

%% Before Run s5 ICA
% Load each exported data and reject trials with big artfacts