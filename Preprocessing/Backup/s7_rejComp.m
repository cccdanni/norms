%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Input: *_EpochRejICA.set
% Output: *__EpochRejICARej
% Update Date: July-15-2021
% Last run date: July-15-2021

% record: 1223 (pre), 1248 (post)

%% Step 6: Manually check ICA components to remove EOG or other artefacts 
% Save the datasets as "...rej.set"
clear; clc;
subNameD = [1248]; 
% subName = [1216:1223] 
workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
%task = {'preLearningImplicit', 'Learning', 'postLearningMemory', 'postLearningImplicit'};
task = {'postLearningImplicit'}

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path

% Load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task)
    curTask = task{iTask};
    
    for iSubject = 1:length(subNameD)
        
        % number to string
        subName = num2str(subNameD(iSubject));
        
        % locate sub folder 
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);
        
        % load Epoched dataset saved in s5
        EEG = pop_loadset('filename', [subName '_' curTask '_EpochRejICA.set'],...
            'filepath',[outputSubjectFolder '/']);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        %eeglab redraw;
        
        % Pop-out IC 
        EEG = pop_iclabel(EEG, 'Default');
        pop_viewprops( EEG, 0, [1:30], {'freqrange', [2 80]}, {}, 1, 'ICLabel' )
        
        % Input reject component
        reject_component_list = [];
        cnt = 1;
        while 1
            x = input('Input an IC number: \n');
            if isempty(x) break; 
            else reject_component_list(cnt) = x; cnt = cnt + 1;
            end
        end
        save (strcat(subName,curTask,'_rejcomp_list'),'reject_component_list');
        
        % reject component
        EEG = pop_subcomp( EEG, reject_component_list , 0);
        
        % save dataset
        EEG.setname = [subName '_' curTask '_EpochRejICARej'];
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',...
            [outputSubjectFolder '/' subName '_' curTask '_EpochRejICARej'],'gui','off');
    end
    
end

% Then run S7 to detect and reject artefact in epoch.

