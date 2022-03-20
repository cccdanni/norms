%% pre-processing EEG data
% s7: Manually check ICA components to remove eye- or muscle-movement
% artifacts
% Project: Social Norms Learning
% Author: Danni Chen 
% Update Date: 2022-02-22
% Input: *_EpochRejICA.set
% Output: *__EpochRejICARej
% Parameter:
%   1) subNameD: Vector, a list of subject number
%   2) task_list
%   3) workFolder


function s7_rejComp(subNameD, task_list, workFolder)

if nargin == 0

    % participant info set up
    subNameD = [1201:1242, 1244:1253];
    
    % directory set up
    workFolder = '/Users/danni/Desktop/Norms/'; % in MAC

    % task set up
    task_list = {'preLearningImplicit','postLearningImplicit','Learning'};
    
end

outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path

% Load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task_list)
    
    curTask = task_list{iTask};
    
    for iSubject = 1:length(subNameD)
        
        % number to string
        subName = num2str(subNameD(iSubject));
        
        % locate sub folder 
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);
        
        if exist (strcat(subName, '_', curTask, '_EpochRejICARej.set'))
            fprintf (strcat('existed: ',  strcat(subName, '_', curTask, '_EpochRejICARej.set'), ' and skipped'));
            fprintf ('\n');
            continue;
        end % just skip this participants
        
        % load Epoched dataset saved in s6 (ICA)
        EEG = pop_loadset('filename', [subName '_' curTask '_EpochRejICA.set'],'filepath',pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        %eeglab redraw;
        
        % Pop-out IC 
        EEG = pop_iclabel(EEG, 'Default');
        pop_viewprops( EEG, 0, 1:30, {'freqrange', [2 80]}, {}, 1, 'ICLabel' )
        
        % Input reject component
        reject_component_list = [];
        cnt = 1;
        while 1
            x = input('Input an IC number: \n');
            if isempty(x) break; 
            else reject_component_list(cnt) = x; cnt = cnt + 1;
            end
        end
        save (strcat(subName,'_', curTask,'_rejcomp_list'),'reject_component_list');
        fprintf(strcat('reject component: ', num2str(reject_component_list)));
        fprintf ('\n');
        
        % reject component
        EEG = pop_subcomp( EEG, reject_component_list , 0);
        
        % save dataset
        EEG.setname = [subName '_' curTask '_EpochRejICARej'];
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',...
            [outputSubjectFolder '/' subName '_' curTask '_EpochRejICARej'],'gui','off');
        
        close all;
    end
    
end

% Then run S7 to detect and reject artefact in epoch.


end
