%% pre-processing EEG data
% s4: epoch in ERPLab
% Project: Social Norms Learning
% Author: Danni Chen 
% Update Date: 2022-02-16 
% Input: *_EpochRej.set
% Output: *_EpochRejICA.set
% Parameter:
%   1) subNameD: Vector, a list of subject number
%   2) task_list
%   3) workFolder
%   4) binFolder
% Note that: Run ICA to remove EOG (be cautious when removing components other than blink artifacts)
% Warning: This step might take quite long time depends on your data size and computer functionality. 
% After this step, load the .set file to manually check ICA components and
% remove EOG or other artefacts. One subject at a time.

function s6_ica(subNameD, task_list, workFolder, sr, lp, hp, notch, linefreq, channelLocationFile, gate)

if nargin == 0

    % participant info set up
    subNameD = [1201:1242, 1244:1253];
    
    % directory set up
    workFolder = '/Users/danni/Desktop/Norms/'; % in MAC

    % task set up
    task_list = {'preLearningImplicit','postLearningImplicit','Learning'};
    
    % parameter set up
    sr = 250;
    lp = 30;
    hp = 0.05;
    notch = 1;
    linefreq = 50;
    
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
        
        if exist (strcat(subName, '_', curTask, '_EpochRejICA.set'))
            fprintf (strcat('existed: ',  strcat(subName, '_', curTask, '_EpochRejICA.set')), ' and skipped');
            fprintf ('\n');
            continue;
        end % just skip this participants
        
        
        if gate == 1
            
            % load Epoched dataset saved in s4
            EEG = pop_loadset('filename', strcat(subName, '_', curTask, '_EpochRej.set'),'filepath',pwd ); % load reject big artifacts ones 
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

            % Use high-pass filter at 1 hz preprocessed data to do ICA (this facilitates ICA caculation)
            EEG = pop_eegfiltnew(EEG, 'locutoff',1,'plotfreqz',0);            
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', strcat(subName, '_EpochRej1Hz'),...
                'overwrite','on','gui','off');
            
            TMP = EEG;
            clear EEG;

        elseif gate == 2
            
            TMP = [];
            TMP = epoch1Hz (subNameD(iSubject), curTask, workFolder, sr, lp, hp, notch, linefreq);
            
        end
        
        % Run ICA. Usually it takes 30 minutes for each subject
        TMP = pop_runica(TMP, 'icatype','runica', 'extended',1, 'pca', TMP.nbchan-1); % 'pca', [num] run how many components 
        [ALLEEG TMP] = eeg_store(ALLEEG, TMP, CURRENTSET); 
        TMP = eeg_checkset( TMP );
        pop_saveset( TMP, [ subName '_' curTask '_EpochRejICA1Hz.set'], pwd);
        floatwrite(TMP.icaweights, [outputSubjectFolder filesep subName '_' curTask '_ICAweight.wts']);
        floatwrite(TMP.icasphere, [outputSubjectFolder filesep subName '_' curTask '_sphere.sph']);

        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[]; TMP = [];  %clear all
                
        % load dataset before highpass 1 Hz (which was saved in s4, 0.1 Hz highpassed)
        EEG = pop_loadset('filename',[ subName '_' curTask '_EpochRej.set'],...
            'filepath',[outputSubjectFolder filesep]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',[subName '_EpochRej.set'],'gui','off');
        
        % apply the ICA defenitions to the 0.1 hz high-passed set
        EEG = pop_editset(EEG, 'icachansind', [],'icaweights',...
            [outputSubjectFolder filesep subName '_' curTask '_ICAweight.wts'],...
            'icasphere',[outputSubjectFolder filesep subName '_' curTask '_sphere.sph']);
        
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        EEG = eeg_checkset(EEG);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',...
            [ subName '_' curTask '_EpochRejICA'],...
            'savenew',[subName '_' curTask '_EpochRejICA'],...
            'overwrite','on','gui','off'); 

    end

end

end