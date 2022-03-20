%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Input: *_EpochRej.set
% Output: *_EpochRejICA.set
% Update Date: July-13-2021
% Last run date: July-13-2021

%% Step 5: Run ICA to remove EOG (be cautious when removing components other than blink artifacts)
% Warning: This step might take quite long time depends on your data size and computer functionality. 
% After this step, load the .set file to manually check ICA components and
% remove EOG or other artefacts. One subject at a time.

clear all; clc;

% subNameD = [1216:1242, 1244:1253]; 
subNameD = 1248;
workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
%task = {'preLearningImplicit', 'Learning', 'postLearningMemory', 'postLearningImplicit'};
task = {'postLearningImplicit'};

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults','Preprocessing'); % set your output path

% Load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

err_cnt = 0;

for iTask = 1:length(task)
    curTask = task{iTask};
    
for iSubject = 1:length(subNameD)
    
    try
        
        % number to string
        subName = num2str(subNameD(iSubject));
        
        % locate sub folder 
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);
        
%         % load Epoched dataset saved in s4
%         EEG = pop_loadset('filename',[ subName '_' curTask '_EpochRej.set'],... %load reject big artifacts ones (DC, 03-17-21)
%             'filepath',[outputSubjectFolder '/']);
%         [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
%         % high-pass filter at 1 hz 
%         EEG = pop_eegfiltnew(EEG, [], 1, [], true, [], 1);
%         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[ subName...
%             'EpochHP1Hz'],'overwrite','on','gui','off'); 

        TMP = [];
        
        % Use high-pass filter at 1 hz preprocessed data to do ICA (this facilitates ICA caculation) 
        TMP = epoch1Hz (subNameD(iSubject), curTask);
        
        % Run ICA. Usually it takes 30 minutes for each subject
        TMP = pop_runica(TMP, 'icatype','runica', 'extended',1, 'pca', TMP.nbchan-1);% 'pca', [num] run how many components 
        [ALLEEG TMP] = eeg_store(ALLEEG, TMP, CURRENTSET); 
        TMP = eeg_checkset( TMP );
        pop_saveset( TMP, [ subName '_' curTask '_EpochRejICA1Hz.set'], pwd);
        floatwrite(TMP.icaweights,[outputSubjectFolder filesep subName '_' curTask '_ICAweight.wts']);
        floatwrite(TMP.icasphere,[outputSubjectFolder filesep subName '_' curTask '_sphere.sph']);
        TMP = pop_loadset('filename',[ subName '_' curTask '_EpochRejICA1Hz.set'],...
                'filepath',[outputSubjectFolder filesep]);

        
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );%clear all
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];  %clear all
        
        % load dataset before highpass 1 Hz (which was saved in s4, 0.1 Hz highpassed)
        EEG = pop_loadset('filename',[ subName '_' curTask '_EpochRej.set'],...
            'filepath',[outputSubjectFolder filesep]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, TMP, 1,'setname',[subName '_EpochRej1Hz.set'],'gui','off');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[subName '_EpochRej.set'],'gui','off');
        
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
 
    catch
        err_cnt = err_cnt + 1;
    end

end

end

eeglab redraw

% Then, manually check ICA components to remove EOG or other artefacts. 
%  Save the data sets as "...rej.set". 
%  Then run S7 to detect and reject artefact in epoch.


eeglab redraw

% Next step (s4): Epoch
