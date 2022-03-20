%% pre-processing EEG data
% s7: Manually check ICA components to remove eye- or muscle-movement
% artifacts Project: Social Norms Learning Author: Danni Chen Update Date:
% 2022-02-22 Input: *_EpochRejICARej.set Output: *_EpochArtDet;
% *_EpochArtRej Parameter:
%   1) subNameD: Vector, a list of subject number
%   2) task_list
%   3) workFolder
%   4) inter_option: boolean, true or false, if you want to interpolate bad
%   channel in the current function
%   5) threshold: 75 or 100, depend on the data quality
% Note:
%   1) One has a suffix of "EpochArtDet" which marks out trials with
%   artifacts but maintain all of them. The other has a suffix of
%   "EpochArtRej" which remove/reject the trials with artifacts. 2) An
%   excel file called "RejectTrialNumber" will be saved under
%   outputParentFolder. You are advised to check this sheet to see how many
%   trials are rejected for each subject. If the number is too large, load
%   the "...EpochArtDet.set" file to manually check the trials marked with
%   artifacts and reject manually. Then save your manual-rejected version
%   as "...EpochArtRej.set" for later analysis.


function s8_rejArt(subNameD, task_list, workFolder, inter_option, threshold, sr, lp, hp, win_size, win_step)

if nargin == 0

    % participant info set up
    subNameD = [1201:1242, 1244:1253];
    
    % directory set up
    workFolder = '/Users/danni/Desktop/Norms/'; % in MAC

    % task set up
    task_list = {'preLearningImplicit','postLearningImplicit','Learning'};
    
    % parameter set-up
    inter_option = true;
    threshold = 100;
    sr = 250;
    lp = 30;
    hp = 0.05;
    win_size = 200;
    win_step = 100;
    
end

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path

%load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task_list)
    
    curTask = task_list{iTask};
    
    
    for iSubject = 1:length(subNameD)

        %number to string
        subName = num2str(subNameD(iSubject));
        savName = strcat(subName, '_', curTask, '_', subName, num2str(sr), 'Hz', num2str(lp), 'LP',  erase(num2str(hp),'0.'), 'HP', 'Notch');

        %locate sub folder
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);
        
        % load existing preprocessed dataset
        EEG = pop_loadset('filename',[subName '_' curTask '_EpochRejICARej.set'],'filepath',pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
        
        % load file saved in s1 before removing bad channels
        EEG = pop_loadset('filename',[ savName '.set'], 'filepath', pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
        
        if inter_option
            % interpolate channels on dataset 1, based on dataset 2(dataset from step 1),
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0);
            EEG = eeg_checkset( EEG );
            EEG = pop_interp(EEG, ALLEEG(2).chanlocs, 'spherical');
        end
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',...
            [subName '_' curTask '_EpochRejICARejInter.set'],'overwrite','on','gui','off');
        
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
        end
        

        EEG = pop_artmwppth( EEG , 'Channel',  1:61, 'Flag',  1, 'Threshold',  threshold, 'Twindow', lwindow,...
            'Windowsize',  win_size, 'Windowstep',win_step ); % GUI: 15-Jul-2021 23:25:26

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off');
        
        EEG = pop_saveset( EEG, 'filename', [ subName '_' curTask '_EpochArtDet'],'filepath', pwd);
        
        RejNum(iSubject,1) = subNameD(iSubject);
        RejNum(iSubject,2) = size(find(EEG.reject.rejmanual==1),2);
        
        % Reject detected artefacts
        EEG = eeg_checkset( EEG );
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, find(EEG.reject.rejmanual==1),0);
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5, 'gui', 'off');
        EEG = pop_saveset( EEG, 'filename', [ subName '_' curTask '_EpochArtRej'],'filepath', pwd);
        
        ALLEEG = []; EEG = []; CURRENTSET = [];
   
    end
    
    cd(outputParentFolder);
    xlswrite([curTask '_RejectTrialNumber_thr' num2str(threshold) 'uV' date], RejNum);
    
    disp(['Completed: ' curTask]);
    

end

end
