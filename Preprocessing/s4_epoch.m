%% pre-processing EEG data
% s4: epoch in ERPLab
% Project: Social Norms Learning
% Author: Danni Chen 
% Update Date: 2022-02-16 
% Input: *_InterReRef.set (from s3); *_binassigned.txt
% Output: *_Epoch.set
% Parameter:
%   1) subNameD: Vector, a list of subject number
%   2) task_list
%   3) workFolder
%   4) binFolder
% Note that: You need to preprare your binlist before running, and store them under 
% binlistFolder.

function s4_epoch(subNameD, task_list, workFolder)

if nargin == 0

    % participant info set up
    subNameD = [1201:1242, 1244:1253];
    
    % directory set up
    workFolder = '/Users/danni/Desktop/Norms/'; % in MAC

    % task set up
    task_list = {'preLearningImplicit','postLearningImplicit','Learning'};
    
end

outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path
binlistFolder = fullfile(workFolder, 'Binlist'); % path to your binlist
binassname = '_binassigned.txt';

%load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task_list)
    
    curTask = task_list{iTask};
    
    epoch_checking = struct;

    for iSubject = 1:length(subNameD)

        %number to string
        subName = num2str(subNameD(iSubject));

        %locate sub folder
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);

        % load existing preprocessed dataset
        EEG = pop_loadset('filename',[subName '_' curTask '_InterReRef.set'],'filepath',pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

        % Trigger/Binlist file name, if you have counterbalancing, set accordingly
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
            case 'postLearningImplicit'
                binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '_implicit.txt'];
                for i = 1:(size(EEG.urevent,2)-624-1)
                   EEG.urevent(i).type = 'practice';
                   EEG.event(i).codelabel = 'practice';
                   EEG.event(i).type = 199;
                end
        end

        % bin name
        binfile = fullfile (binlistFolder,binlistfile);
        outputeventfile = fullfile (outputSubjectFolder, strcat(subName, '_', curTask, binassname));
        inputeventfile  = fullfile (outputSubjectFolder, strcat(subName, '_', curTask, '_eventlist.txt'));
        
        % create eventlist
        EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, ...
            'BoundaryString', { 'boundary' }, 'Eventlist', inputeventfile); 
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName 'InterReRef.set'],'gui','off');

        % Assign Binlister    
        EEG  = pop_binlister( EEG , 'BDF', binfile,... 
            'ExportEL', outputeventfile,... 
            'ImportEL', inputeventfile,... 
            'IndexEL',  1, 'SendEL2', 'All', 'Voutput', 'EEG' ); 

        % extract bin-based epochs
        % set the epoch length and baseline correction window
        switch curTask
            case 'preLearningImplicit'
                EEG = pop_epochbin( EEG , [-1000.0  2000.0],  'none'); 
            case 'Learning'
                EEG = pop_epochbin( EEG , [-1000.0  4000.0],  'none'); 
            case 'postLearningImplicit'
                EEG = pop_epochbin( EEG , [-1000.0  2000.0],  'none'); 
        end

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'savenew',...
            fullfile (outputSubjectFolder, strcat(subName, '_', curTask, '_Epoch.set')), 'gui','off');

        % epoch checking history
        epoch_checking(iSubject).subjectID = subName;
        epoch_checking(iSubject).task = curTask;
        epoch_checking(iSubject).trialsperbin = EEG.EVENTLIST.trialsperbin;
    end
    
    save (fullfile(outputParentFolder, strcat(curTask, '_epoch_checking.mat')), 'epoch_checking');

end

end