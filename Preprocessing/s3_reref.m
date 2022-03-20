%% pre-processing EEG data
% s3: remove bad channels
% Project: Social Norms Learning
% Author: Danni Chen 
% Update Date: 2022-02-16 
% Input: dataset from previous step
% Function: 
%   1) Remove bad channels saved in s2;
%   2) Then implement interpolation before re-referece to common average;
%   3) re-reference
%   4) Remove interpolated channels before epoch & ICA.
% Parameter:
%   1) subNameD: Vector, a list of subject number
%   2) task_list
%   2) detect_method: a) pop_rejchan b) pop_clean_rawdata c)manual
%   3) sr: sampling rate
%   4) lp: low-pass filter
%   5) hp: high-pass filter
%   6) workFolder

function s3_reref(subNameD, task_list, detect_method, sr, lp, hp, workFolder)

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

    % detect_method
    detect_method = "manual";
    
end

outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path

% load EEGlab

for iTask = 1:length(task_list)
    
    curTask = task_list{iTask};
      
    %% load bad channel table (optional)
    if strcmp (detect_method, "manual")
        
        % Set up the Import Options and import the data
        opts = spreadsheetImportOptions("NumVariables", 5);

        % Specify sheet and range
        opts.Sheet = "Sheet1";
        opts.DataRange = "A2:E157";

        % Specify column names and types
        opts.VariableNames = ["SubjectID", "Task", "badChanNo", "badChanLabel", "Note"];
        opts.VariableTypes = ["double", "categorical", "categorical", "categorical", "categorical"];

        % Specify variable properties
        opts = setvaropts(opts, ["Task", "badChanNo", "badChanLabel", "Note"], "EmptyFieldRule", "auto");

        % Import the data
        badChanTable = readtable(fullfile (outputParentFolder, "badChanManually.xlsx"), opts, "UseExcel", false);

        % Clear temporary variables
        clear opts
        
        % for certain task
        badChanTable = badChanTable (find (badChanTable.Task == curTask),:); 
    end
    
    %% process EEG data
    for iSubject = 1:length(subNameD)
        
        % clear EEG
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; cnt = 1; badChan = [];
        
        % number to string
        subName = num2str(subNameD(iSubject));
        savName = strcat(subName, '_', curTask, '_', subName, num2str(sr), 'Hz', num2str(lp), 'LP',  erase(num2str(hp),'0.'), 'HP', 'NotchBadSec');
        setName = strcat(subName, num2str(sr), 'Hz', num2str(lp), 'LP',  erase(num2str(hp),'0.'), 'HP', 'NotchBadSec');
        
        % locate sub folder
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);

        % load file saved in s1 after removing bad section
        EEG = pop_loadset('filename',strcat(savName, '.set'),'filepath',pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, cnt );    
        cnt = cnt + 1;
        
        % load badSection (must have this file in the same subject folder)
        % load([ subName 'badSection.mat']);
        % badSection = TMPREJ(:,1:2);        
        % OR you could load dataset with bad section already been rejected
        
        % remove bad section
        % EEG = eeg_eegrej( EEG, badSection);
        % [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, cnt,'setname',...
        %     [ subName '250Hz01BP30BadSec']...
        %     ,'overwrite','on','gui','off');
        % cnt = cnt + 1;
        % EEG = eeg_checkset( EEG );
        
        % load bad channel
        if strcmp (detect_method, "manual")
            badChan = badChanTable (find (badChanTable.SubjectID == subNameD(iSubject)),'badChanNo'); 
            badChan = str2num(char(table2array(badChan)));
        elseif strcmp (detect_method, "pop_rejchan")
            load([curTask '_' subName '_badChanBadSec.mat']);
        elseif strcmp (detect_method, "pop_clean_rawdata")
            load([curTask '_' subName '_badChanCleanChanBadSec.mat']);
        else
            error ('no such function');
        end
        
        % deal with bad channel  
        if ~isnan (badChan)
            
            % remove bad channel
            EEG = pop_select(EEG, 'nochannel', badChan);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, cnt,'setname',...
                strcat (setName, 'BadChan'),'overwrite','on','gui','off');
            cnt_rm = cnt;
            cnt = cnt + 1;
            
            % interpolate channels on dataset removed bad channel, based on dataset 1 (dataset from step 1),
            EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, cnt,'setname',...
                strcat(subName, '_Inter.set'),'overwrite','on','gui','off');
            
            % re-referece to common average
            EEG = pop_reref( EEG, []);
            
            % remove bad channel again, then save set before next step
            EEG = pop_select(EEG, 'nochannel', badChan);
            
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname',...
                strcat(subName, '_', curTask, '_InterReRef.set'),...
                'savenew', strcat(subName, '_', curTask, '_InterReRef.set'),'gui','off'); 
            
        else
            
            % re-referece to common average
            EEG = pop_reref( EEG, []);
            
            % save set before next step
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname',...
                strcat(subName, '_', curTask, '_InterReRef.set'),...
                'savenew', strcat(subName, '_', curTask, '_InterReRef.set'),'gui','off'); 
        end

    end
end

end