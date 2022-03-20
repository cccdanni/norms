%% pre-processing EEG data
% Step 8: Detect Incorrect Resposne in Epoch in ERPlab
% Project: Social Norms Learning 
% Author: Danni Chen 
% Update Date: 2022-02-22 
% Input: *_EpochArtRej
% Parameter:
%   1) subNameD: Vector, a list of subject number
%   2) task_list
%   3) workFolder

function s9_rejFalseResp(subNameD, task_list, workFolder)

if nargin == 0

    % participant info set up
    subNameD = [1205:1242,1244:1253]; 
    
    % directory set up
    workFolder = '/home/chendanni/Documents/Norms/analysis/'; % lab server

    % task set up
    task_list = {'postLearningImplicit', 'preLearningImplicit'};
    
    
end

cd (workFolder);

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults','Preprocessing'); % set your output path

% read Behavior Response
preLearningBehaResp = readtable('BehaData/PreLearningImplicitBehaData2021-02-23.csv');
postLearningBehaResp = readtable('BehaData/PostLearningImplicitBehaData2021-02-23.csv');

%load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
rmN = [];

for iTask = 1:length(task_list)
    
    curTask = task_list{iTask};
    
    for iSubject = 1:length(subNameD)
        
        %number to string
        subName = num2str(subNameD(iSubject));
        
        %locate sub folder
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        cd(outputSubjectFolder);
        
        % load existing preprocessed dataset        
        EEG = pop_loadset('filename',[subName '_' curTask '_Epoch.set'],'filepath',pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
        
        % special treatment
        switch curTask
            case 'preLearningImplicit'
                switch subName
                    case '1206'
                        rmN = 610; 
                    case '1238'
                        rmN = 329; 
                end
            case 'postLearningImplicit'
                switch subName
                    case '1211'
                        rmN = 132;
                end
        end
        
        if ~isempty(rmN)
            EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
            EEG = pop_rejepoch( EEG, rmN, 0);
        end
        
        % find incorrect response trials 
        if strcmp(curTask, 'preLearningImplicit')
            thisBeha = preLearningBehaResp(find(preLearningBehaResp.SubjectID == subNameD(iSubject)),:);
        elseif strcmp(curTask, 'postLearningImplicit')
             thisBeha = postLearningBehaResp(find(postLearningBehaResp.SubjectID == subNameD(iSubject)),:);
        end
        thisBeha = thisBeha(:, {'triggerID', 'photoKey_corr', 'BlockRep_thisN','Block_thisN'});
        thisBeha = table2array(thisBeha);
        
        % EEG trials
        
        thisEpoch = EEG.epoch;
        thisEpochTrans = nan (2, length(thisEpoch));
        for i = 1: length(thisEpoch)
            
            this_eventcodelabel = thisEpoch(i).eventcodelabel;
            this_eventbepoch    = thisEpoch(i).eventbepoch;
            
            if isdoublep (this_eventbepoch)
                thisEpochTrans(1, i) = str2double(regexp(this_eventcodelabel{1,1},'\d+$','match'));
                thisEpochTrans(2, i) = this_eventbepoch;
            elseif strcmp(class(this_eventbepoch), 'cell')
                if ~strcmp (this_eventcodelabel{1,1}, 'boundary')
                    thisEpochTrans(1, i) = str2double(regexp(this_eventcodelabel{1,1},'\d+$','match'));
                    thisEpochTrans(2, i) = this_eventbepoch{1,1};
                else 
                    thisEpochTrans(1, i) = str2double(regexp(this_eventcodelabel{1,2},'\d+$','match'));
                    thisEpochTrans(2, i) = this_eventbepoch{1,2};
                end
            end
            
            this_eventcodelabel = []; this_eventbepoch = [];
          
        end
        thisEpochTrans = thisEpochTrans';
        
        
        A = thisEpochTrans(:, 1)';
        Au = unique(A);
        for k = 1:numel(Au)
            [~,ic] = find(A==Au(k));
            [~,ia,aval] = find(A(ic));
            B(ic) = ia;
        end
        thisEpochTrans(:, 3) = B;
        A = []; Au =[]; B = []; ic = []; ia=[]; aval=[]; k = [];
        thisEpochTrans(:, 3) = thisEpochTrans(:, 3)- 1; % block ID 
        thisEpochTrans2 = thisEpochTrans(thisEpochTrans(:,1) <99, :);
        
        for k = 1:size(thisEpochTrans2, 1)-1
            if (thisEpochTrans2(k, 3))>(thisEpochTrans2(k+1,3))
                thisEpochTrans2(k+1, 3) = thisEpochTrans2(k, 3);
            end
        end
        
        
        for k = 1:size(thisEpochTrans2,1)
            thistrigger = thisEpochTrans2(k, 1);
            thisblock   = thisEpochTrans2(k, 3);
            thisEpochTrans2(k, 4) = thisBeha(thisBeha(:,1)==thistrigger & thisBeha(:,3)==thisblock, 2);
        end
        
        EEG = []; ALLEEG = []; CURRENTSET = [];
        
        % load preprocessed 
        EEG = pop_loadset('filename',[subName '_' curTask '_EpochArtRej.set'],'filepath',pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 2 );
        
        thisEpoch = EEG.epoch;
        allEEGTrial = nan (2, length(thisEpoch));
        for i = 1: length(thisEpoch)
            
            this_eventcodelabel = thisEpoch(i).eventcodelabel;
            this_eventbepoch    = thisEpoch(i).eventbepoch;
            
            if isdoublep (this_eventbepoch)
                allEEGTrial(1, i) = str2double(regexp(this_eventcodelabel{1,1},'\d+$','match'));
                allEEGTrial(2, i) = this_eventbepoch;
            elseif strcmp(class(this_eventbepoch), 'cell')
                if ~strcmp (this_eventcodelabel{1,1}, 'boundary')
                    allEEGTrial(1, i) = str2double(regexp(this_eventcodelabel{1,1},'\d+$','match'));
                    allEEGTrial(2, i) = this_eventbepoch{1,1};
                else 
                    allEEGTrial(1, i) = str2double(regexp(this_eventcodelabel{1,2},'\d+$','match'));
                    allEEGTrial(2, i) = this_eventbepoch{1,2};
                end
            end
            
            this_eventcodelabel = []; this_eventbepoch = [];
          
        end
        allEEGTrial = allEEGTrial';
        
        
        % link all trials need removal
        objectTrials = allEEGTrial(allEEGTrial(:,1)==99, 2)';
        wrongTrials  = thisEpochTrans2(thisEpochTrans2(:,4)==0, 2)';
        removeTrials = [objectTrials wrongTrials];
           
        if ~isempty(rmN)
            removeTrials = [removeTrials rmN];
        end
        
        
        [tf, indec] = ismember(removeTrials, allEEGTrial(:,2));

        % Reject incorrect response trials and object viewing trials
        EEG = eeg_checkset( EEG );
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, indec,0);
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'savenew',...
            [outputSubjectFolder filesep subName '_' curTask '_EpochArtRej_FalseRespRej'],'gui','off');
        
        % save record
        thisEpochTrans2 = thisEpochTrans2( ismember(thisEpochTrans2(:,2), allEEGTrial(:,2)), :);
        save (strcat(subName, '_', curTask, '_EpochInfo.mat'), 'thisEpochTrans2');
        
        % save record general
        RejNum(iSubject,1) = subNameD(iSubject); % subject ID
        RejNum(iSubject,2) = size(EEG.epoch,2);  % remained trials
        RejNum(iSubject,3) = sum(strcmp({EEG.event.type},'B9(s99)')); % remained object trials
        RejNum(iSubject,4) = sum(~strcmp({EEG.event.type},'B9(s99)'));% remained non-object trials
        RejNum(iSubject,5) = sum(~strcmp({EEG.event.type},'B9(s99)')) > 240; % if accept the subject
        RejNum(iSubject,6) = 60 - sum([EEG.event.bini]==1); % removed trials bin = 1
        RejNum(iSubject,7) = 60 - sum([EEG.event.bini]==2); % removed trials bin = 2
        RejNum(iSubject,8) = 60 - sum([EEG.event.bini]==3); % removed trials bin = 3
        RejNum(iSubject,9) = 60 - sum([EEG.event.bini]==4); % removed trials bin = 4
        RejNum(iSubject,10) = 60 - sum([EEG.event.bini]==5); % removed trials bin = 5
        RejNum(iSubject,11) = 60 - sum([EEG.event.bini]==6); % removed trials bin = 6
        RejNum(iSubject,12) = 60 - sum([EEG.event.bini]==7); % removed trials bin = 7
        RejNum(iSubject,13) = 60 - sum([EEG.event.bini]==8); % removed trials bin = 8
        
        
        EEG = []; ALLEEG = []; ALLCOM = []; rmN =[];
        thisEpoch = []; thisEpochTrans = []; thisBeha = []; thisEpochTrans = []; 
        thisEpochTrans2 = []; objectTrials = []; wrongTrials = []; removeTrials = []; 
        allEEGTrials = []; tf = []; indec = [];
   
    end
    
    cd(outputParentFolder);
    T = array2table (RejNum, ...
        'VariableNames', {'SubjectID', 'Trials', 'Object_Trials', 'Nonobject_Trials',...
        'Accept_TF', 'Removed_bin1', 'Removed_bin2', 'Removed_bin3', 'Removed_bin4',...
        'Removed_bin5', 'Removed_bin6', 'Removed_bin7', 'Removed_bin8'})
    writetable(T, [curTask 'RejectTrialNumber_FalseResp' date, '.xls']);
    
    disp(['Completed: ' curTask]);

end

end