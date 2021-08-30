%% Step 8: Detect Incorrect Resposne in Epoch in ERPlab

clear; clc;

subNameD = [1205:1242,1244:1253]; 
workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
%task = {'preLearningImplicit', 'Learning', 'postLearningMemory', 'postLearningImplicit'};
task = {'preLearningImplicit','postLearningImplicit'};

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults','Preprocessing'); % set your output path

% read Behavior Response
preLearningBehaResp = readtable('BehaData/PreLearningImplicitBehaData2021-02-23.csv');
postLearningBehaResp = readtable('BehaData/PostLearningImplicitBehaData2021-02-23.csv');

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
        EEG = pop_loadset('filename',[subName '_' curTask '_EpochArtRej.set'],'filepath',outputSubjectFolder);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
        
        % find incorrect response trials 
        if strcmp(curTask, 'preLearningImplicit')
            thisBeha = preLearningBehaResp(find(preLearningBehaResp.SubjectID == subNameD(iSubject)),:);
        elseif strcmp(curTask, 'postLearningImplicit')
             thisBeha = postLearningBehaResp(find(postLearningBehaResp.SubjectID == subNameD(iSubject)),:);
        end
        
        allEEGTrial = nan(1, length(EEG.epoch));
        
        for i = 1: length(EEG.epoch)  
            thisEpoch = EEG.epoch(i).eventbepoch;
            
            if length(thisEpoch)>1   
                thisEpoch = thisEpoch(1); 
            end
            
            if strcmp(class(thisEpoch), 'cell')
                allEEGTrial(i) = cell2mat(thisEpoch);
            elseif strcmp(class(thisEpoch), 'double')
                allEEGTrial(i) = thisEpoch;
            end
            
        end
        
        thisWrongTrial = find(thisBeha.photoKey_corr == 0);
        
        [tf, indec] = ismember(thisWrongTrial, allEEGTrial);
        
        % Reject incorrect response trials
        EEG = eeg_checkset( EEG );
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, indec,0);
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'savenew',...
            [outputSubjectFolder filesep subName '_' curTask '_EpochArtRej_FalseRespRej'],'gui','off');
        
        RejNum(iSubject,1) = subNameD(iSubject);
        RejNum(iSubject,2) = size(EEG.event,2);
        RejNum(iSubject,3) = sum(strcmp({EEG.event.type},'B9(s99)'));
        RejNum(iSubject,4) = sum(~strcmp({EEG.event.type},'B9(s99)'));
        RejNum(iSubject,5) = sum(~strcmp({EEG.event.type},'B9(s99)')) > 240;
   
    end
    
    cd(outputParentFolder);
    xlswrite([curTask 'RejectTrialNumber_FalseResp' date], RejNum);
    
    disp(['Completed: ' curTask]);

end

%% DONE!