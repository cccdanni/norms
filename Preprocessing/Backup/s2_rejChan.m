%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Input: *_250Hz05HP30LPNotch.set
% Output: 1) "pop_rejchan":*__badChan.mat; "pop_clean_rawdata":*_badChanCleanChan.mat
% Update Date: July-12-2021
% Last run date: July-12-2021

%% Step 2: Automatically Reject Bad Channels 

clear all; clc;

subNameD = [1201:1242, 1244:1253]; 
workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)

rawFolder = fullfile(workFolder, 'EEGRawData'); 
% task = {'preLearningImplicit', 'Learning', 'postLearningMemory', 'postLearningImplicit'};
task = {'preLearningImplicit'};

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path
channelLocationFile = '/home/chendanni/MATLAB/toolbox/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';% Channel location file(path according to your file location)
% change accordingly to your eeglab version and location

detect_method = "pop_rejchan";
%detect_method = "pop_clean_rawdata";

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task)
    
    curTask = task{iTask};
    
    for iSubject = 1:length(subNameD)
    
    % number to string
    subName = num2str(subNameD(iSubject));
    
    outputSubjectFolder = fullfile(outputParentFolder ,subName);
    cd(outputSubjectFolder);
    
    % Read Data
    EEG = pop_loadset('filename',[ subName '_' curTask '_250Hz05HP30LPNotch.set'],...
        'filepath',[outputSubjectFolder '/']);
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
    
    % Bad channel: mark down bad channels while manually rejecting data. Type
    % in channel names (surrounded by '') before running s2.
    
    
    % Detect Bad Channel
   
    if detect_method == "pop_rejchan"
            % extract bin-based epochs
            switch curTask
                case 'preLearningImplicit'
                    EEG = pop_epoch( EEG, [] , [-0.2           1], 'newname', '1219250Hz05HP30LPNotch epochs', 'epochinfo', 'yes'); % set the epoch length and baseline correction window
                case 'Learning'
                    EEG = pop_epoch( EEG, [] , [-0.3           3], 'newname', '1219250Hz05HP30LPNotch epochs', 'epochinfo', 'yes'); % set the epoch length and baseline correction window
                case 'postLearningMemory'
                    EEG = pop_epoch( EEG, [] , [-0.2           1], 'newname', '1219250Hz05HP30LPNotch epochs', 'epochinfo', 'yes'); % set the epoch length and baseline correction window
                case 'postLearningImplicit'
                    EEG = pop_epoch( EEG, [] , [-0.2           1], 'newname', '1219250Hz05HP30LPNotch epochs', 'epochinfo', 'yes'); % set the epoch length and baseline correction window
            end 
        [EEG, indelec] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan] ,'threshold',5,'norm','on','measure','kurt');
        badChan = indelec;
        save([outputSubjectFolder '/' curTask '_' subName '_badChan.mat'],'badChan');
    elseif detect_method == "pop_clean_rawdata"
            EEG_cleanchan = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,...
                'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
            indelec_original = [EEG.chanlocs.urchan];
            indelec_cleaned = [EEG_cleanchan.chanlocs.urchan];
            indelec = indelec_original(setdiff(indelec_original, indelec_cleaned));
            badChan = indelec;
            save([outputSubjectFolder '/' curTask '_' subName '_badChanCleanChan.mat'],'badChan');
    end
    
end

end

% examine bad channels:

for iTask = 1:length(task)
    
    curTask = task{iTask};
    thisBadChanTask = {};
    
    for iSubject = 1:length(subNameD)
        
        % number to string
        subName = num2str(subNameD(iSubject));
        
        outputSubjectFolder = fullfile(outputParentFolder ,subName);
        cd(outputSubjectFolder);
        
        % Read Data
        
        if detect_method == "pop_rejchan"
            thisbadChanList = load([outputSubjectFolder '/' curTask '_' subName '_badChan.mat']);
            thisBadChanTask{iSubject, 1} = subNameD(iSubject);
            thisBadChanTask{iSubject, 2} = thisbadChanList.badChan;
            thisBadChanTask{iSubject, 3} = length(thisbadChanList.badChan); 
            cd(outputParentFolder);
            % save([curTask,'_badChanRejAll.mat'],'thisBadChanTask');
            T = cell2table (thisBadChanTask, 'VariableNames', {'SubjectID','RemoveChannels','nbchans'});
            writetable(T, [curTask,'_badChanRejAll.csv']);
        elseif detect_method == "pop_clean_rawdata"
            thisbadChanList = load([outputSubjectFolder '/' curTask '_' subName '_badChanCleanChan.mat']);
            thisBadChanTask{iSubject, 1} = subNameD(iSubject);
            thisBadChanTask{iSubject, 2} = thisbadChanList.badChan;
            thisBadChanTask{iSubject, 3} = length(thisbadChanList.badChan); 
            cd(outputParentFolder);
            % save([curTask,'_badChanRejAllCleanChan.mat'],'thisBadChanTask');
            T = cell2table (thisBadChanTask, 'VariableNames', {'SubjectID','RemoveChannels','nbchans'});
            writetable(T, [curTask,'_badChanRejAllCleanChan.csv']);
        end
    
    end
    
end

% NOTE:
% Should further manually check channel and revise the badChanRejAll.mat
% accordingly.
