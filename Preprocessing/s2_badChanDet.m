%% pre-processing EEG data
% s2: automatically detect bad channels
% Project: Social Norms Learning
% Author: Danni Chen 
% Update Date: 2022-02-16 
% Input: dataset from previous step
% Output: 1) "pop_rejchan":*__badChan.mat; "pop_clean_rawdata":*_badChanCleanChan.mat
% NOTE: Should further manually check channel accordingly.
% Parameter:
%   1) subNameD: Vector, a list of subject number
%   2) detect_method: a) pop_rejchan b) pop_clean_rawdata
%   3) sr: sampling rate
%   4) lp: low-pass filter
%   5) hp: high-pass filter
%   6) workFolder
% History
%   1) 2022-02-16: load dataset after removing bad section


function s2_badChanDet(subNameD, task_list, detect_method, sr, lp, hp, workFolder)

if nargin == 0

    % participant info set up
    subNameD = [1201:1242, 1244:1253];
    
    % directory set up
    workFolder = '/home/chendanni/Documents/Norms/analysis/';

    % task set up
    task_list = {'preLearningImplicit','postLearningImplicit','Learning'};
    
    % parameter set up
    sr = 250;
    lp = 30;
    hp = 0.05;
    notch = 1;
    linefreq = 50;

    % detect_method
    detect_method = "pop_rejchan";
    
end

outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing'); % set your output path
    
% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


%% detect bad channel
for iTask = 1:length(task_list)

    curTask = task_list{iTask};

    for iSubject = 1:length(subNameD)
    
        % number to string
        subName = num2str(subNameD(iSubject));
        savName = strcat(subName, '_', curTask, '_', subName, num2str(sr), 'Hz', num2str(lp), 'LP',  erase(num2str(hp),'0.'), 'HP', 'NotchBadSec');
    
        outputSubjectFolder = fullfile(outputParentFolder ,subName);
        cd(outputSubjectFolder);

        % Read Data
        EEG = pop_loadset('filename',strcat(savName, '.set'), 'filepath',pwd);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );

        % Detect Bad Channel

        if detect_method == "pop_rejchan"
            
            [EEG, indelec] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan] ,'threshold',5,'norm','on','measure','kurt');
            badChan = indelec;
            save([outputSubjectFolder '/' curTask '_' subName '_badChanBadSec.mat'],'badChan');
        
        elseif detect_method == "pop_clean_rawdata"
            
            EEG_cleanchan = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,...
                'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
            indelec_original = [EEG.chanlocs.urchan];
            indelec_cleaned = [EEG_cleanchan.chanlocs.urchan];
            indelec = indelec_original(setdiff(indelec_original, indelec_cleaned));
            badChan = indelec;
            save([outputSubjectFolder '/' curTask '_' subName '_badChanCleanChanBadSec.mat'],'badChan');
            
        end

    end

end

%% examine bad channel

for iTask = 1:length(task_list)

    curTask = task_list{iTask};
    thisBadChanTask = {};

    for iSubject = 1:length(subNameD)

        % number to string
        subName = num2str(subNameD(iSubject));

        outputSubjectFolder = fullfile(outputParentFolder ,subName);
        cd(outputSubjectFolder);

        % Read Data

        if detect_method == "pop_rejchan"
            
            thisbadChanList = load([outputSubjectFolder '/' curTask '_' subName '_badChanBadSec.mat']);
            thisBadChanTask{iSubject, 1} = subNameD(iSubject);
            thisBadChanTask{iSubject, 2} = thisbadChanList.badChan;
            thisBadChanTask{iSubject, 3} = length(thisbadChanList.badChan); 
            cd(outputParentFolder);
            T = cell2table (thisBadChanTask, 'VariableNames', {'SubjectID','RemoveChannels','nbchans'});
            writetable(T, [curTask,'_badChanRejAllBadSec.csv']);
            
        elseif detect_method == "pop_clean_rawdata"
            
            thisbadChanList = load([outputSubjectFolder '/' curTask '_' subName '_badChanCleanChanBadSec.mat']);
            thisBadChanTask{iSubject, 1} = subNameD(iSubject);
            thisBadChanTask{iSubject, 2} = thisbadChanList.badChan;
            thisBadChanTask{iSubject, 3} = length(thisbadChanList.badChan); 
            cd(outputParentFolder);
            T = cell2table (thisBadChanTask, 'VariableNames', {'SubjectID','RemoveChannels','nbchans'});
            writetable(T, [curTask,'_badChanRejAllCleanChanBadSec.csv']);
            
        end

    end

end


end
