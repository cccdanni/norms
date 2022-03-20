%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Input: *_InterReRef.set (remove bad channel); *_binassigned.txt
% Output: *_Epoch.set
% Update Date: July-13-2021
% Last run date: July-13-2021

%% Step 5: Reject Bad Epoch

clear;clc;

parentFolder = "/home/chendanni/Documents/Norms/analysis";
pastFolder = fullfile(parentFolder, 'EEGAnalysisResults_PreprocessingVersion1');
currFolder = fullfile(parentFolder, 'EEGAnalysisResults', 'Preprocessing');

subNameD = [1201:1242,1244:1253]; 

% task = {'preLearningImplicit', 'postLearningImplicit'};
task = {'preLearningImplicit','postLearningImplicit'};


%% Extract previous record 
% cd (pastFolder);
% 
% for i = 1:length (task)
%     
%     curTask = task{i};
% 
%     for sub_nr = 1:length(subNameD)
%         thisFolder  = fullfile(pastFolder, num2str(subNameD(sub_nr)));
%         cd (thisFolder);
%         subName = num2str(subNameD(sub_nr));
%     
%         thisEpoch = pop_loadset('filename',[ subName '_' curTask '_Epoch.set']);
%         thisEpochRej = pop_loadset('filename',[ subName '_' curTask '_EpochRej.set']);
%         thisEpochEvent = thisEpoch.event; 
%         thisEpochRejEvent = thisEpochRej.event;
%         thisEpochEvent = [thisEpochEvent.bepoch];
%         thisEpochRejEvent = [thisEpochRejEvent.bepoch];
%         thisRejectedEpoch = thisEpochEvent(find(~ismember(thisEpochEvent, thisEpochRejEvent)))
%         
%         thisOutputFolder = fullfile(currFolder, num2str(subNameD(sub_nr)));
%         cd (thisOutputFolder);
%         save(strcat(subName, '_', curTask, '_rejectedEpoch.mat'), 'thisRejectedEpoch');
%         
%         clear thisEpoch thisEpochRej thisEpochEvent thisEpochRejEvent thisRejectedEpoch thisFolder thisOutputFolder;
%     
%     end
% 
% end

%% Remove bad epochs 
cd (currFolder);

for i = 1:length (task)
    
    curTask = task{i};

    for sub_nr = 1:length(subNameD)
        
        thisFolder  = fullfile(currFolder, num2str(subNameD(sub_nr)));
        cd (thisFolder);
        subName = num2str(subNameD(sub_nr));
       
        EEG = pop_loadset('filename',[ subName '_' curTask '_Epoch.set']);
        load ([ subName '_' curTask '_rejectedEpoch.mat']);
        EEG = pop_selectevent( EEG, 'omitevent',thisRejectedEpoch ,'deleteevents','off','deleteepochs','on','invertepochs','off');
        pop_saveset( EEG, [ subName '_' curTask '_EpochRej.set'], pwd);
        
    end
end