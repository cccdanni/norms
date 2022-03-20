%% pre-processing EEG data
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and BrainVision .bv Data Loading Plugin
% Author: Danni Chen 
% Update Date: 2022-01-25

%% basic information set up

clear;clc;

% participant info set up
subNameD = [1205:1242, 1244:1253];
%subNameD = fliplr(subNameD);

% directory set up
workFolder = '/home/chendanni/Documents/Norms/analysis/'; % lab server
%workFolder = '/Users/danni/Desktop/Norms'; % my MAC 
rawFolder  = fullfile (workFolder, 'EEGRawData');
outputParentFolder = fullfile (workFolder, 'EEGAnalysisResults', 'Preprocessing');
scriptFolder = fullfile (workFolder, 'MyScripts');
addpath (genpath(scriptFolder));

% channel location file
channelLocationFile = '/home/chendanni/MATLAB/toolbox/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';% Channel location file(path according to your file location)

% task set up
task_list = {'preLearningImplicit','postLearningImplicit','Learning'};
% task_list = {'Learning'};

% parameter set up
sr = 250;
lp = 30;
hp = 0.05;
notch = 1;
linefreq = 50;
inter_option = true;
threshold = 100;
win_size = 200;
win_step = 100;

%% s1: De-sampling, Filtering, Remove Channels
% s1_fil (subNameD, task_list, sr, hp, lp, notch, linefreq,  workFolder);
% before running next step, load dataset saved from step 1, 
% and then mark down bad sections and bad channels. 
% 'detectBadSection' is of particular importance here. 

%% s2: detect channel (optional)
% s2_badChanDet(subNameD, task_list, 'pop_rejchan', sr, lp, hp, workFolder); 
% s2_badChanDet(subNameD, task_list, 'pop_clean_rawdata', sr, lp, hp, workFolder);
% The automatic bad channel detection is unsatisfactory for my case. Too
% many channels get rejected. 

%% s3: remove big artifact and rereference
% remove bad sections (optional, if you have remove in the s2), 
% remove bad channels, interpolate, reref, remove bad channels again.
% re-reference is optional here, you could re-reference after remove
% artifacts
% s3_reref(subNameD, task_list, 'manual', sr, lp, hp, workFolder);

%% s4: epoch (longer than needed)
% s4_epoch(subNameD, task_list, workFolder);
% Before Run s5 ICA, Load each exported data and reject trials with big artfacts
% 'detectBadEpoch' is of particular importance here. 

%% s5: manually detect and remove bad epochs
% you should have '_EpochRej.set' & '_BadEpoch.mat' after this step

%% s6: ICA 
% s6_ica(subNameD, task_list, workFolder, sr, lp, hp, notch, linefreq, channelLocationFile, 1);
%  Then, manually check ICA components to remove EOG or other artefacts. 
%  Save the data sets as "...rej.set". 
%  Then run S7 to detect and reject component in epoch.
s6_ica(1222, task_list, workFolder, sr, lp, hp, notch, linefreq, channelLocationFile, 1);

%% s7: remove components (with reference of ICLabel)
% s7_rejComp(subNameD, task_list, workFolder);
s7_rejComp(1222, task_list, workFolder);

%% s8: reject epoch with remaining big artifacts and interpolate
s8_rejArt(subNameD, task_list, workFolder, inter_option, threshold, sr, lp, hp, win_size, win_step);

%% s9: automatically remove bad epochs and unused epochs (optional)
s9_rejFalseResp(subNameD, {'preLearningImplicit', 'postLearningImplicit'}, workFolder)
% unused epochs in my case - incorrect response for pre- and post implicit
% perception tasks

%% the data have been cleaned! Now we could move on to the fancy stuff..
%% Good Luck!