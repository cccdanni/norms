eeglab;

% for each subject
task = 'preLearningImplicit';
subjectID = '1222';
subjectName = strcat(subjectID, '_', task, '_Epoch.set');
subjectFolder = fullfile ('/home/chendanni/Documents/Norms/analysis/EEGAnalysisResults/Preprocessing', subjectID);
EEG = pop_loadset('filename',subjectName,'filepath',subjectFolder);
eeglab redraw
pop_eegplot( EEG, 1, 1, 1);

% reject bad epoch(s) 

% save new dataset
subjectNameNew = strcat(subjectID, '_', task, '_EpochRej');
EEG = pop_editset(EEG, 'setname', subjectNameNew, 'run', []);
EEG = pop_saveset( EEG, 'filename',strcat(subjectID, '_', task, '_EpochRej'),'filepath',subjectFolder);
eeglab redraw

% save reject epoch
ori = {ALLEEG(1).epoch.eventbepoch};
ori_trans = [];
rmm = {ALLEEG(2).epoch.eventbepoch};
rmm_trans = [];
for i = 1:length(ori)
    tmp = ori{i};
    if iscell(tmp)
        tmp = cell2mat(tmp);
        ori_trans = [ori_trans tmp(1)];
    else
        ori_trans = [ori_trans tmp];
    end
end
for i = 1:length(rmm)
    tmp = rmm{i};
    if iscell(tmp)
        tmp = cell2mat(tmp);
        rmm_trans = [rmm_trans tmp(1)];
    else
        rmm_trans = [rmm_trans tmp];
    end
end
badEpoch = setdiff(ori_trans, rmm_trans);

% badEpoch = [];

save(strcat(subjectFolder, '/', subjectID, '_', task, '_BadEpoch'), 'badEpoch');
dir(subjectFolder)

clear badEpoch ori ori_trans rmm rmm_trans;
ALLEEG = []; EEG = []; ALLCOM = [];
