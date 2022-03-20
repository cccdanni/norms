% for each subject

eeglab
task = 'Learning';
subjectID = '1237';
subjectName = strcat(subjectID, '_', task, '_', subjectID, '250Hz30LP05HPNotch.set');
subjectFolder = fullfile ('/Users/danni/Desktop/Norms/EEGAnalysisResults/Preprocessing', subjectID);
EEG = pop_loadset('filename',subjectName,'filepath',subjectFolder);
eeglab redraw
pop_eegplot( EEG, 1, 1, 1);

% reject bad section 

% save new dataset
subjectNameNew = strcat(subjectID, '250Hz30LP05HPNotchBadSec');
EEG = pop_editset(EEG, 'setname', subjectNameNew, 'run', []);
EEG = pop_saveset( EEG, 'filename',strcat(subjectID, '_', task, '_', subjectNameNew),'filepath',subjectFolder);

% save reject section
save(strcat(subjectFolder, '/', subjectID, '_', task, '_BadSec'), 'TMPREJ');
clear TMPREJ;
dir(subjectFolder)

% check bad channel
