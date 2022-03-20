% decoding step 1: transform the data
% edited from below source: 

% binorgEEG script:
% Use *right after* Extracting Bin-Based Epochs
% Creates bin-organized .mat files for Decoding
% Writes .mat files to current working directory
% DC: binwise_data: 16 Orientations; Every box: Chans X Timepoints X trials
% DC: n_trials: trials nr in every box

% Place continuous BE datasets (.set & .fdt) in current working directory
% Examples of these kind of data can be found on box: 
% https://ucdavis.box.com/s/jrd9zus448nkc0wjkjf1x5nxbc7oxwg8

% by Aaron M. Simmons, University of California, Davis

clear all;
close all; 


parentfolder = '/home/chendanni/Documents/Norms/analysis/';  
savefolder = fullfile (parentfolder, 'EEGAnalysisResults', 'Decoding', 'Transdata');
if ~exist (savefolder)
    mkdir (savefolder);
end
cd (savefolder) ;


subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238]; % 1201 - 1204: HK Local; 1214, 1215 & 1238 - visually detected bad subjects
subNameD = subNameD(~ismember(subNameD, badsubNameD));
subject_list = num2cell(subNameD); %DC: Has been epoched
numsubjects = length(subject_list);

task = 'preLearningImplicit';

for s = 1:numsubjects
    subject = subject_list{s};
    subjectfolder = fullfile(parentfolder, 'EEGAnalysisResults', 'Preprocessing', num2str(subject_list{s})); %loc of file
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, num2str(subject));
    
    %Initialize EEG
    eeglab;
    
    %load data
    EEG = pop_loadset('filename', strcat ( num2str(subject),'_', task, '_EpochArtRej_FalseRespRej.set'), 'filepath', subjectfolder); 
    %DC: has been assigned bin
    
    %use binorg-EEG function (this function is only on ERPlab V8.01)
    binepEEG_to_binorgEEG(EEG, ['Decoding_BE_', task, '_', num2str(subject)]); 
    % Produces bin-organized .mat file and outputs file into current working directory
    % it will be called "Decoding_BE_XXX" as specified in string filename
   
    close all;
    clear EEG;
    
end

