% Project: Social Norms Learning (Time resolved)
% Function: [Prelearning] RSA analysis (Time-resolved) 
% Author: Danni Chen 
% Update Date: Apr-29-2021

clear;clc;

project_folder='/home/chendanni/Documents/Norms/data/';
data_folder = fullfile(project_folder, 'EEGAnalysisResults_PreprocessingVersion1');
toolbox_folder = '/home/chendanni/MATLAB/toolbox';
cd(project_folder);

%% add toolboxes
addpath (fullfile(project_folder,'additional_functions')); 
addpath (fullfile(toolbox_folder,'fieldtrip-20210330'));

%% RSA for Prelearning
path_in=fullfile(project_folder,'EEGAnalysisResults_PreprocessingVersion1');
path_out = fullfile(project_folder,'RSA','data',strcat('all_trials_pre_sr',num2str(sr),'window','200','slide','10'));
if ~exist(path_out)
    mkdir(path_out)
end

% define the participants - 
subs = [1204:1242, 1244:1253];
badsubs = [1228, 1237, 1239, 1229];
subs = setdiff(subs, badsubs);
all_subs = {};
for i = 1:length(subs)
    all_subs (i) =  { num2str(subs(i)) } ;
end

cd(path_in);

stim_nr = 80;

% define new Sampling rate
sr=100; % in Hz
step=win/(1/sr)+1;

channels='all';

%define start and end of item window
t_start= 0;
t_end=0.996;
tois=t_start:1/sr:t_end;
t1=t_start:slide:(t_end-win);
t2=t1+win;
ind_t1=1:slide/(1/sr):((numel(tois)-win/(1/sr)));
ind_t2=ind_t1+win/(1/sr);
n_bins=numel(t1);

for bin_nr = 1:n_bins
    % create rdm
    allsubject_neuralrdm           = nan (stim_nr, stim_nr, numel(subs));
    allsubject_behardm             = nan (stim_nr, stim_nr, numel(subs));
    allsubject_coceprdm_conflict   = nan (stim_nr, stim_nr, numel(subs));
    allsubject_coceprdm_inconflict = nan (80, 80, numel(subs));
    
    for sub_nr = 1:length(subs)
        
    end
    
end