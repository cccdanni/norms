%% ReadME
% Preprocessing EEG Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: July-12-2021
% Last run date: July-13-2021   
% Input: raw BrainVision data
% Output: *_250Hz05HP30LPNotch.set

% Mostly revised based on EEG Preprocessing Scripts of Hu Lab
% Dependency: EEGlab, ERPlab, and ANT .cnt Data Loading Plugin
% Authors: Yao Ziqing, Lin Xuanyi, Hu Xiaoqing
% Create Date: 2018.07
% Update Date: 2020.07

%% s1: De-sampling, Filtering, Remove Channels

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

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iTask = 1:length(task)
    curTask = task{iTask};
    
    for iSubject = 1:length(subNameD)
        
        % number to string
        subName = num2str(subNameD(iSubject));
        
        % creat folder for each subject
        if ~exist(fullfile(outputParentFolder, subName))
            mkdir(outputParentFolder,subName)
        end
        
        outputSubjectFolder = fullfile(outputParentFolder,subName);
        
        % load data from the rawFolder: for ANT.cnt: loadeep_v4; for other file
        % type, manually load data and check the function by 'eegh'
        curRawFolder = fullfile(rawFolder, ['sub', subName]);
        cd (curRawFolder);
        curRawFile = dir(['SN_Study1_' curTask '_sub', subName,'*.vhdr']);
        curRawFile = curRawFile.name;
        
        cd(outputSubjectFolder);
        
        EEG = pop_loadbv(curRawFolder, curRawFile); % pop_loadbv for .vhdr 
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname', [subName],'gui','off');
        
        % down sample to 250Hz: optional
        EEG = pop_resample( EEG, 250); 
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250Hz'],'overwrite','on','gui','off'); %Not saved
        
        %low pass 30 for ERP
        %EEG = pop_eegfiltnew(EEG, [], 30, 110, 0, [], 1);
        EEG = pop_eegfiltnew(EEG, 'hicutoff',30)
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250Hz30LP'],'overwrite','on','gui','off'); 
        
        %high-pass filter 
        %EEG = pop_eegfiltnew(EEG, [], 0.1, 33000, 1, [], 1);
        %Changed to eeglab 19 version 
        EEG = pop_eegfiltnew(EEG, 'locutoff',0.05)
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250Hz30LP05HP'],'overwrite','on','gui','off'); 
        
        %Notch filter 50Hz(for HK); filter function from ERPlab
        %CleanLine instead - DC
        %EEG  = pop_basicfilter( EEG,  [1:EEG.nbchan] , 'Boundary', 'boundary', 'Cutoff',  50, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180, 'RemoveDC', 'on' ); 
        EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:63] ,'computepower',1,'linefreqs',50,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250Hz30LP005HPNotch'],'overwrite','on','gui','off'); 
        
        % load channel location file
        EEG = eeg_checkset( EEG );
        EEG = pop_chanedit(EEG, 'lookup',channelLocationFile);
        
        % exclude channels not used
        EEG = eeg_checkset( EEG );
        
        if EEG.nbchan == 88
            EEG = pop_select( EEG,'nochannel',{EEG.chanlocs(65:88).labels 'M1' 'M2' 'EOG'  });
        elseif EEG.nbchan == 66
            EEG = pop_select( EEG,'nochannel',{'M1' 'M2' 'EOG' 'VEOG' 'HEOG'}); 
        elseif EEG.nbchan == 64
            EEG = pop_select( EEG,'nochannel',{'M1' 'M2' 'EOG'});
        elseif EEG.nbchan == 63
            EEG = pop_select( EEG, 'nochannel', {'M1' 'M2'});
        end
        
        EEG = eeg_checkset( EEG );
        
        % Change saved file names if filter parameters are difference
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,...
            'setname',[subName '250Hz05HP30LPNotch'],...
            'savenew',[subName '_' curTask '_250Hz05HP30LPNotch'],'overwrite','on','gui','off'); 
    
    end
end

eeglab redraw

%% Next step(s2), before running s2 script, load dataset saved from this step,
%  and then "Reject continuous data by eye". Reject periods of rest and
%  obvious artefacts. Mark down bad channels.

