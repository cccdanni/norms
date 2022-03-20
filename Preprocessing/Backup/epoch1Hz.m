function EEG = epoch1Hz(subjectID, curTask)
% Function: apply same parameters except high-passed 1 Hz
% Arg:
%   Input: subjectID, curTask
%   Output: EEG

subNameD = subjectID; 
workFolder = '/home/chendanni/Documents/Norms/analysis/';
outputFolder = fullfile(workFolder, 'EEGAnalysisResults', 'Preprocessing');
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
channelLocationFile = '/home/chendanni/MATLAB/toolbox/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';% Channel location file(path according to your file location)

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% s1: De-sampling, Filtering, Remove Channels

% number to string
subName = num2str(subNameD);

% load data from the rawFolder: for ANT.cnt: loadeep_v4; for other file
% type, manually load data and check the function by 'eegh'
curRawFolder = fullfile(rawFolder, ['sub', subName]);
cd (curRawFolder);
curRawFile = dir(['SN_Study1_' curTask '_sub', subName,'*.vhdr']);
curRawFile = curRawFile.name;
EEG = pop_loadbv(curRawFolder, curRawFile); % pop_loadbv for .vhdr 

outputSubjectFolder = fullfile(outputFolder, subName);
cd(outputSubjectFolder);

% down sample to 250Hz: optional
EEG = pop_resample( EEG, 250); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250Hz'],'overwrite','on','gui','off'); %Not saved
        
%low pass 30 for ERP
EEG = pop_eegfiltnew(EEG, 'hicutoff',30)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250HzLP30'],'overwrite','on','gui','off'); %Not saved
        
%high-pass filter 
EEG = pop_eegfiltnew(EEG, 'locutoff',1) % For create a 1 Hz high-passed date
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250HzLP30HP1'],'overwrite','on','gui','off'); %Not saved
        
%Notch filter 50Hz(for HK); filter function from ERPlab
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:63] ,'computepower',1,'linefreqs',50,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[subName '250HzLP30HP1Notch50'],'overwrite','on','gui','off'); %Not saved

% load channel location file
EEG = pop_chanedit(EEG, 'lookup',channelLocationFile);
        
% exclude channel
if EEG.nbchan == 88
    EEG = pop_select( EEG,'nochannel',{EEG.chanlocs(65:88).labels 'M1' 'M2' 'EOG'  });
elseif EEG.nbchan == 66
    EEG = pop_select( EEG,'nochannel',{'M1' 'M2' 'EOG' 'VEOG' 'HEOG'}); 
elseif EEG.nbchan == 64
    EEG = pop_select( EEG,'nochannel',{'M1' 'M2' 'EOG'});
elseif EEG.nbchan == 63
    EEG = pop_select( EEG, 'nochannel', {'M1' 'M2'});
end
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 1 );
        
%% s2: remove bad channels, rereference, remove interpolated channels

load([curTask '_' subName '_badChan.mat']);
    
% remove bad channel
EEG = pop_select(EEG, 'nochannel', badChan);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',...
        [ subName '250HzLP30HP1Notch50BadChan']...
        ,'overwrite','on','gui','off');

% interpolate channels on dataset 3, based on dataset 1(dataset from step 1),
EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',[subName 'RMBadChanInter'],'overwrite','on','gui','off'); %Not saved

% Re-referece to common average
EEG = pop_reref( EEG, []);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',[subName 'RMBadChanInterReref'],'overwrite','on','gui','off'); %Not saved

% Remove bad channel again, then save set before next step
EEG = pop_select(EEG, 'nochannel', badChan);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname',[subName 'InterRerefRM'],'overwrite','on','gui','off'); %Not saved

%% s3: Epoch for ERP in ERPlab
binlistFolder = '/home/chendanni/Documents/Norms/analysis/Binlist'; % path to your binlist
binassname = '_binassigned.txt';

% Trigger/Binlist file name, if you have counterbalancing, set
% accordingly
iSubject = 1;
switch curTask
    case 'preLearningImplicit' 
        binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '_implicit.txt'];
        for i = 1:(size(EEG.urevent,2)-624-1)
            EEG.urevent(i).type = 'practice';
            EEG.event(i).codelabel = 'practice';
            EEG.event(i).type = 199;
        end
    case 'Learning'
        binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '.txt'];
    case 'postLearningMemory'
        binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '.txt'];
    case 'postLearningImplicit'
        binlistfile = ['SNBin' num2str(mod(subNameD(iSubject)-1200,14)) '_implicit.txt'];
        for i = 1:(size(EEG.urevent,2)-624-1)
            EEG.urevent(i).type = 'practice';
            EEG.event(i).codelabel = 'practice';
            EEG.event(i).type = 199;
        end
end

% Assign Binlister
binfile = strcat(binlistFolder,filesep,binlistfile);
outputeventfile = [outputSubjectFolder filesep subName '_' curTask binassname '_1Hz'];
inputeventfile = [outputSubjectFolder filesep subName '_' curTask '_eventlist.txt'];
EEG  = pop_binlister( EEG , 'BDF', binfile, 'ExportEL',...
    outputeventfile, 'ImportEL',...
    inputeventfile, 'IndexEL',  1, 'SendEL2', 'EEG&Text', 'Voutput', 'EEG'  );
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); 

% extract bin-based epochs
switch curTask
    case 'preLearningImplicit'
        EEG = pop_epochbin( EEG , [-200.0  1000.0],  'pre'); % set the epoch length and baseline correction window
    case 'Learning'
        EEG = pop_epochbin( EEG , [-300.0  3000.0],  'pre'); % set the epoch length and baseline correction window
    case 'postLearningMemory'
        EEG = pop_epochbin( EEG , [-200.0  1000.0],  'pre'); % set the epoch length and baseline correction window
    case 'postLearningImplicit'
        EEG = pop_epochbin( EEG , [-200.0  1000.0],  'pre'); % set the epoch length and baseline correction window
end
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'setname',[subName 'Epoch'],'overwrite','on','gui','off'); %Not saved

% reject bad epochs
load ([ subName '_' curTask '_rejectedEpoch.mat']);
EEG = pop_selectevent( EEG, 'omitevent',thisRejectedEpoch ,'deleteevents','off','deleteepochs','on','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'setname',[subName 'EpochRej'],'overwrite','on','gui','off'); %Not saved
pop_saveset( EEG, [ subName '_' curTask '_EpochRej1Hz.set'], pwd);


end
