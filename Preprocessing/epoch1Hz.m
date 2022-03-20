function EEG = epoch1Hz(subjectID, curTask, rawFolder, workFolder, sr, lp, hp, notch, linefreq, channelLocationFile)

% Function: apply same parameters except high-passed 1 Hz
% Arg:
%   Input: subjectID, curTask
%   Output: EEG

if nargin == 0

    % participant info set up
    subjectID = 1201;
    
    % directory set up
    workFolder = '/Users/danni/Desktop/Norms/'; % in MAC
    rawFolder  = fullfile (workFolder, 'EEGRawData');

    % task set up
    curTask = 'preLearningImplicit';
    
    % parameter set up
    sr = 250;
    lp = 30;
    hp = 1;
    notch = 1;
    linefreq = 50;
    
    % channel location file
    channelLocationFile = '/home/chendanni/MATLAB/toolbox/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';% Channel location file(path according to your file location)
    
end

outputParentFolder = fullfile(workFolder, 'EEGAnalysisResults', 'Preprocessing');

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

ALLEEG = []; EEG = []; CURRENTSET = [];
subName = num2str(subjectID);
outputSubjectFolder = fullfile(outputParentFolder,subName);
 
%% s1: De-sampling, Filtering, Remove Channels

% load data from the rawFolder: for ANT.cnt: loadeep_v4; for other file
% type, manually load data and check the function by 'eegh'
curRawFolder = fullfile(rawFolder, ['sub', subName]);
cd (curRawFolder);
curRawFile = dir(['SN_Study1_' curTask '_sub', subName,'*.vhdr']);
            
if subNameD(iSubject)==1239 & strcat (curTask, 'Learning')
    EEG = pop_loadset('filename','SN_Study1_Learning_sub1239.set','filepath',char(curRawFolder));
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
else
    curRawFile = curRawFile.name;
    EEG = pop_loadbv(curRawFolder, curRawFile); % pop_loadbv for .vhdr 
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname', setname,'gui','off');
end
            
cd (outputSubjectFolder);

setname = subName;
            
% down sample to SR: optional, default is 250Hz
setname = strcat(setname, num2str(sr), 'Hz');
EEG = pop_resample( EEG, sr); 

% low pass filter to lp for ERP: optional, default is 30
setname = strcat(setname, num2str(lp), 'LP');
EEG = pop_eegfiltnew(EEG, 'hicutoff',lp);

% high pass filter to hp for ERP: optional, default is 1 Hz (for ICA)
setname = strcat(setname, '1HP');
EEG = pop_eegfiltnew(EEG, 'locutoff',hp);

% Notch filter 50Hz(for HK); filter function from ERPlab
setname = strcat(setname, 'Notch');
if notch 
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:63] ,'computepower',1,'linefreqs',linefreq,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
end

% load channel location file
EEG = pop_chanedit(EEG, 'lookup',channelLocationFile);
EEG = eeg_checkset( EEG );

% exclude channels not used
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

%% s3: remove big artifact and rereference
% remove bad section
badSecFile = strcat(subName, '_', task, '_BadSec.mat');
if exist (badSecFile)
    load(badSecFile);
    badSection = TMPREJ(:,1:2);        
    EEG = eeg_eegrej( EEG, badSection);
end
clear TMPREJ;

% remove bad channels & intepolate 
EEGrej = pop_loadset('filename',[ subName '_' curTask '_EpochRej.set'],'filepath',pwd);
rejected = EEGrej.chanlocs;
rejected = [rejected(:).urchan];
original = TMP.chanlocs;
original = [original(:).urchan];
badChan = setdiff(rejected, original);

if ~isnan (badChan)

    % remove bad channel
    EEG = pop_select(EEG, 'nochannel', badChan);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',...
        strcat (setName, 'BadChan'),'overwrite','on','gui','off');
    cnt_rm = cnt;
    cnt = cnt + 1;

    % interpolate channels on dataset removed bad channel, based on dataset 1 (dataset from step 1),
    EEG = pop_interp(EEG, ALLEEG(1).chanlocs, 'spherical');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',...
        strcat(subName, '_Inter.set'),'overwrite','on','gui','off');

    % re-referece to common average
    EEG = pop_reref( EEG, []);

    % remove bad channel again, then save set before next step
    EEG = pop_select(EEG, 'nochannel', badChan);

else

    % re-referece to common average
    EEG = pop_reref( EEG, []);

end

%% s4: epoch (longer than needed)
% bin name
binfile = fullfile (binlistFolder,binlistfile);
outputeventfile = fullfile (outputSubjectFolder, strcat(subName, '_', curTask, binassname));
inputeventfile  = fullfile (outputSubjectFolder, strcat(subName, '_', curTask, '_eventlist2.txt'));

% Assign Binlister    
EEG  = pop_binlister( EEG , 'BDF', binfile,... 
    'ExportEL', outputeventfile,... 
    'ImportEL', inputeventfile,... 
    'IndexEL',  1, 'SendEL2', 'All', 'Voutput', 'EEG' ); 

% extract bin-based epochs
% set the epoch length and baseline correction window
switch curTask
    case 'preLearningImplicit'
        EEG = pop_epochbin( EEG , [-1000.0  2000.0],  'none'); 
    case 'Learning'
        EEG = pop_epochbin( EEG , [-1000.0  4000.0],  'none'); 
    case 'postLearningImplicit'
        EEG = pop_epochbin( EEG , [-1000.0  2000.0],  'none'); 
end

%% s5: manually detect and remove bad epochs
badEpochFile = strcat(subName, '_', task, '_BadEpoch.mat');
if exist (badSecFile)
    load(badEpochFile);
    EEG = pop_rejepoch( EEG, badEpoch, 0);
end

end
