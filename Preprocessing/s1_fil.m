%% pre-processing EEG data
% s1: De-sampling, Filtering, Remove Channels
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and BrainVision .bv Data Loading Plugin
% Author: Danni Chen 
% Update Date: 2022-01-25
% Input: raw BrainVision data
% Output: *_250Hz05HP30LPNotch.set

function s1_fil (subNameD, task, sr, hp, lp, notch, linefreq, workFolder)

if nargin == 0
    
    % participant info set up
    subNameD = [1201:1242, 1244:1253];
    
    % directory set up
    workFolder = '/home/chendanni/Documents/Norms/analysis/';

    % task set up
    task = {'preLearningImplicit','postLearningImplicit','Learning'};
    
    % parameter set up
    sr = 250;
    lp = 30;
    hp = 0.05;
    notch = 1;
    linefreq = 50;
   
end

% folder set-up
rawFolder  = fullfile (workFolder, 'EEGRawData');
outputParentFolder = fullfile (workFolder, 'EEGAnalysisResults', 'Preprocessing');

% channel location file
channelLocationFile = '/home/chendanni/MATLAB/toolbox/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';% Channel location file(path according to your file location)
 

% load EEGlab
[ALLEEG EEG CURRENTSET] = eeglab;

    for iTask = 1:length(task)

        curTask = task{iTask};

        for iSubject = 1:length(subNameD)
        % parfor Execute for loop in parallel on workers in parallel pool
            
            ALLEEG = []; EEG = []; CURRENTSET = [];

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
            
            if subNameD(iSubject)==1239 & strcat (curTask, 'Learning')
                EEG = pop_loadset('filename','SN_Study1_Learning_sub1239.set','filepath',char(curRawFolder));
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            else
                curRawFile = curRawFile.name;
                EEG = pop_loadbv(curRawFolder, curRawFile); % pop_loadbv for .vhdr 
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname', setname,'gui','off');
            end
            
            setname = subName;
            cd(outputSubjectFolder);
            
            % down sample to SR: optional, default is 250Hz
            setname = strcat(setname, num2str(sr), 'Hz');
            EEG = pop_resample( EEG, sr); 
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',setname,'overwrite','on','gui','off'); %Not saved

            %low pass filter to lp for ERP: optional, default is 30
            setname = strcat(setname, num2str(lp), 'LP');
            EEG = pop_eegfiltnew(EEG, 'hicutoff',lp);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',setname,'overwrite','on','gui','off'); 

            %high pass filter to hp for ERP: optional, default is 0.05
            setname = strcat(setname, erase(num2str(hp),'0.'), 'HP');
            EEG = pop_eegfiltnew(EEG, 'locutoff',hp);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',setname,'overwrite','on','gui','off'); 

            %Notch filter 50Hz(for HK); filter function from ERPlab
            setname = strcat(setname, 'Notch');
            if notch 
                EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:63] ,'computepower',1,'linefreqs',linefreq,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',setname,'overwrite','on','gui','off'); 

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

            % Change saved file names if filter parameters are difference
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,...
                'setname',setname,...
                'savenew',strcat(subName,'_',curTask,'_', setname),'overwrite','on','gui','off'); 

        end
        
    end
    
end