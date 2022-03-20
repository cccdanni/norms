% ERSP analysis
% step 2: time-frequency analysis
% Project: Social Norms Learning
% Dependency: FieldTrip
% Author: Danni Chen 
% Update Date: 2022-03-01


clear all; clc; close all;
%% basic information setup

% Subject info
subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238]; % 1201 - 1204: HK Local; 1214, 1215 & 1238 - visually detected bad subjects
subNameD = setdiff(subNameD, badsubNameD);
task                = {'preLearningImplicit', 'postLearningImplicit'};

% folder info
workFolder          = '/home/chendanni/Documents/Norms/analysis/';
ERSPParentFolder    = fullfile (workFolder, 'EEGAnalysisResults', 'ERSP', 'Results');
ScriptsFolder       = fullfile(workFolder, 'MyScripts');
FieldTripFolder     = '/home/chendanni/MATLAB/toolbox/fieldtrip-20210330';
addpath ( genpath(ScriptsFolder) );
addpath ( genpath(FieldTripFolder) );

% parameter setup
bltw = [-0.200 0.000]; % in secs
toiw = [-0.200 1.000]; % in secs
bindescr = {'ingroup-higher', 'ingroup-lower', 'ingroup-consistent',...
        'outgroup-higher', 'outgroup-lower', 'outgroup-consistent'};

for iTask = 1:length(task)

    curTask = task{iTask};

    avg_in_high = {};
    avg_in_low  = {};
    avg_in_con  = {};
    avg_out_high = {};
    avg_out_low  = {};
    avg_out_con  = {};

    taskERSPFolder = fullfile (ERSPParentFolder, strcat(curTask,'Avg'));
    
    
    for iSubject = 1:length(subNameD)
        
        % load processed EEG data
        subName = num2str(subNameD(iSubject));
        subjectERSPfolder = fullfile (ERSPParentFolder, subName);
        cd (subjectERSPfolder);
        
        data = [];
        load (strcat(subName, '_', curTask, '_processed.mat'));
        
        
        for ibin = 1:length(bindescr)
            
            % define trial
            nTrials = find (data.trialinfo.bini ==  ibin)';
            
            % time-frequency analysis
            cfg = [];
            cfg.trials = nTrials;
            cfg.method = 'mtmconvol'; % implements multitaper time-frequency transformation based on multiplication in the frequency domain
            cfg.output = 'pow'; % return power-spectrm
            cfg.channel = 'all';
            cfg.foi = [1:30];
            cfg.taper = 'hanning';
            cfg.t_ftimwin = ones(length(cfg.foi),1).*0.2;
            % length of time window (in second), here I am using a fix time window (0.200 secs)
            % could also use a changing time window, i.e., nCycle./cfg.foi
            cfg.toi = -1:0.01:2; % time window slide from -1 to 2 secs in steps of 0.01 secs (10ms)
            freq = ft_freqanalysis (cfg, data);

            % baseline correction
            cfg.baseline = [-0.2 0];
            cfg.baselinetype = 'relative'; % data = data ./ meanVals;
            cfg.parameter = 'powspctrm';
            freq = ft_freqbaseline(cfg, freq);

            % save plotting
            cfg = [];
            cfg.baseline = 'no'; % already done it in previous code
            cfg.xlim = [-0.2 1];
            cfg.maskstyle = 'saturation';
            cfg.channel = 'Fz';
            cfg.layout = 'EEG1005.lay';
            cfg.showlabels = 'yes';
            cfg.colormap = '*RdBu';
            ft_singleplotTFR (cfg, freq);
            saveas(gcf,strcat(subName, '-', curTask,'-', bindescr{ibin}, '-Fz'),'png');
            close all;
            
            % save data
            curbin = bindescr{ibin};
            save(strcat(subName, '-', curTask,'-', bindescr{ibin},'.mat'), 'freq', 'curbin');

            % save in a big cell
            switch ibin
                case 1
                    avg_in_high{iSubject} = freq;
                case 2
                    avg_in_low{iSubject} = freq;
                case 3
                    avg_in_con{iSubject} = freq;
                case 4
                    avg_out_high{iSubject} = freq;
                case 5
                    avg_out_low{iSubject} = freq;
                case 6
                    avg_out_con{iSubject} = freq;
            end
            freq = [];

        end
        
        % time-frequency analysis, keep trials
        cfg = [];
        cfg.keeptrials = 'yes';
        cfg.method = 'mtmconvol'; % implements multitaper time-frequency transformation based on multiplication in the frequency domain
        cfg.output = 'pow'; % return power-spectrm
        cfg.channel = 'all';
        cfg.foi = [1:30];
        cfg.taper = 'hanning';
        cfg.t_ftimwin = ones(length(cfg.foi),1).*0.2;
        % length of time window (in second), here I am using a fix time window (0.200 secs)
        % could also use a changing time window, i.e., nCycle./cfg.foi
        cfg.toi = -1:0.01:2; % time window slide from -1 to 2 secs in steps of 0.05 secs (50ms)
        freq = ft_freqanalysis (cfg, data);
        
        % baseline correction
        cfg.baseline = [-0.2 0];
        cfg.baselinetype = 'relative'; % data = data ./ meanVals;
        cfg.parameter = 'powspctrm';
        freq = ft_freqbaseline(cfg, freq);

        % save data
        save(strcat(subName, '-', curTask,'-allfreq.mat'), 'freq');
        freq = [];

    end
    
    cd (taskERSPFolder);
    
    for ibin = 1:length(bindescr)
    
        switch ibin
            case 1
                avg_all = avg_in_high;
            case 2
                avg_all = avg_in_low;
            case 3
                avg_all = avg_in_con;
            case 4
                avg_all = avg_out_high;
            case 5
                avg_all = avg_out_low;
            case 6
                avg_all = avg_out_con;
        end
            
       
        % grand average
        cfg = [];
        grandavg = ft_freqgrandaverage(cfg, avg_in_high{:});
        curBin = bindescr{ibin};
        save (strcat( curTask,'-', curBin, '-grandavg.mat'), 'avg_all', 'grandavg', 'curBin', 'subNameD')
        
        % plot 
        cfg = [];
        cfg.baseline = 'no'; % already done it in previous code
        cfg.xlim = [-0.2 1];
        cfg.maskstyle = 'saturation';
        cfg.layout = 'EEG1005.lay';
        cfg.showlabels = 'yes';
        cfg.colormap = '*RdBu';
        figure;
        ft_multiplotTFR (cfg, grandavg);
        close all;

    end
      
    
end
