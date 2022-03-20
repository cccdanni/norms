%% ReadME
% Plot ERP Brainwaves
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: 8/30/2021

%% basic setting

% Subject info
clear all; clc;
subNameD = [1205:1242,1244:1253];
badsubNameD   = [1201:1204, 1214, 1215, 1238, 1210, 1211, 1223, 1239]; 
% 1201 - 1204: HK Local; 
% 1214, 1215 & 1238 - visually detected bad subjects
% 1210, 1211, 1223, 1239 - removed more than 100 trials
subNameD = subNameD(~ismember(subNameD, badsubNameD));

% folder info
workFolder = '/home/chendanni/Documents/Norms/analysis/';
task = {'Learning'};
ERPParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'ERP'); % set your output path
PreprocessParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing');
ScriptsFolder = '/home/chendanni/Documents/Norms/analysis/MyScripts/ERP';
FieldTripFolder = '/home/chendanni/MATLAB/toolbox/fieldtrip-20210330';
addpath ( genpath(ScriptsFolder) );
addpath ( genpath(FieldTripFolder) );
addpath ( genpath("/home/chendanni/Documents/Norms/analysis/MyScripts/additional_functions/boundedline-pkg-master") );


% ROI info
ROI = struct;
ROIs = {'FC','Ce','CP','LP','RP', 'LOT','ROT','OT'}; 
for i = 1:8 ROI(i).ROIName = ROIs{i}; end
ROI(1).ChanList = {'Fz', 'FCz', 'F1', 'F2', 'FC1', 'FC2'}; %FC
ROI(2).ChanList = {'Cz', 'C1', 'C2'}; %Ce
ROI(3).ChanList = {'CP1', 'CP2', 'Pz', 'P1', 'P2'}; %CP
ROI(4).ChanList = {'P3', 'P5', 'CP3', 'CP5'}; %LP
ROI(5).ChanList = {'P4', 'P6', 'CP4', 'CP6'}; %RP
ROI(6).ChanList = {'T7','TP7','P7','PO7'}; %LOT
ROI(7).ChanList = {'T8','TP8','P8','PO8'}; %ROT
ROI(8).ChanList = {'T7','TP7','P7','PO7','T8','TP8','P8','PO8'}; %OT


% other setting
fs = 250;
admtw = 0.05;

% type 
typelist = {'Mean', 'AdaptiveMean', 'Peak'};

cl=colormap(parula(50));
mycmap = [150 195 125; 243 210 102; 169 184 198]/255;
igcmap = [153 153 153; 40 120 181; 154 201 219]/255;
ogcmap = [153 153 153; 200 36 35; 248 172 140]/255;

bl = [-300 0];

%% Plot ERP
cnt = 1

for iTask = 1:length(task)
    
    curTask = task{iTask};
    
    % component of interest info
    switch curTask
        case 'Learning'
            bltw = [-0.30 0];
            tw   = [-300 2996];
            complist = {'n170','lpc','n400','frn','epn','lpp'};
            twlist = [0.12 0.22;...
              0.30 0.80;...
              0.30 0.50;...
              0.25 0.45;...
              0.23 0.28;...
              1.00 1.60];
            roilist = {{'LOT','ROT','OT'},...
               {'CP','FC','LP','RP'},...
               {'Ce','CP'},...
               {'FC'}...
               {'LOT','ROT','OT'}...
               {'LP', 'RP', 'CP', 'FC'}};
        case 'preLearningImplicit'
            bltw = [-0.20 0];
            tw   = [-200 996];
            complist = {'n170','lpc','n400','frn','epn'};
            twlist = [0.12 0.22;...
              0.30 0.80;...
              0.30 0.50;...
              0.25 0.45;...
              0.23 0.28];
            roilist = {{'LOT','ROT','OT'},...
               {'CP','FC','LP','RP'},...
               {'Ce','CP'},...
               {'FC'}...
               {'LOT','ROT','OT'}};
        case 'postLearningImplicit'
            bltw = [-0.20 0];
            tw   = [-200 996];
            complist = {'n170','lpc','n400','frn','epn'};
            twlist = [0.12 0.22;...
              0.30 0.80;...
              0.30 0.50;...
              0.25 0.45;...
              0.23 0.28];
            roilist = {{'LOT','ROT','OT'},...
               {'CP','FC','LP','RP'},...
               {'Ce','CP'},...
               {'FC'}...
               {'LOT','ROT','OT'}};
    end


    ERPName = strcat (curTask, '_GrandAvg_ManuRej', num2str(bl(1)), '.erp'); 
    ERPFolder = [ERPParentFolder '/' curTask 'ERPlabGrandAvg'];
    cd (ERPFolder);
    ERP = pop_loaderp('filename', ERPName, 'filepath', ERPFolder);
    allChan = char({ERP(:).chanlocs.labels});
    
    idx_start = nearest(ERP.times, tw(1));
    idx_end   = nearest(ERP.times, tw(2));
    
    
    %% ingroup: higher, lower v.s. consistent 
    in_higher_data = ERP.bindata(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Ingroup_Higher')));
    in_lower_data = ERP.bindata(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Ingroup_Lower ')));
    in_consistent_data = ERP.bindata(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Ingroup_Consistent')));
    in_higher_sem = ERP.binerror(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Ingroup_Higher')));
    in_lower_sem = ERP.binerror(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Ingroup_Lower ')));
    in_consistent_sem = ERP.binerror(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Ingroup_Consistent')));
    
    for iROI = 1:length(ROI)
        
        thisROIName = char(ROIs(iROI));
        thisChanName = char(ROI(iROI).ChanList);
        if length(thisROIName) < 3 
            thisROIName = pad(thisROIName, 3);
        end
        
        thisROWs         = find(ismember(allChan, thisROIName, 'rows'));
        tl_in_higher     = squeeze(in_higher_data(thisROWs, :));
        tl_in_lower      = squeeze(in_lower_data(thisROWs, :));
        tl_in_consistent = squeeze(in_consistent_data(thisROWs, :));
        sem_in_higher     = squeeze(in_higher_sem(thisROWs, :));
        sem_in_lower      = squeeze(in_lower_sem(thisROWs, :));
        sem_in_consistent = squeeze(in_consistent_sem(thisROWs, :));
        
        xtime = squeeze(ERP.times(idx_start:idx_end));
        figure();
        p = plot (xtime, tl_in_consistent,...
            xtime, tl_in_higher,...
            xtime, tl_in_lower,'--',...
            'LineWidth', 2)
        p(1).Color = igcmap(1,:)
        p(2).Color = igcmap(2,:)
        p(3).Color = igcmap(3,:)
        
        
        legend('consistent', 'higher','lower')

        % Add all our previous improvements:
        xlabel('Time (ms)');
        ylabel('Amplitude');
        title(strcat(curTask, '-', thisROIName, '-ingroup'));
        ax = gca();
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        ax.TickDir = 'out';
        box off;
        
        saveas(gcf,strcat('noSE_', curTask, thisROIName, 'IGERP', num2str(bl(1)), date()),'epsc');
        saveas(gcf,strcat('noSE_', curTask, thisROIName, 'IGERP', num2str(bl(1)), date()),'png');
        cnt = cnt + 1;
        close all;
        
    end
    
    %% outgroup: higher, lower v.s. consistent
    out_higher_data = ERP.bindata(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Outgroup_Higher ')));
    out_lower_data = ERP.bindata(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Outgroup_Lower')));
    out_consistent_data = ERP.bindata(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Outgroup_Consistent ')));
    out_higher_sem = ERP.binerror(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Outgroup_Higher ')));
    out_lower_sem = ERP.binerror(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Outgroup_Lower')));
    out_consistent_sem = ERP.binerror(:,idx_start:idx_end,find(strcmp(ERP.bindescr,'Outgroup_Consistent ')));
    
    for iROI = 1:length(ROI)
        
        thisROIName = char(ROIs(iROI));
        thisChanName = char(ROI(iROI).ChanList);
        if length(thisROIName) < 3 
            thisROIName = pad(thisROIName, 3);
        end
        
        thisROWs         = find(ismember(allChan, thisROIName, 'rows'));
        tl_out_higher     = squeeze(out_higher_data(thisROWs, :));
        tl_out_lower      = squeeze(out_lower_data(thisROWs, :));
        tl_out_consistent = squeeze(out_consistent_data(thisROWs, :));
        sem_out_higher     = squeeze(out_higher_sem(thisROWs, :));
        sem_out_lower      = squeeze(out_lower_sem(thisROWs, :));
        sem_out_consistent = squeeze(out_consistent_sem(thisROWs, :));
        
        xtime = squeeze(ERP.times(idx_start:idx_end));
        figure();
        p = plot (xtime, tl_out_consistent,...
            xtime, tl_out_higher,...
            xtime, tl_out_lower,'--',...
            'LineWidth', 2)
        p(1).Color = ogcmap(1,:)
        p(2).Color = ogcmap(2,:)
        p(3).Color = ogcmap(3,:)
        legend('consistent', 'higher','lower')

        % Add all our previous improvements:
        xlabel('Time (ms)');
        ylabel('Amplitude');
        title(strcat(curTask, '-', thisROIName, '-outgroup'));
        ax = gca();
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        ax.TickDir = 'out';
        box off;
        
        saveas(gcf,strcat('noSE_',curTask, thisROIName, 'OGERP', num2str(bl(1)), date()),'epsc');
        saveas(gcf,strcat('noSE_',curTask, thisROIName, 'OGERP', num2str(bl(1)), date()),'png');
        cnt = cnt + 1;
        close all;
    end
        
end