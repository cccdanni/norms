%% ReadME
% Evaluate ERP Data Quality & Extract Data 
% Project: Social Norms Learning
% Dependency: EEGlab (2019 version), ERPlab, and ANT .cnt Data Loading Plugin
% Author: Danni Chen 
% Update Date: Jan-14-2021


clear all; clc;

subNameD = [1205:1220]; 
badsubNameD   = [1201:1204, 1228, 1237, 1239, 1229]; % 1201 - 1204: Local, 1228, 129, 1237, 1239 (remained less than 50%)
subNameD = subNameD(~ismember(subNameD, badsubNameD));

workFolder = '/home/chendanni/Documents/Norms/analysis/';
cd(workFolder)
rawFolder = fullfile(workFolder, 'EEGRawData'); 
%task = {'preLearningImplicit', 'Learning', 'postLearningImplicit'};
task = {'preLearningImplicit', 'postLearningImplicit'}

% change outputParentFolder before ERP OR ERSP
outputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'ERP'); % set your output path
inputParentFolder = fullfile(workFolder,'EEGAnalysisResults', 'Preprocessing')
% change accordingly to your eeglab version and location

ROIList = struct;
ROIs = {'FC','Ce','CP','LP','RP','BiP','LOT','ROT','OT'}; 
for i = 1:8 ROIList(i).ROIName = ROIs{i}; end
ROIList(1).ChanList = {'Fz', 'FCz', 'F1', 'F2', 'FC1', 'FC2'}; %FC
ROIList(2).ChanList = {'Cz', 'C1', 'C2'}; %Ce
ROIList(3).ChanList = {'CP1', 'CP2', 'Pz', 'P1', 'P2'}; %CP
ROIList(4).ChanList = {'P3', 'P5', 'CP3', 'CP5'}; %LP
ROIList(5).ChanList = {'P4', 'P6', 'CP4', 'CP6'}; %RP
ROIList(6).ChanList = {'T7','TP7','P7','PO7'}; %LOT
ROIList(7).ChanList = {'T8','TP8','P8','PO8'}; %ROT
ROIList(8).ChanList = {'T7','TP7','P7','PO7','T8','TP8','P8','PO8'}; %OT

fs = 250;
bltw = 0.2;
tktw = 1;

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iROI = 1:length(ROIList)
    
    thisROI = ROIList(iROI).ROIName;
    thisROIChans = ROIList(iROI).ChanList;
    
    for iTask = 1:length(task)
        
        ERPAll = struct; 
        
        curTask = task{iTask};
        
        ERPsetList = [curTask thisROI 'GrandAvgAll_epoch.txt']; % create a txt file to save the subject level ERP sets to be involved in grand averge
        
        cnt = 1;
        
        for iSubject = 1:length(subNameD)
            
            %number to string
            subName = num2str(subNameD(iSubject));
            
            %locate sub folder
            outputSubjectFolder = [outputParentFolder subName];
            
            % load averaged ERP
            ERP = pop_loaderp( 'filename', [subName '_' curTask thisROI '_ERPAvg.erp'],  'filepath', outputSubjectFolder );
            
            % average ERP
            thisERP = mean(ERP.bindata,1);
            
            % baseline correction 
            blmean = mean(thisERP(1:fs*bltw))
            thisERP = thisERP-blmean;
            indec = 1:length(thisERP);
            
            % N170 (LOT, ROT, OT; 120-220ms)
            ERPN170 = thisERP(:, ((bltw + 0.12) * fs) : ((bltw + 0.22) *fs), :);
            indecN170 = indec(1, ((bltw + 0.12) * fs) : ((bltw + 0.22) *fs));
            admtwN170 = (round(0.05 * fs)-1) /2;
            
            % LPC (CP, FC; 300 - 800ms)
            ERPLPC = thisERP(:, ((bltw + 0.3) * fs) : ((bltw + 0.8) *fs), :);
            indecLPC = indec(1, ((bltw + 0.3) * fs) : ((bltw + 0.8) *fs));
            admtwLPC = (round(0.1 * fs)-1) /2;
            
            % N400 (N400, Ce, CP, 300 - 500ms)
            ERPN400 = thisERP(:, ((bltw + 0.3) * fs) : ((bltw + 0.5) *fs), :);
            indecN400 = indec(1, ((bltw + 0.3) * fs) : ((bltw + 0.5) *fs));
            admtwN400 = (round(0.1 * fs)-1) /2;
            
            for cond = 1:9
               
                ERPAll(cnt).SubjectID = subName;
                ERPAll(cnt).Cond         = ERP.bindescr(cond);
                ERPAll(cnt).ROI            = thisROI;
                ERPAll(cnt).ROIChan   = thisROIChans; 
                ERPAll(cnt).task   = curTask; 
                
                ERPAll(cnt).N170_Peak = max(ERPN170(:,:,cond));
                thisPeak = find(ERPN170(:,:,cond) == max(ERPN170(:,:,cond)));
                ERPAll(cnt).N170_AdaptiveMean = mean(thisERP(:, (indecN170(thisPeak)-admtwN170): (indecN170(thisPeak)+admtwN170),cond));
                ERPAll(cnt).N170_Mean = mean(ERPN170(:,:,cond));
                
                ERPAll(cnt).LPC_Peak = max(ERPLPC(:,:,cond));
                thisPeak = find(ERPLPC(:,:,cond) == max(ERPLPC(:,:,cond)));
                ERPAll(cnt).LPC_AdaptiveMean = mean(thisERP(:, (indecLPC(thisPeak)-admtwLPC): (indecLPC(thisPeak)+admtwLPC),cond));
                ERPAll(cnt).LPC_Mean = mean(ERPLPC(:,:,cond));
                
                ERPAll(cnt).N400_Peak = max(ERPN400(:,:,cond));
                thisPeak = find(ERPN400(:,:,cond) == max(ERPN400(:,:,cond)));
                ERPAll(cnt).N400_AdaptiveMean = mean(thisERP(:, (indecN400(thisPeak)-admtwN400): (indecN400(thisPeak)+admtwN400),cond));
                ERPAll(cnt).N400_Mean = mean(ERPN400(:,:,cond));
                
                cnt = cnt + 1
            
            end
        
        end
        
        outputERPFolder = [outputParentFolder '/' curTask 'ERPlabGrandAvg'];
        cd(outputERPFolder);
        
        ERP = pop_loaderp( 'filename', [curTask thisROI '_GrandAvg_ManuRej.erp'],  'filepath', outputERPFolder );
        
        % average ERP
        thisERP = squeeze(mean(ERP.bindata,1));
        thisError = squeeze(mean(ERP.binerror,1));
        
        writetable(struct2table(ERPAll), ['ERPAll' curTask thisROI date '.xlsx']);
        save(['ERPAll' curTask thisROI date '.mat'], 'thisERP', 'thisError');
        
    end

end

