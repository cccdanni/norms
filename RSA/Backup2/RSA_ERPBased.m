% This is a script for conducting RSA analysis for Norms Learning
% pipeline that utilizes a nested bin-epoched data structure.

% Author: Danni Chen, Ziqing Yao & Xiaoqing Hu (The University of Hong Kong)
% Update date: 11/15/2021

function RSA_ERPBased(subs, work_dir) 

%% Subject List: 
if nargin ==0
    
    subs = [1205:1242,1244:1253];
    badsubs   = [1201:1204, 1220, 1229, 1238, 1239, 1249]; % 1201 - 1204: Local, 1229, 1238, 1239,1249 (reject epochs > 100)
    subs = setdiff (subs, badsubs);
    
    work_dir = '/home/chendanni/Documents/Norms/analysis/';
    cd (work_dir);
    
    task = "postLearningImplicit";
    erp  = "full_range";
     
end

nSubs = length(subs); %number of subjects


%% Subject Filename(s)

% Load Subject Datafiles 
data_folder = fullfile (work_dir, 'EEGAnalysisResults', 'Preprocessing'); % set directory of data set (default: pwd)
toolbox_folder = '/home/chendanni/MATLAB/toolbox';
fName = strcat('/Decoding_BE_', task, '_'); % subject filename (subject # appended at end)

% Save RSA results
saveLocation = fullfile (pwd, 'EEGAnalysisResults','RSA');
if ~exist (saveLocation)
    mkdir (saveLocation);
end

% define the participants 
all_subs = {};
for i = 1:length(subs)
    all_subs (i) =  { num2str(subs(i)) } ;
end

%% define RSA parameters
% define sliding window
% win = 0.02; % in sec 
% slide = 0.01;
win = 1;
slide = 0;

% define new Sampling rate
sr=100; % in Hz
step=win/(1/sr)+1;

%define start and end of item window
t_start= 0;
t_end=0.996;
tois=t_start:1/sr:t_end;
t1=t_start:slide:(t_end-win);
t2=t1+win;
ind_t1=1:slide/(1/sr):((numel(tois)-win/(1/sr)));
ind_t2=ind_t1+win/(1/sr);
n_bins=numel(t1);

channels='all';

%% define path_in & path_out

path_in  = data_folder;
savName = strcat ('/', 'Results_ERPbasedRSA_', task ,'_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_'); 
path_out = fullfile(saveLocation, 'Results',strcat('ERPbasedRSA_', task, '_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide)));
mkdir(path_out)


%% add toolboxes
addpath (fullfile(work_dir,'MyScripts','additional_functions'));
addpath (fullfile(toolbox_folder,'fieldtrip-20210330'));

%% First step: correlate for each subject in all trials

cd (path_in);

for n = 1:numel(all_subs)
    
    sel_sub = all_subs{n};
    sel_sub_folder = fullfile( data_folder, sel_sub );
    
    
    %% EEG data loading & processing
    % load in the preprocessed data of each participant
    EEG = pop_loadset('filename',char(strcat( sel_sub, '_', task, '_EpochArtRej_FalseRespRej.set')),'filepath',fullfile(sel_sub_folder,'/'));
    data = eeglab2fieldtrip ( EEG, 'preprocessing', 'none' );
    events = EEG.event;
    
    % resample to sr Hz
    cfg=[];
    cfg.resamplefs=sr;
    cfg.detrend='no';
    data=ft_resampledata(cfg,data);
    
    % select channels if needed
    if ~strcmp(channels, 'all')
        cfg=[];
        cfg.channel=channels;
        data=ft_preprocessing(cfg,data);
    end

    % prepare task data (after stimuli onset)
    cfg=[];
    cfg.latency=[t_start t_end+1/sr];
    data=ft_selectdata( cfg, data );
    
    % remove unused trials if needed, here we removed object trials and
    % duplicated epochs    
    epochid = vertcat(events.epoch)';
    [uniqueid i j] = unique(epochid, 'first');
    indexdups = find( not (ismember(1:numel(epochid), i)));
    idex = 1:numel(epochid);
    idex = idex(~ismember (idex, indexdups) );
    events = events(idex); % clean up repeated epochs
    
    object_trials = find(arrayfun(@(events) ismember(9, events.bini), events));
    face_trials = find(arrayfun(@(events) ~ismember(9, events.bini), events));
    
    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=face_trials;
    data=ft_timelockanalysis(cfg, data);
    n_trials=numel(face_trials); % keep only face trials
    
    % z trans across trials
    mean_trials = repmat(nanmean(data.trial,1),n_trials,1,1);
    std_trials  = repmat(nanstd(data.trial,1),n_trials,1,1);
    data.trial  = (data.trial-mean_trials)./std_trials;
    clear n_trials mean_trials std_trials

    %% ERP data reorganizing
    % average trials by face_id
    codelabel = unique( [data.trialinfo.codelabel].', 'rows', 'stable'); 
    codeseq = [];
    
    for code_nr = 1: size(codelabel, 2)
        this_code = codelabel{code_nr};
        this_code = str2double ( this_code(2:end) );
        codeseq(code_nr,1) = this_code;
        codeseq(code_nr,2) = code_nr;
    end
    
    [~, idx] = sort(codeseq(:,1)); 
    codeseq  = codeseq (idx, :); % label sequence arranging 
    
    trialinfo = table2struct( data.trialinfo );

    erpdata = struct('time',{},'label', {}, 'elec', {}, 'avg', {}, 'var', {}, 'dof', {}, 'dimord', {}, 'cfg', {}, 'codelabel', {});
    for code_nr = 1:size(codelabel, 2)
        clear this_data this_code this_trials
        
        this_code = strcat('s', num2str(codeseq(code_nr, 1))); % this trial label
        this_trials = find(arrayfun(@(trialinfo) strcmp(this_code, trialinfo.codelabel), trialinfo));
        
        cfg=[];
        cfg.trials=this_trials;
        this_data = ft_timelockanalysis(cfg, data);
        this_data.codelabel = this_code;
       
        erpdata(code_nr)= this_data;
    end % re-arranging the erpdata into 80 rows, each row indicates one 
    
    erpdata_AVG = zeros(size(erpdata(1).label, 1),  size(erpdata(1).time, 2), size(erpdata,2));
    for code_nr = 1:size(erpdata, 2)
        this_erp = erpdata(code_nr).avg;
        erpdata_AVG(:, :, code_nr) = squeeze(this_erp);
    end  % average the ERP of each image
    
    erpdata_RESHAPE = reshape (erpdata_AVG, size(erpdata_AVG, 1) * size(erpdata_AVG, 2), size (erpdata_AVG, 3)); % reshape the ERP into (nChannels * nTp) * nBins
    
    
    %% Data correlation
    % correlation
    corr_res = corr(erpdata_RESHAPE, 'type', 'Spearman'); % now we would have one corrlation with nBins * nBins
   
    % fisher z transform correlations
    corr_res_ztrans=  0.5.*log(((ones(size(corr_res))+corr_res)./(ones(size(corr_res))-corr_res)));
    
    %% save data 
    save(fullfile(path_out, strcat(savName, all_subs{n},'.mat')),...
        'corr_res','corr_res_ztrans','erpdata','erpdata_AVG','erpdata_RESHAPE');
    
    clear corr_trials corr_cue_enc_trial trialinfo corr_cue_enc items corr_res cfg code_nr codelabel codeseq data EEG erpdata erpdata2 erpdata3 events face_trials object_trials idx i this_code this_data this erp this_trials
end % end of each subject 

%% ----------------------------------------------------------------------
%% Second step: Correlate with RDM

this_rsa_folder = path_out;
cd( this_rsa_folder );

allsubject_neuralrdm = zeros (80, 80, numel(subs));
allsubject_behardm = zeros (80, 80, numel(subs));
allsubject_coceprdm_conflict = zeros (80, 80, numel(subs));
allsubject_coceprdm_inconflict = zeros (80, 80, numel(subs));

allsubject = struct;

savName = strcat ('Results_ERPbasedRSA_', task ,'_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_'); 

if strcmp (task, 'preLearningImplicit')
    behadata = readmatrix (fullfile(work_dir, 'BehaData', 'PreLearningExplicitBehaData2021-02-23.csv'));
elseif strcmp (task, 'postLearningImplicit')
    behadata = readmatrix (fullfile(work_dir, 'BehaData', 'PostLearningExplicitBehaData2021-02-23.csv'));
end
behadata = behadata (:,[8,30,81]);

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

cnt = 1;

for sub_nr = 1:numel(subs)
    thissub = num2str(subs(sub_nr));
    
    cd(this_rsa_folder);
    load (strcat(savName, thissub,'.mat')); % corr_res_ztrans, erpdata, erpdata_RESHAPE
    clear corr_res; % clear raw correlation coefficient matrix in case of mistakes
    
    corr_res = 1 - corr_res_ztrans; % simialrity to RDM
    
    sel_sub_folder = fullfile( data_folder, thissub );
    EEG = pop_loadset('filename',char(strcat( thissub, '_', task, '_EpochArtRej_FalseRespRej.set')),'filepath',fullfile(sel_sub_folder,'/'));
   
    thisbini = [];
    thiscode = [];
    for i = 1: size(EEG.event, 2)
        thisbini(i) = EEG.event(i).bini;
        thislabel = EEG.event(i).codelabel;
        thiscode(i) = str2double(thislabel(2:end));
    end
    uniquecode = unique(thiscode, 'first');
    thisbini = thisbini (uniquecode);
    thiscode = thiscode (uniquecode);
    thisall  = [thisbini;thiscode]'; % Code and Corresponding Bin
    
    thisbinlist = zeros(8, 10);
    for i = 1:8 % 8 bins
        x = ceil((min(thiscode(find(thisbini == i))))/10);
        x1 = ((x-1)*10+1):x*10;
        thisbinlist(i,:) = x1;
    end % check each bin includes which figures
    
    newseq = reshape(thisbinlist',1,80);
    neural_rdm = zeros (size(corr_res,1), size(corr_res,2));
    for i = 1: size(corr_res,1)
        for j = 1:size(corr_res,2)
            neural_rdm(i,j) = corr_res (newseq(i),newseq(j));
        end
    end % neural_rdm, rearrange corr_res by bin
    allsubject_neuralrdm (:,:, sub_nr) = neural_rdm;
    
    conceptual_rdm_conflict = ones (size(corr_res,1), size(corr_res,2));
    for i = 1: size(corr_res,1)
        for j = 1:size(corr_res,2)
            if j>=i 
                conceptual_rdm_conflict (i, j ) = nan;
            elseif (ismember(j, [1:10, 31:40])) & (ismember(i, [71:80]))
                    conceptual_rdm_conflict (i, j ) = 0;
            end
        end
    end % conceptual rdm 1: both ingroup / outgroup higher is similar to morphed faces
    allsubject_coceprdm_conflict (:,:, sub_nr) = conceptual_rdm_conflict;
    
    conceptual_rdm_inconflict = ones (size(corr_res,1), size(corr_res,2));
    for i = 1: size(corr_res,1)
        for j = 1:size(corr_res,2)
            if j>=i 
                conceptual_rdm_inconflict (i, j ) = nan;
            elseif (ismember(j, [1:10])) & (ismember(i, [71:80]))
                    conceptual_rdm_inconflict (i, j ) = 0;
            end
        end
    end % conceptual rdm 2: only ingroup higher is similar to morphed faces
    allsubject_coceprdm_inconflict (:,:, sub_nr) = conceptual_rdm_inconflict;
    
    behavior_rdm = zeros (size(corr_res,1), size(corr_res,2));
    thisbeha = behadata(find(behadata(:,3) == subs(sub_nr)),:);
    for i = 1:size(corr_res,1)
        for j = 1:size(corr_res,2)
            rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
            rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
            behavior_rdm(i,j) = abs(rating_i - rating_j);
        end
    end % behavioral rdm
    allsubject_behardm (:,:,sub_nr) = behavior_rdm;
    
    
    for i = 1:6
        allsubject(cnt).SubjectID = subs(sub_nr);
        allsubject(cnt).Session   = i; 
        
        switch i 
            case 1 
                allsubject(cnt).Group = 1; allsubject(cnt).Compare = 1; 
            case 2 
                allsubject(cnt).Group = 1; allsubject(cnt).Compare = 2; 
            case 3 
                allsubject(cnt).Group = 1; allsubject(cnt).Compare = 3;
            case 4 
                allsubject(cnt).Group = 2; allsubject(cnt).Compare = 1; 
            case 5 
                allsubject(cnt).Group = 2; allsubject(cnt).Compare = 2; 
            case 6 
                allsubject(cnt).Group = 2; allsubject(cnt).Compare = 3; 
        end
        
        thiscondNeural = neural_rdm (71:80, ((i-1)*10+1):i*10);
        thiscondBeha   = behavior_rdm (71:80, ((i-1)*10+1):i*10);
        allsubject(cnt).meanNeuralRDM2Morph = mean(thiscondNeural, 'all');
        allsubject(cnt).meanBehaRDM2Morph   = mean(thiscondBeha, 'all');
        
        thiscondNeural = reshape(thiscondNeural, size(thiscondNeural,1)*size(thiscondNeural,2), 1 );
        thiscondBeha = reshape(thiscondBeha, size(thiscondBeha,1)*size(thiscondBeha,2), 1);
        allsubject(cnt).NeuralBehaCorr      = corr(thiscondNeural, thiscondBeha, 'type','Spearman');
        
        cnt = cnt +1;
    end

    cd(this_rsa_folder);
    save(strcat(savName,thissub, '_trans.mat'),...
        'neural_rdm', 'thisbinlist',...
        'conceptual_rdm_conflict', 'conceptual_rdm_inconflict', 'behavior_rdm');
    clear thissub corr_res corr_res_update thisbinlist newseq;
end

save(strcat(savName, 'all_trans.mat'), 'allsubject_behardm', 'allsubject_neuralrdm',....
    'allsubject_coceprdm_inconflict', 'allsubject_coceprdm_conflict', 'allsubject');

%% ----------------------------------------------------------------------
%% Third step: Visualization Neural RDM
allsubject_neuralrdm_mean = mean(allsubject_neuralrdm, 3);

CustomXLabels = string ([1:8]);
for i = 1:8 CustomXLabels(i) = ""; end
CustomXLabels(1) = 'In-Higher'
CustomXLabels(2) = 'In-Lower'
CustomXLabels(3) = 'In-Consistent'
CustomXLabels(4) = 'Out-Higher'
CustomXLabels(5) = 'Out-Lower'
CustomXLabels(6) = 'Out-Consistent'
CustomXLabels(7) = 'Control'
CustomXLabels(8) = 'Morph'

figure;
set (gcf, 'Position', [300 300 800 700]);
imagesc (allsubject_neuralrdm_mean)
colorbar
title (strcat(savName, ' Neural Dissimilarity'))
set(gca,'xtick',[5:10:80],'xticklabel',CustomXLabels,...
    'ytick',[5:10:80],'yticklabel',CustomXLabels)
xtickangle(45)
hold on 
for i = [0:10:80]
    for j = [0:10:80]
        rectangle ('Position', [i+.5 j+.5  10 10], 'LineWidth', 2);
    end
end

saveas(gcf,strcat(savName, '_all_trans.jpg'), 'jpg')
close all;


end