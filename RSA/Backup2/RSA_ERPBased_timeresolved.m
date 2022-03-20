% This is a script for conducting RSA analysis for Norms Learning
% pipeline that utilizes a nested bin-epoched data structure.

% Author: Danni Chen, Ziqing Yao & Xiaoqing Hu (The University of Hong Kong)
% Update date: 11/15/2021

function RSA_ERPBased_timeresolved(subs, work_dir) 

%% Subject List: 
if nargin ==0
    
    subs = [1205:1242,1244:1253];
    badsubs   = [1201:1204, 1220, 1229, 1238, 1239, 1249]; % 1201 - 1204: Local, 1229, 1238, 1239,1249 (reject epochs > 100)
    subs = setdiff (subs, badsubs);
   
    work_dir = '/home/chendanni/Documents/Norms/analysis/';
    cd (work_dir);
    
    task = "postLearningImplicit";
    erp  = "full_range";
    
    stat_test = "yes"; 
    perm_test = "no";
    ztrans = "no"; 
    
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
win = 0.1;
slide = 0.01;

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
savName = strcat ('/', 'Results_ERPbasedRSA_', task ,'_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 
path_out = fullfile(saveLocation, 'Results',strcat('ERPbasedRSA_', task, '_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
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
    if (strcmp (ztrans, 'yes'))
        mean_trials = repmat(nanmean(data.trial,1),n_trials,1,1);
        std_trials  = repmat(nanstd(data.trial,1),n_trials,1,1);
        data.trial  = (data.trial-mean_trials)./std_trials;
        clear n_trials mean_trials std_trials 
    end
    % As we applied a z-transformation to the EEG activity at every time point 
    % across all trials of all conditions in order to remove ERP componentsprior 
    % to calculating correlations (Fellner et al., 2020, Current Biology)

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
    end % re-arranging the erpdata into 80 rows, each row indicates one erp (has been averaged by fieldtrip) 
    
    erpdata_AVG = zeros(size(erpdata(1).label, 1),  size(erpdata(1).time, 2), size(erpdata,2));
    for code_nr = 1:size(erpdata, 2)
        this_erp = erpdata(code_nr).avg;
        erpdata_AVG(:, :, code_nr) = squeeze(this_erp);
    end  % extract erp data from erpdata and make a (channel*timepoints*stimuli) matrix
    
    
    %% time-resolved ERP and RSA
    erpdata_AVG_bins_all = nan (size (erpdata_AVG,1), step, size (erpdata_AVG,3), n_bins);
    erpdata_RESHAPE_bins_all = nan (size(erpdata_AVG, 1) * step, size (erpdata_AVG, 3), n_bins);
    corr_res_all = nan ( size (erpdata_AVG, 3), size (erpdata_AVG, 3), n_bins);
    corr_res_ztrans_all = nan ( size (erpdata_AVG, 3), size (erpdata_AVG, 3), n_bins);
    
    for nbin = 1:n_bins
        
        erpdata_AVG_thisbin = erpdata_AVG(:, ind_t1(nbin):ind_t2(nbin), :);
        erpdata_RESHAPE_thisbin = reshape (erpdata_AVG_thisbin, size(erpdata_AVG_thisbin, 1) * size(erpdata_AVG_thisbin, 2), size (erpdata_AVG_thisbin, 3)); % reshape the ERP into (nChannels * nTp) * nBins
        
        %% Data correlation
        % correlation
        corr_res_thisbin = corr(erpdata_RESHAPE_thisbin, 'type', 'Spearman'); % now we would have one corrlation with nBins * nBins
        % fisher z transform correlations
        corr_res_ztrans_thisbin =  0.5.*log(((ones(size(corr_res_thisbin))+corr_res_thisbin)./(ones(size(corr_res_thisbin))-corr_res_thisbin)));
        
        %% save data 
        erpdata_AVG_bins_all(:,:,:,nbin) = erpdata_AVG_thisbin;
        erpdata_RESHAPE_bins_all(:,:,nbin) = erpdata_RESHAPE_thisbin;
        corr_res_all(:,:,nbin) = corr_res_thisbin;
        corr_res_ztrans_all(:,:,nbin) = corr_res_ztrans_thisbin;

        clear erpdata_AVG_thisbin erpdata_RESHAPE_thisbin corr_res_thisbin corr_res_ztrans_thisbin;
    end
    
    

    %% save data 
    save(fullfile(path_out, strcat(savName, all_subs{n},'.mat')),...
        'corr_res_all','corr_res_ztrans_all','erpdata_AVG_bins_all','erpdata_RESHAPE_bins_all');
    
    clear corr_trials corr_cue_enc_trial trialinfo corr_cue_enc items corr_res cfg code_nr codelabel codeseq data EEG erpdata erpdata2 erpdata3 events face_trials object_trials idx i this_code this_data this erp this_trials
end % end of each subject 

%% ----------------------------------------------------------------------
%% Second step: Correlate with Conceptual and Behavioral RDM

this_rsa_folder = path_out;
cd( this_rsa_folder );

allsubject_neuralrdm_allbins = zeros (80, 80, n_bins, numel(subs));
allsubject_behardm_allbins = zeros (80, 80, n_bins, numel(subs));
allsubject_coceprdm_conflict_allbins = zeros (80, 80, n_bins, numel(subs));
allsubject_coceprdm_inconflict_allbins = zeros (80, 80, n_bins, numel(subs));

allsubject = struct;

savName = strcat ('Results_ERPbasedRSA_', task ,'_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 

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
    clear corr_res_all; % clear raw correlation coefficient matrix in case of mistakes
    
    corr_res_all = 1 - corr_res_ztrans_all; % simialrity to RDM
    
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
    
    conceptual_rdm_conflict = ones (size(corr_res_all,1), size(corr_res_all,2));
    for i = 1: size(corr_res_all,1)
        for j = 1:size(corr_res_all,2)
            if j>=i 
                conceptual_rdm_conflict (i, j ) = nan;
            elseif (ismember(j, [1:10, 31:40])) & (ismember(i, [71:80]))
                    conceptual_rdm_conflict (i, j ) = 0;
            end
        end
    end % conceptual rdm 1: both ingroup / outgroup higher is similar to morphed faces
    
    conceptual_rdm_inconflict = ones (size(corr_res_all,1), size(corr_res_all,2));
    for i = 1: size(corr_res_all,1)
        for j = 1:size(corr_res_all,2)
            if j>=i 
                conceptual_rdm_inconflict (i, j ) = nan;
            elseif (ismember(j, [1:10])) & (ismember(i, [71:80]))
                    conceptual_rdm_inconflict (i, j ) = 0;
            end
        end
    end % conceptual rdm 2: only ingroup higher is similar to morphed faces
    
    behavior_rdm = zeros (size(corr_res_all,1), size(corr_res_all,2));
    thisbeha = behadata(find(behadata(:,3) == subs(sub_nr)),:);
    for i = 1:size(corr_res_all,1)
        for j = 1:size(corr_res_all,2)
            rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
            rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
            behavior_rdm(i,j) = abs(rating_i - rating_j);
        end
    end % behavioral rdm

    for nbin = 1:n_bins
        corr_res = corr_res_all(:,:,nbin)
        
        neural_rdm = zeros (size(corr_res,1), size(corr_res,2));
        for i = 1: size(corr_res,1)
            for j = 1:size(corr_res,2)
                neural_rdm(i,j) = corr_res (newseq(i),newseq(j));
            end
        end % neural_rdm, rearrange corr_res by bin
        allsubject_neuralrdm_allbins (:, :, nbin, sub_nr) = neural_rdm;
        allsubject_behardm_allbins (:, :, nbin, sub_nr) = behavior_rdm;
        allsubject_coceprdm_conflict_allbins (:, :, nbin, sub_nr) = conceptual_rdm_conflict;
        allsubject_coceprdm_inconflict_allbins (:, :, nbin, sub_nr) = conceptual_rdm_inconflict;
        
    end % end of bins
        
 
end % end of subjects

save(strcat(savName, 'trans.mat'), 'allsubject_neuralrdm_allbins', 'allsubject_behardm_allbins',....
    'allsubject_coceprdm_conflict_allbins', 'allsubject_coceprdm_inconflict_allbins');

%% ----------------------------------------------------------------------

%% Conduct statical analysis

if strcmp (stat_test, 'yes')
    
    %% Data Prep. 
    
    load (strcat(savName, 'trans.mat'));
    % allsubject_behardm_allbins: nItems x nItems x nBins x nSubject
    % allsubject_conceprdm_conflict_allbins: nItems x nItems x nBins x nSubject
    % allsubject_conceprdm_inconflict_allbins: nItems x nItems x nBins x nSubject
    % allsubject_neuralrdm_allbin: nItems x nItems x nBins x nSubject
    
    n_bins = size (allsubject_neuralrdm_allbins, 3);
    
    
    allsubject_neuralrdm_allbin_AVG = nan (size (allsubject_neuralrdm_allbins, 4), 8, n_bins);
    % nSubject x 7 (1 + 6 nConds) x nBin
    allsubject_neuralbeha_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralbeha_P = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralconp_all_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralconp_all_P = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralconp_in_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralconp_in_P = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    % nSubject x Bin
    
    
    for nSub = 1: size (allsubject_neuralrdm_allbins, 4)
        
        for nBin = 1: size (allsubject_neuralrdm_allbins, 3)
            
            behardm     = allsubject_behardm_allbins (:, :, nBin, nSub);
            conprdm_all = allsubject_coceprdm_conflict_allbins (:, :, nBin, nSub);
            conprdm_in  = allsubject_coceprdm_inconflict_allbins (:, :, nBin, nSub);
            neuralrdm   = allsubject_neuralrdm_allbins (:, :, nBin, nSub);
            
            %% Extracting the neural dissimilarity only
            for nCond = 1: 6
                allsubject_neuralrdm_allbin_AVG (nSub, nCond, nBin) = mean (neuralrdm(((nCond-1)*10 +1):nCond*10, 71:80) , 'all');
            end
            allsubject_neuralrdm_allbin_AVG (nSub, 7, nBin) = nSub;
            allsubject_neuralrdm_allbin_AVG (nSub, 8, nBin) = nBin;
            
            %% Correlation between neural rdm and conceptual rdm (in and out conflict)
            conprdm_VECTOR = matrix2vector (conprdm_all, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            [neuralconp_all_corr, neuralconp_all_corr_p] = corr(neuralrdm_VECTOR', conprdm_VECTOR', 'type','Spearman');
            allsubject_neuralconp_all_R (nSub, nBin) = neuralconp_all_corr;
            allsubject_neuralconp_all_P (nSub, nBin) = neuralconp_all_corr_p;
            
            %% Correlation between neural rdm and conceptual rdm (in conflict only)
            conprdm_VECTOR = matrix2vector (conprdm_in, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            [neuralconp_in_corr, neuralconp_in_corr_p]  = corr(neuralrdm_VECTOR', conprdm_VECTOR', 'type','Spearman');
            allsubject_neuralconp_in_R (nSub, nBin) = neuralconp_in_corr;
            allsubject_neuralconp_in_P (nSub, nBin) = neuralconp_in_corr_p;
            
            %% Correlation between neural rdm and behavioral rdm 
            behardm = matrix2vector (behardm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            [neuralbeha_corr, neuralbeha_corr_p]  = corr(neuralrdm_VECTOR', behardm', 'type','Spearman');
            allsubject_neuralbeha_R (nSub, nBin) = neuralbeha_corr;
            allsubject_neuralbeha_P(nSub, nBin) = neuralbeha_corr_p;
            
        end
        
    end
    
    Results_NeuralConpAll = struct;
    Results_NeuralConpIn  = struct;
    Results_NeuralBeha    = struct; 
    Results_Diff          = struct;
    
    %% Significant Analysis 
    for nBin = 1:n_bins
        
        %% Correlation between neural rdm and conceptual rdm (in and out conflict)
        [Results_NeuralConpAll(nBin).H, Results_NeuralConpAll(nBin).P] = ttest (allsubject_neuralconp_all_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralConpAll(nBin).R = mean (allsubject_neuralconp_all_R(:,nBin));
        Results_NeuralConpAll(nBin).SE = std (allsubject_neuralconp_all_R(:,nBin)) / length (allsubject_neuralconp_all_R(:,nBin));
        
        %% Correlation between neural rdm and conceptual rdm (in conflict only)
        [Results_NeuralConpIn(nBin).H, Results_NeuralConpIn(nBin).P] = ttest (allsubject_neuralconp_in_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralConpIn(nBin).R = mean (allsubject_neuralconp_in_R(:,nBin));
        Results_NeuralConpIn(nBin).SE = std (allsubject_neuralconp_in_R(:,nBin)) / length (allsubject_neuralconp_in_R(:,nBin));
   
        %% Correlation between neural rdm and behavioral rdm 
        [Results_NeuralBeha(nBin).H, Results_NeuralBeha(nBin).P] = ttest (allsubject_neuralbeha_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralBeha(nBin).R = mean (allsubject_neuralbeha_R(:,nBin));
        Results_NeuralBeha(nBin).SE = std (allsubject_neuralbeha_R(:,nBin)) / length (allsubject_neuralbeha_R(:,nBin));
        
        %% Differences among different conditions
        thisBinNeuralRDM = allsubject_neuralrdm_allbin_AVG(:,:,nBin);
        
        [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,2));
        Results_Diff(nBin).inHL_H = H;
        Results_Diff(nBin).inHL_P = P;
        Results_Diff(nBin).inHL_T = STATS.tstat; % ingroup higher vs. lower

        [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,3));
        Results_Diff(nBin).inHC_H = H;
        Results_Diff(nBin).inHC_P = P;
        Results_Diff(nBin).inHC_T = STATS.tstat; % ingroup higher vs. consistent 
        
        [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,2), thisBinNeuralRDM(:,3));
        Results_Diff(nBin).inLC_H = H;
        Results_Diff(nBin).inLC_P = P;
        Results_Diff(nBin).inLC_T = STATS.tstat; % ingroup lower vs. consistent 
        
        [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,5));
        Results_Diff(nBin).outHL_H = H;
        Results_Diff(nBin).outHL_P = P;
        Results_Diff(nBin).outHL_T = STATS.tstat; % outgroup higher vs. lower
        
        [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,6));
        Results_Diff(nBin).outHC_H = H;
        Results_Diff(nBin).outHC_P = P;
        Results_Diff(nBin).outHC_T = STATS.tstat; % outgroup higher vs. consistent 
        
        [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,5), thisBinNeuralRDM(:,6));
        Results_Diff(nBin).outLC_H = H;
        Results_Diff(nBin).outLC_P = P;
        Results_Diff(nBin).outLC_T = STATS.tstat; % outgroup lower vs. consistent 
        
        Results_Diff(nBin).inH_M = mean (thisBinNeuralRDM(:,1));
        Results_Diff(nBin).inH_SE = std (thisBinNeuralRDM(:,1)) / length(thisBinNeuralRDM(:,1));
        Results_Diff(nBin).inL_M = mean (thisBinNeuralRDM(:,2));
        Results_Diff(nBin).inL_SE = std (thisBinNeuralRDM(:,2)) / length(thisBinNeuralRDM(:,2));
        Results_Diff(nBin).inC_M = mean (thisBinNeuralRDM(:,3));
        Results_Diff(nBin).inC_SE = std (thisBinNeuralRDM(:,3)) / length(thisBinNeuralRDM(:,3));
        Results_Diff(nBin).outH_M = mean (thisBinNeuralRDM(:,4));
        Results_Diff(nBin).outH_SE = std (thisBinNeuralRDM(:,4)) / length(thisBinNeuralRDM(:,4));
        Results_Diff(nBin).outL_M = mean (thisBinNeuralRDM(:,5));
        Results_Diff(nBin).outL_SE = std (thisBinNeuralRDM(:,5)) / length(thisBinNeuralRDM(:,5));
        Results_Diff(nBin).outC_M = mean (thisBinNeuralRDM(:,6));
        Results_Diff(nBin).outC_SE = std (thisBinNeuralRDM(:,6)) / length(thisBinNeuralRDM(:,6));
        
    end
    
    
    save(strcat(savName, 'TransRESULTS.mat'),...
        'Results_Diff', 'Results_NeuralBeha',....
        'Results_NeuralConpAll', 'Results_NeuralConpIn');
    
    %% Visualization 
    
    tm = 1: n_bins; 
    rAVG = [Results_NeuralConpAll.R]; 
    seAVG = [Results_NeuralConpAll.SE]; 
    hAVG = [Results_NeuralConpAll.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural All Conceptual Correlation.png'));
    close all;
    
    tm = 1: n_bins; 
    rAVG = [Results_NeuralConpIn.R]; 
    seAVG = [Results_NeuralConpIn.SE]; 
    hAVG = [Results_NeuralConpIn.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural Ingroup Conceptual Correlation.png'));
    close all;
    
    writetable(struct2table(Results_Diff), strcat(savName,'Difference.xlsx'))

    

end


end
