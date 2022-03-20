% This is a script for conducting RSA analysis for Norms Learning

% Author: Danni Chen, Ziqing Yao & Xiaoqing Hu (The University of Hong Kong)
% Update date: 2022-3-02

% the conceptual model include:
%   1) completely no-learning model: #1 self rating
%   2) partially learning (ingroup): #1 self rating & group rating for
%   ingroup member
%   3) partially learning (outgroup): #1 self rating & group rating for
%   outgroup member
%   4) completely learning model: group rating
%   5) final rating model: #2 self rating
%   6) final rating model while controlling self-rating
%   7) final rating model while controlling group-rating

function RSA_cal_stim60_beha(subs, work_dir) 

clear;clc;

%% Subject List: 
if nargin ==0
    
    subs = [1205:1242,1244:1253];
    badsubs   = [1214, 1215, 1238]; 
    subs = setdiff (subs, badsubs);
    
    work_dir = '/home/chendanni/Documents/Norms/analysis/';
    cd (work_dir);
    
    task = "postLearningImplicit";
    erp  = "full_range";
    
    stat_test = "yes"; 
    perm_test = "no";
    ztrans = "yes"; 
    average_erp = "no";
    exclude = "exclude_object_control_morph_faces";
    n_all = 80;
    n_cond = 6; % this is when we are excluding control & morph faces
    bin_list = [1,2,3,4,5,6];
    
end

nSubs = length(subs); %number of subjects


%% Subject Filename(s)

% Load Subject Datafiles 
data_folder = fullfile (work_dir, 'EEGAnalysisResults', 'ERP'); % set directory of data set (default: pwd)
toolbox_folder = '/home/chendanni/MATLAB/toolbox';

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
win = 0.20;
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

%% define path_in & 

path_in  = data_folder;
savName = strcat ('/', 'Results_RSAratingnan_', task ,'_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 
path_out = fullfile(saveLocation, 'Results',strcat('RSAratingnan_', task, '_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
mkdir(path_out)


%% add toolboxes
addpath (genpath(fullfile(work_dir,'MyScripts','additional_functions')));
addpath (fullfile(work_dir,'MyScripts','RSA'));
addpath (fullfile(toolbox_folder,'fieldtrip-20210330'));

%% First step: correlate for each subject in all trials

cd (path_in);

for n = 1:numel(all_subs)
    
    TStart = tic;
    
    sel_sub = all_subs{n};
    sel_sub_folder = fullfile( data_folder, sel_sub );
    
    fprintf (strcat ("Running: Subject ", sel_sub, "\n")); 
    
    
    %% EEG data loading & processing
    % load in the preprocessed data of each participant
    EEG = pop_loadset('filename',char(strcat( sel_sub, '_', task, '_rmbl.set')),'filepath',fullfile(sel_sub_folder,'/'));
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
    
    face_trials = find(arrayfun(@(events) (~ismember(9, events.bini)) & (~ismember(7, events.bini))  & (~ismember(8, events.bini)), events)); 
    % exclude both object trials, control, and morph faces trials
    
    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=face_trials;
    data=ft_timelockanalysis(cfg, data); % keep only experimental face trials, exclude both object trials and control faces trials
    n_trials=numel(face_trials); % number of trials
    n_stim=length (unique ([data.trialinfo.codelabel])); % number of stimuli
    
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
    erpdata_AVG = permute (data.trial, [2, 3, 1]);
    
    %% time-resolved ERP and RSA
    erpdata_AVG_bins_all     = nan (size (erpdata_AVG,1), step, size (erpdata_AVG,3), n_bins); % channel * step * trial * bin
    erpdata_RESHAPE_bins_all = nan (size(erpdata_AVG, 1) * step, size (erpdata_AVG, 3), n_bins); % [channel * step] * trial * bin
    corr_res_all             = nan ( size (erpdata_AVG, 3), size (erpdata_AVG, 3), n_bins); % trial * trial * bin
    corr_res_ztrans_all      = nan ( size (erpdata_AVG, 3), size (erpdata_AVG, 3), n_bins); % trial * trial * bin
    
    for nbin = 1:n_bins
        
        erpdata_AVG_thisbin = erpdata_AVG(:, ind_t1(nbin):ind_t2(nbin), :);
        erpdata_RESHAPE_thisbin = reshape (erpdata_AVG_thisbin, size(erpdata_AVG_thisbin, 1) * size(erpdata_AVG_thisbin, 2), size (erpdata_AVG_thisbin, 3)); % reshape the ERP into (nChannels * nTp) * nTrials
        
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
    
    %% reshape the trial-based correlation matrix to stimuli-based correlation matrix
    corr_res_stim = nan (n_stim, n_stim, n_bins);
    corr_res_ztrans_stim = nan (n_stim, n_stim, n_bins); 
    trialinfo = data.trialinfo;
    triallabel = [string(trialinfo.codelabel)];
    
    stim_list = [];
    for tmp_i = 1:n_stim
        stim_list = [stim_list; sscanf(triallabel(tmp_i),"s%d")];
    end
    stim_list = sort (stim_list)';
    
    for i_stim = 1:n_all 
        for j_stim = 1:n_all
            
            i_trials = find (strcmp (strcat('s', num2str(i_stim)) , triallabel))';
            j_trials = find (strcmp (strcat('s', num2str(j_stim)) , triallabel))';
            corr_ij = corr_res_all (i_trials, j_trials, :);
            corr_z_ij = corr_res_ztrans_all (i_trials, j_trials, :);
            
            for i_bin = 1:n_bins
                i_corr_ij = corr_ij(:,:,i_bin);
                i_corr_z_ij = corr_z_ij (:,:,i_bin);
                if (i_stim == j_stim)
                    i_corr_ij = tril (i_corr_ij, -1);
                    i_corr_ij(i_corr_ij == 0) = nan;
                    i_corr_z_ij = tril (i_corr_z_ij, -1);
                    i_corr_z_ij(i_corr_z_ij == 0) = nan;
                    corr_res_stim (i_stim, j_stim, i_bin) = nanmean (nanmean (i_corr_ij));
                    corr_res_ztrans_stim (i_stim, j_stim, i_bin) = nanmean (nanmean (i_corr_z_ij));
                else
                    corr_res_stim (i_stim, j_stim, i_bin) = nanmean (nanmean (i_corr_ij));
                    corr_res_ztrans_stim (i_stim, j_stim, i_bin) = nanmean (nanmean (i_corr_z_ij));
                end
                clear  i_corr_ij i_corr_z_ij; 
            end
            
            clear i_trials j_trials corr_ij corr_z_ij;

        end
    end
    
    %% save data 
    save(fullfile(path_out, strcat(savName, all_subs{n},'.mat')),...
        'corr_res_stim', 'corr_res_ztrans_stim','trialinfo');
    
    clear corr_trials corr_cue_enc_trial trialinfo corr_cue_enc items corr_res cfg code_nr codelabel codeseq data EEG erpdata erpdata2 erpdata3 events face_trials object_trials idx i this_code this_data this erp this_trials

    toc(TStart);
    
end % end of each subject 

%% ----------------------------------------------------------------------
%% Second step: Correlate with Conceptual and Behavioral RDM

this_rsa_folder = path_out;
cd( this_rsa_folder );

allsubject_neuralrdm_allbins              = zeros (n_stim, n_stim, n_bins, numel(subs));
allsubject_nolearning_allbins             = zeros (n_stim, n_stim, n_bins, numel(subs)); % #1 self rating
allsubject_iglearning_allbins             = zeros (n_stim, n_stim, n_bins, numel(subs));
allsubject_oglearning_allbins            = zeros (n_stim, n_stim, n_bins, numel(subs));
allsubject_learnt_allbins                 = zeros (n_stim, n_stim, n_bins, numel(subs)); % group rating
allsubject_final_allbins                  = zeros (n_stim, n_stim, n_bins, numel(subs)); % #2 self rating

allsubject = struct;

savName = strcat ('Results_RSAratingnan_', task ,'_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 

selfrating1 = readmatrix (fullfile(work_dir, 'BehaData', 'PreLearningExplicitBehaData2021-02-23.csv'));
selfrating2 = readmatrix (fullfile(work_dir, 'BehaData', 'PostLearningExplicitBehaData2021-02-23.csv'));
grouprating = readmatrix (fullfile(work_dir, 'BehaData', 'LearningBehaData2021-06-23.csv'));

selfrating1 = selfrating1 (:,[8,30,81]); % figureID (same as triggerID); self rating; subjectID
selfrating2 = selfrating2 (:,[8,30,81]); % figureID (same as triggerID); self rating; subjectID 
grouprating = grouprating (:,[8,12,15,25]); % group rating; subjectID; triggerID; block
grouprating = grouprating (grouprating(:,4)==0, :); 

cnt = 1;

for sub_nr = 1:numel(subs)
    
    thissub = num2str(subs(sub_nr));
    
    cd(this_rsa_folder);
    load (strcat(savName, thissub,'.mat')); % corr_res_stim, corr_res_ztrans_stim, trialinfo
    clear corr_res_stim; % clear raw correlation coefficient matrix in case of mistakes
    
    corr_res_all = 1 - corr_res_ztrans_stim; % simialrity to RDM
    
    sel_sub_folder = fullfile( data_folder, thissub );
    EEG = pop_loadset('filename',char(strcat( thissub, '_', task, '_rmbl.set')),'filepath',fullfile(sel_sub_folder,'/'));
   
    thisbini = [];
    thiscode = [];
    for i = 1: size(EEG.event, 2)
        thisbini(i) = EEG.event(i).bini;
        thislabel = EEG.event(i).codelabel;
        thiscode(i) = str2double(thislabel(2:end));
    end
    uniquecode = unique(thiscode, 'first');
    uniquecode = uniquecode(~isnan(uniquecode));
    thisbini = thisbini (uniquecode);
    thiscode = thiscode (uniquecode);
    thisall  = [thisbini;thiscode]'; % Code and Corresponding Bin
    
    controlFigure = sort(thisall(thisall(:,1)==7,2))'
    controlthres  = max(controlFigure);

    thisbinlist = zeros(n_cond, 10);
    for i = 1:n_cond % n_cond bins
        x = ceil((min(thiscode(thisbini == bin_list(i))))/10);
        x1 = ((x-1)*10+1):x*10;
        thisbinlist(i,:) = x1;
    end % check each bin includes which figures
    newseq = reshape(thisbinlist',1,n_cond*10);
    
    
    % 1) completely no-learning model with #1 self rating
    nolearning_rdm = nan (n_cond*10, n_cond*10);
    thisbeha = selfrating1(selfrating1(:,3) == subs(sub_nr),:);
    for i = 1: n_cond*10
        for j = 1: n_cond*10
            if j>=i 
                nolearning_rdm (i, j) = nan;
            elseif (ceil(i/10) == ceil (j/10))
                nolearning_rdm (i, j) = nan;
            else
                rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
                rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
                nolearning_rdm(i,j) = abs(rating_i - rating_j);
            end
        end
    end 

    % 2) partially learning (ingroup) #1 self rating & group rating for
    % ingroup conditions
    iglearning_rdm = nan (n_cond*10, n_cond*10);
    thisbeha  = selfrating1(selfrating1(:,3) == subs(sub_nr),:);
    thisgroup = grouprating (grouprating(:,2) == subs(sub_nr),:);
    for i = 1: n_cond*10
        for j = 1: n_cond*10
            % for i 
            if ismember(i, [1:30])
                rating_i = thisgroup(find(thisgroup(:,3)==newseq(i)), 1);
            elseif ismember(i, [31:60])
                rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
            end
            % for j 
            if ismember(j, [1:30])
                rating_j = thisgroup(find(thisgroup(:,3)==newseq(j)), 1);
            elseif ismember (j, [31:60]);
                rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
            end            
            % construct rdm
            if j>=i 
                iglearning_rdm (i, j) = nan;
            elseif (ceil(i/10) == ceil (j/10))
                iglearning_rdm (i, j) = nan;
            elseif ismember(i, [31:60]) && ismember(j, [1:30])
                iglearning_rdm (i, j) = nan;
            elseif ismember(i, [1:30]) && ismember(j, [31:60])
                iglearning_rdm (i, j) = nan;            
            else
                iglearning_rdm(i,j) = abs(rating_i - rating_j);
            end
        end
    end 
    
    % 3) partially learning (outgroup) #1 self rating & group rating for
    % outgroup conditions
    oglearning_rdm = nan (n_cond*10, n_cond*10);
    thisbeha  = selfrating1(selfrating1(:,3) == subs(sub_nr),:);
    thisgroup = grouprating (grouprating(:,2) == subs(sub_nr),:);
    for i = 1: n_cond*10
        for j = 1: n_cond*10
            % for i 
            if ismember(i, [31:60])
                rating_i = thisgroup(find(thisgroup(:,3)==newseq(i)), 1);
            elseif ismember(i, [1:30])
                rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
            end
            % for j 
            if ismember(j, [31:60])
                rating_j = thisgroup(find(thisgroup(:,3)==newseq(j)), 1);
            elseif ismember (j, [1:30]);
                rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
            end            
            % construct rdm
            if j>=i 
                oglearning_rdm (i, j) = nan;
            elseif (ceil(i/10) == ceil (j/10))
                oglearning_rdm (i, j) = nan;
            elseif ismember(i, [31:60]) && ismember(j, [1:30])
                oglearning_rdm (i, j) = nan;
            elseif ismember(i, [1:30]) && ismember(j, [31:60])
                oglearning_rdm (i, j) = nan;
            else
                oglearning_rdm(i,j) = abs(rating_i - rating_j);
            end
        end
    end 
    
    % 4) completely learnt model with group rating
    learnt_rdm = nan (n_cond*10, n_cond*10);
    thisgroup = grouprating (grouprating(:,2) == subs(sub_nr),:);
    for i = 1: n_cond*10
        for j = 1: n_cond*10
            if j>=i 
                learnt_rdm (i, j) = nan;
            elseif (ceil(i/10) == ceil (j/10))
                learnt_rdm (i, j) = nan;
            else
                rating_i = thisgroup(find(thisgroup(:,3)==newseq(i)), 1);
                rating_j = thisgroup(find(thisgroup(:,3)==newseq(j)), 1);
                learnt_rdm(i,j) = abs(rating_i - rating_j);
            end
        end
    end 
    
    % 5) final model with #2 self rating
    final_rdm = nan (n_cond*10, n_cond*10);
    thisbeha = selfrating2(selfrating2(:,3) == subs(sub_nr),:);
    for i = 1: n_cond*10
        for j = 1: n_cond*10
            if j>=i 
                final_rdm (i, j) = nan;
            elseif (ceil(i/10) == ceil (j/10))
                final_rdm (i, j) = nan;
            else
                rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
                rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
                final_rdm(i,j) = abs(rating_i - rating_j);
            end
        end
    end 

    
    for nbin = 1:n_bins
        
        corr_res = corr_res_all(:,:,nbin);
   
        neural_rdm = zeros (n_cond*10, n_cond*10);
        for i = 1: n_cond*10
            for j = 1: n_cond*10
                neural_rdm(i,j) = corr_res (newseq(i),newseq(j));
            end
        end % neural_rdm, rearrange corr_res by bin
        
        allsubject_neuralrdm_allbins (:, :, nbin, sub_nr)  = neural_rdm;
        allsubject_nolearning_allbins (:, :, nbin, sub_nr) = nolearning_rdm;
        allsubject_iglearning_allbins (:, :, nbin, sub_nr) = iglearning_rdm;
        allsubject_oglearning_allbins (:, :, nbin, sub_nr) = oglearning_rdm;
        allsubject_learnt_allbins (:, :, nbin, sub_nr)     = learnt_rdm;
        allsubject_final_allbins (:, :, nbin, sub_nr)      = final_rdm;
     
    end % end of bins
        
end % end of subjects

RSA_heatmap_plot (squeeze(mean (allsubject_neuralrdm_allbins(:,:,1,:), 4)), 'Neural-RDM');
RSA_heatmap_plot (squeeze(mean (allsubject_nolearning_allbins(:,:,1,:), 4)), 'Nolearning-Conceptual-RDM');
RSA_heatmap_plot (squeeze(mean (allsubject_iglearning_allbins(:,:,1,:), 4)), 'Partiallearning-Ingroup-RDM');
RSA_heatmap_plot (squeeze(mean (allsubject_oglearning_allbins(:,:,1,:), 4)), 'Partiallearning-Outgroup-RDM');
RSA_heatmap_plot (squeeze(mean (allsubject_learnt_allbins(:,:,1,:), 4)), 'Learnt-Conceptual-RDM');
RSA_heatmap_plot (squeeze(mean (allsubject_final_allbins(:,:,1,:), 4)), 'Final-Beha-RDM');


save(strcat(savName, 'trans.mat'), 'allsubject_neuralrdm_allbins',...
    'allsubject_nolearning_allbins','allsubject_iglearning_allbins',...
    'allsubject_oglearning_allbins','allsubject_learnt_allbins',...
    'allsubject_learnt_allbins', 'allsubject_final_allbins');

%% ----------------------------------------------------------------------

%% Conduct statical analysis

if strcmp (stat_test, 'yes')
    
    %% Data Prep. 
    
    load (strcat(savName, 'trans.mat'));
    
%     allsubject_neuralrdm_allbins (:, :, nbin, sub_nr)  nItems x nItems x nBins x nSubject
%     allsubject_nolearning_allbins (:, :, nbin, sub_nr) nItems x nItems x nBins x nSubject
%     allsubject_iglearning_allbins (:, :, nbin, sub_nr) nItems x nItems x nBins x nSubject
%     allsubject_oglearning_allbins (:, :, nbin, sub_nr) nItems x nItems x nBins x nSubject
%     allsubject_learnt_allbins (:, :, nbin, sub_nr)     nItems x nItems x nBins x nSubject
%     allsubject_final_allbins (:, :, nbin, sub_nr)      nItems x nItems x nBins x nSubject
    
    n_bins = size (allsubject_neuralrdm_allbins, 3);
    
    allsubject_neuralrdm_allbin_AVG = nan (size (allsubject_neuralrdm_allbins, 4), n_cond, n_bins);
    % nSubject x 6 nConds) x nBin
    allsubject_neuralnolearning_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuraliglearning_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuraloglearning_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neurallearnt_R     = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralfinal_R      = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralfinal_cSelf_R      = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins); % c: control self-rating
    allsubject_neuralfinal_cGroup_R     = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins); % c: control group-rating
    
    allsubject_neuralnolearning_P = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuraliglearning_P = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuraloglearning_P = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neurallearnt_P     = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralfinal_P      = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
    allsubject_neuralfinal_cSelf_P      = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins); % c: control self-rating
    allsubject_neuralfinal_cGroup_P     = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins); % c: control group-rating
    % nSubject x Bin
    
    
    for nSub = 1: size (allsubject_neuralrdm_allbins, 4)
        
        for nBin = 1: size (allsubject_neuralrdm_allbins, 3)
            
            neuralrdm      = allsubject_neuralrdm_allbins (:, :, nBin, nSub);
            nolearning_rdm = allsubject_nolearning_allbins (:, :, nBin, nSub);
            iglearning_rdm = allsubject_iglearning_allbins (:, :, nBin, nSub);
            oglearning_rdm = allsubject_oglearning_allbins (:, :, nBin, nSub);
            learnt_rdm     = allsubject_learnt_allbins (:, :, nBin, nSub);
            final_rdm      = allsubject_final_allbins (:, :, nBin, nSub);
            
            %% Correlation between neural rdm and 1) completely no-learning model #1 - self rating
            behardm = matrix2vector (nolearning_rdm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            neuralrdm_VECTOR(isnan(behardm)) = []; % remove same condition 
            behardm(isnan(behardm)) = [];
            [neuralbeha_corr, neuralbeha_corr_p]  = corr(neuralrdm_VECTOR', behardm', 'type','Spearman');
            allsubject_neuralnolearning_R(nSub, nBin)  = neuralbeha_corr;
            allsubject_neuralnolearning_P(nSub, nBin)   = neuralbeha_corr_p;
            
            %% Correlation between neural rdm and 2) partially learning model (ingroup)
            behardm = matrix2vector (iglearning_rdm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            neuralrdm_VECTOR(isnan(behardm)) = []; % remove same condition 
            behardm(isnan(behardm)) = [];
            [neuralbeha_corr, neuralbeha_corr_p]  = corr(neuralrdm_VECTOR', behardm', 'type','Spearman');
            allsubject_neuraliglearning_R(nSub, nBin)  = neuralbeha_corr;
            allsubject_neuraliglearning_P(nSub, nBin)   = neuralbeha_corr_p;
            
            %% Correlation between neural rdm and 3) partially learning model (outgroup)
            behardm = matrix2vector (oglearning_rdm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            neuralrdm_VECTOR(isnan(behardm)) = []; % remove same condition 
            behardm(isnan(behardm)) = [];
            [neuralbeha_corr, neuralbeha_corr_p]  = corr(neuralrdm_VECTOR', behardm', 'type','Spearman');
            allsubject_neuraloglearning_R(nSub, nBin)  = neuralbeha_corr;
            allsubject_neuraloglearning_P(nSub, nBin)   = neuralbeha_corr_p;        
        
            %% Correlation between neural rdm and 4) completely learnt model 
            behardm = matrix2vector (learnt_rdm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            neuralrdm_VECTOR(isnan(behardm)) = []; % remove same condition 
            behardm(isnan(behardm)) = [];
            [neuralbeha_corr, neuralbeha_corr_p]  = corr(neuralrdm_VECTOR', behardm', 'type','Spearman');
            allsubject_neurallearnt_R(nSub, nBin)  = neuralbeha_corr;
            allsubject_neurallearnt_P(nSub, nBin)   = neuralbeha_corr_p;     
  
            %% Correlation between neural rdm and 5) final data #2 self rating
            behardm = matrix2vector (final_rdm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            neuralrdm_VECTOR(isnan(behardm)) = []; % remove same condition 
            behardm(isnan(behardm)) = [];
            [neuralbeha_corr, neuralbeha_corr_p]  = corr(neuralrdm_VECTOR', behardm', 'type','Spearman');
            allsubject_neuralfinal_R(nSub, nBin)  = neuralbeha_corr;
            allsubject_neuralfinal_P(nSub, nBin)  = neuralbeha_corr_p;     
            
            %% Correlation between neural rdm and final data #2 self rating, while controlling #1 selfrating
            behardm = matrix2vector (final_rdm, "lower");
            contrdm = matrix2vector (nolearning_rdm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            neuralrdm_VECTOR(isnan(behardm)) = []; % remove same condition 
            contrdm(isnan(behardm)) = [];
            behardm(isnan(behardm)) = [];
            [neuralbeha_corr, neuralbeha_corr_p]  = partialcorr (neuralrdm_VECTOR', behardm', contrdm', 'rows','complete', 'type','Spearman');
            allsubject_neuralfinal_cSelf_R(nSub, nBin)  = neuralbeha_corr;
            allsubject_neuralfinal_cSelf_P(nSub, nBin)  = neuralbeha_corr_p;     
        
            %% Correlation between neural rdm and final data #2 self rating, while controlling grouprating
            behardm = matrix2vector (final_rdm, "lower");
            contrdm = matrix2vector (learnt_rdm, "lower");
            neuralrdm_VECTOR = matrix2vector (neuralrdm, "lower");
            neuralrdm_VECTOR(isnan(behardm)) = []; % remove same condition 
            contrdm(isnan(behardm)) = [];
            behardm(isnan(behardm)) = [];
            [neuralbeha_corr, neuralbeha_corr_p]  = partialcorr (neuralrdm_VECTOR', behardm', contrdm', 'rows','complete', 'type','Spearman');
            allsubject_neuralfinal_cGroup_R(nSub, nBin)  = neuralbeha_corr;
            allsubject_neuralfinal_cGroup_P(nSub, nBin)  = neuralbeha_corr_p; 
            
        end
        
    end
    
    
    Results_NeuralNolearning = struct;
    Results_NeuralIglearning = struct;
    Results_NeuralOglearning = struct;
    Results_NeuralLearnt     = struct;
    Results_NeuralFinal      = struct;
    Results_NeuralFinal_cSelf      = struct;
    Results_NeuralFinal_cGroup     = struct;
    
    %% Significant Analysis 
    for nBin = 1:n_bins
        
        %% Correlation between neural rdm and 1) nolearning conceptual RDM 
        [Results_NeuralNolearning(nBin).H, Results_NeuralNolearning(nBin).P] = ttest (allsubject_neuralnolearning_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralNolearning(nBin).R = mean (allsubject_neuralnolearning_R(:,nBin));
        Results_NeuralNolearning(nBin).SE = std (allsubject_neuralnolearning_R(:,nBin)) / sqrt(length (allsubject_neuralnolearning_R(:,nBin)));
        
        %% Correlation between neural rdm and 2) iglearning conceptual RDM 
        [Results_NeuralIglearning(nBin).H, Results_NeuralIglearning(nBin).P] = ttest (allsubject_neuraliglearning_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralIglearning(nBin).R = mean (allsubject_neuraliglearning_R(:,nBin));
        Results_NeuralIglearning(nBin).SE = std (allsubject_neuraliglearning_R(:,nBin)) / sqrt(length (allsubject_neuraliglearning_R(:,nBin)));
                
        %% Correlation between neural rdm and 3) oglearning conceptual RDM 
        [Results_NeuralOglearning(nBin).H, Results_NeuralOglearning(nBin).P] = ttest (allsubject_neuraloglearning_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralOglearning(nBin).R = mean (allsubject_neuraloglearning_R(:,nBin));
        Results_NeuralOglearning(nBin).SE = std (allsubject_neuraloglearning_R(:,nBin)) / sqrt(length (allsubject_neuraloglearning_R(:,nBin)));

        %% Correlation between neural rdm and 4) learnt conceptual RDM 
        [Results_NeuralLearnt(nBin).H, Results_NeuralLearnt(nBin).P] = ttest (allsubject_neurallearnt_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralLearnt(nBin).R = mean (allsubject_neurallearnt_R(:,nBin));
        Results_NeuralLearnt(nBin).SE = std (allsubject_neurallearnt_R(:,nBin)) / sqrt(length (allsubject_neurallearnt_R(:,nBin)));

        %% Correlation between neural rdm and 5) final conceptual RDM 
        [Results_NeuralFinal(nBin).H, Results_NeuralFinal(nBin).P] = ttest (allsubject_neuralfinal_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralFinal(nBin).R = mean (allsubject_neuralfinal_R(:,nBin));
        Results_NeuralFinal(nBin).SE = std (allsubject_neuralfinal_R(:,nBin)) / sqrt(length (allsubject_neuralfinal_R(:,nBin)));
         
        %% Correlation between neural rdm and final conceptual RDM and control 6) nolearning RDM
        [Results_NeuralFinal_cSelf(nBin).H, Results_NeuralFinal_cSelf(nBin).P] = ttest (allsubject_neuralfinal_cSelf_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralFinal_cSelf(nBin).R = mean (allsubject_neuralfinal_cSelf_R(:,nBin));
        Results_NeuralFinal_cSelf(nBin).SE = std (allsubject_neuralfinal_cSelf_R(:,nBin)) / sqrt(length (allsubject_neuralfinal_cSelf_R(:,nBin)));

        %% Correlation between neural rdm and final conceptual RDM and control 7) learnt RDM
        [Results_NeuralFinal_cGroup(nBin).H, Results_NeuralFinal_cGroup(nBin).P] = ttest (allsubject_neuralfinal_cGroup_R(:,nBin), 0, 'tail', 'right');
        Results_NeuralFinal_cGroup(nBin).R = mean (allsubject_neuralfinal_cGroup_R(:,nBin));
        Results_NeuralFinal_cGroup(nBin).SE = std (allsubject_neuralfinal_cGroup_R(:,nBin)) / sqrt(length (allsubject_neuralfinal_cGroup_R(:,nBin)));
                   
    end
    
    
    save(strcat(savName, 'TransRESULTS.mat'),...
        'Results_NeuralNolearning', 'Results_NeuralIglearning',....
        'Results_NeuralOglearning', 'Results_NeuralLearnt','Results_NeuralFinal',...
        'Results_NeuralFinal_cSelf','Results_NeuralFinal_cGroup');
    
    %% Visualization 
    
    tm = 1: n_bins; 
    rAVG  = [Results_NeuralNolearning.R]; 
    seAVG = [Results_NeuralNolearning.SE]; 
    hAVG  = [Results_NeuralNolearning.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural-Nolearning.png'));
    close all;
    
    tm = 1: n_bins; 
    rAVG  = [Results_NeuralIglearning.R]; 
    seAVG = [Results_NeuralIglearning.SE]; 
    hAVG  = [Results_NeuralIglearning.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural-Ingroup-learning.png'));
    close all;
    
    tm = 1: n_bins; 
    rAVG  = [Results_NeuralOglearning.R]; 
    seAVG = [Results_NeuralOglearning.SE]; 
    hAVG  = [Results_NeuralOglearning.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural-Outgroup-learning.png'));
    close all;
    
    tm = 1: n_bins; 
    rAVG  = [Results_NeuralLearnt.R]; 
    seAVG = [Results_NeuralLearnt.SE]; 
    hAVG  = [Results_NeuralLearnt.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural-Learnt.png'));
    close all;
    
    tm = 1: n_bins; 
    rAVG  = [Results_NeuralFinal.R]; 
    seAVG = [Results_NeuralFinal.SE]; 
    hAVG  = [Results_NeuralFinal.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural-Final.png'));
    close all;
    
    tm = 1: n_bins; 
    rAVG  = [Results_NeuralFinal_cSelf.R]; 
    seAVG = [Results_NeuralFinal_cSelf.SE]; 
    hAVG  = [Results_NeuralFinal_cSelf.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural-Final-controlSelfRating#1.png'));
    close all;
    
    tm = 1: n_bins; 
    rAVG  = [Results_NeuralFinal_cGroup.R]; 
    seAVG = [Results_NeuralFinal_cGroup.SE]; 
    hAVG  = [Results_NeuralFinal_cGroup.H];
    chancelvl = 0;
    figure(1)
    cl=colormap(parula(50));
    plot(tm, rAVG); %plot
    text (tm(hAVG==1), ones(length(hAVG(hAVG==1)),1)'*(0), '*', 'Color', 'red', 'FontSize', 14);
    boundedline(tm,rAVG,seAVG,'cmap',cl(42,:),'alpha','transparency',0.35)
    line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
    saveas(gcf,strcat(savName,'Neural-Final-controlGroupRating.png'));
    close all;
    
end


end
