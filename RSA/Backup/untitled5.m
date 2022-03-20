% This is a script for conducting Statistical analysis with RSA results 
% Author: Danni Chen, Ziqing Yao & Xiaoqing Hu (The University of Hong Kong)
% Update date: 11/28/2021

%% define basic information 
subs = [1205:1242,1244:1253];
badsubs   = [1201:1204, 1220, 1229, 1238, 1239, 1249]; % 1201 - 1204: Local, 1229, 1238, 1239,1249 (reject epochs > 100)
subs = setdiff (subs, badsubs);
nSubs = length(subs); %number of subjects

work_dir = '/home/chendanni/Documents/Norms/analysis/';
cd (work_dir);

erp  = "full_range";

stat_test = "yes"; 
perm_test = "no";
ztrans = "yes"; 
average_erp = "no";
permutation = "no";

% Load Subject Datafiles 
data_folder = fullfile (work_dir, 'EEGAnalysisResults', 'Preprocessing'); % set directory of data set (default: pwd)
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

%% set up permutation parameter 
rand = 1000;
p_cluster = 0.05;
p_first   = 0.05;
contrast='ttest';    
load(fullfile(work_dir,'MyScripts','additional_functions','jet_grey2.mat'));

%% define path_in & path_out

path_in  = data_folder;
savName = strcat ('Results_RSA_diffImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 

path_pre  = fullfile(saveLocation, 'Results',strcat('RSA_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
path_post = fullfile(saveLocation, 'Results',strcat('RSA_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
path_diff = fullfile(saveLocation, 'Stats',strcat('RSA_diffImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
if ~exist (path_diff)
    mkdir(path_diff)
end

%% difference between post-learning and pre-learning

RDM_pre = load ( strcat (path_pre, '/', 'Results_RSA_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_trans.mat') );
RDM_post = load ( strcat (path_post, '/', 'Results_RSA_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_trans.mat') );
% allsubject_behardm_allbins: nItems x nItems x nBins x nSubject
% allsubject_conceprdm_conflict_allbins: nItems x nItems x nBins x nSubject
% allsubject_conceprdm_inconflict_allbins: nItems x nItems x nBins x nSubject
% allsubject_neuralrdm_allbin: nItems x nItems x nBins x nSubject
RDM_diff.allsubject_neuralrdm_allbins = RDM_post.allsubject_neuralrdm_allbins - RDM_pre.allsubject_neuralrdm_allbins;
RDM_diff.allsubject_behardm_allbins   = RDM_post.allsubject_behardm_allbins - RDM_pre.allsubject_behardm_allbins;
RDM_diff.allsubject_coceprdm_conflict_allbins    = RDM_post.allsubject_coceprdm_conflict_allbins;
RDM_diff.allsubject_coceprdm_inconflict_allbins  = RDM_post.allsubject_coceprdm_inconflict_allbins;
RDM_diff.allsubject_coceprdm_outconflict_allbins = RDM_post.allsubject_coceprdm_outconflict_allbins;
save (strcat (path_diff, '/', 'Results_RSA_diffImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_trans.mat'),...
    'RDM_pre','RDM_post','RDM_diff');

%% comparing higher vs. lower, lower vs. consistent, higher vs. consistent

n_bins = size (RDM_diff.allsubject_neuralrdm_allbins, 3);

allsubject_neuralrdm_allbin_Pre_AVG  = nan (size (RDM_diff.allsubject_neuralrdm_allbins, 4), 8, n_bins);
allsubject_neuralrdm_allbin_Post_AVG = nan (size (RDM_diff.allsubject_neuralrdm_allbins, 4), 8, n_bins);
allsubject_neuralrdm_allbin_Diff_AVG = nan (size (RDM_diff.allsubject_neuralrdm_allbins, 4), 8, n_bins);

for nSub = 1: size (RDM_pre.allsubject_neuralrdm_allbins, 4)

    for nBin = 1: size (RDM_pre.allsubject_neuralrdm_allbins, 3)
        
        %% pre-learning 
        neuralrdm   = RDM_pre.allsubject_neuralrdm_allbins (:, :, nBin, nSub);
        for nCond = 1: 6
            allsubject_neuralrdm_allbin_Pre_AVG (nSub, nCond, nBin) = mean (neuralrdm(((nCond-1)*10 +1):nCond*10, 71:80) , 'all');
        end
        allsubject_neuralrdm_allbin_Pre_AVG (nSub, 7, nBin) = nSub;
        allsubject_neuralrdm_allbin_Pre_AVG (nSub, 8, nBin) = nBin;

        %% post-learning 
        neuralrdm   = RDM_post.allsubject_neuralrdm_allbins (:, :, nBin, nSub);
        for nCond = 1: 6
            allsubject_neuralrdm_allbin_Post_AVG (nSub, nCond, nBin) = mean (neuralrdm(((nCond-1)*10 +1):nCond*10, 71:80) , 'all');
        end
        allsubject_neuralrdm_allbin_Post_AVG (nSub, 7, nBin) = nSub;
        allsubject_neuralrdm_allbin_Post_AVG (nSub, 8, nBin) = nBin;
        
        %% difference
        neuralrdm   = RDM_diff.allsubject_neuralrdm_allbins (:, :, nBin, nSub);
        for nCond = 1: 6
            allsubject_neuralrdm_allbin_Diff_AVG (nSub, nCond, nBin) = mean (neuralrdm(((nCond-1)*10 +1):nCond*10, 71:80) , 'all');
        end
        allsubject_neuralrdm_allbin_Diff_AVG (nSub, 7, nBin) = nSub;
        allsubject_neuralrdm_allbin_Diff_AVG (nSub, 8, nBin) = nBin;
    end

end

Results_Pre  = struct;
Results_Post  = struct;
Results_Diff  = struct;

%% Significant Analysis 

contrast_pairs = [1 2; 1 3; 3 2; 4 5; 4 6; 6 1];
contrast_names = {'in_HL', 'in_HC', 'in_CL','in_HL', 'out_HC', 'out_CL'}

% pre-learning
Results_TMP = struct; 
for iContrast = 1:6
    AVG = 
    
end


inHL_AVG = squeeze (allsubject_neuralrdm_allbin_Pre_AVG(:,1,:) - allsubject_neuralrdm_allbin_Pre_AVG(:,2,:));
[inHL_Results.H, inHL_Results.P, inHL_Results.CI, inHL_Results.STATS] = ttest (inHL_AVG);
inHC_AVG = squeeze (allsubject_neuralrdm_allbin_Pre_AVG(:,1,:) - allsubject_neuralrdm_allbin_Pre_AVG(:,3,:));
[inHC_Results.H, inHC_Results.P, inHC_Results.CI, inHC_Results.STATS] = ttest (inHC_AVG);
inCL_AVG = squeeze (allsubject_neuralrdm_allbin_Pre_AVG(:,3,:) - allsubject_neuralrdm_allbin_Pre_AVG(:,2,:));
[inCL_Results.H, inCL_Results.P, inCL_Results.CI, inCL_Results.STATS] = ttest (inCL_AVG);

for nBin = 1:n_bins

    thisBinNeuralRDM = allsubject_neuralrdm_allbin_Pre_AVG(:,:,nBin);
    
    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inHL_H = H;
    Results_TMP(nBin).inHL_P = P;
    Results_TMP(nBin).inHL_T = STATS.tstat; % ingroup higher vs. lower

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inHC_H = H;
    Results_TMP(nBin).inHC_P = P;
    Results_TMP(nBin).inHC_T = STATS.tstat; % ingroup higher vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,2), thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inLC_H = H;
    Results_TMP(nBin).inLC_P = P;
    Results_TMP(nBin).inLC_T = STATS.tstat; % ingroup lower vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outHL_H = H;
    Results_TMP(nBin).outHL_P = P;
    Results_TMP(nBin).outHL_T = STATS.tstat; % outgroup higher vs. lower

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outHC_H = H;
    Results_TMP(nBin).outHC_P = P;
    Results_TMP(nBin).outHC_T = STATS.tstat; % outgroup higher vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,5), thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outLC_H = H;
    Results_TMP(nBin).outLC_P = P;
    Results_TMP(nBin).outLC_T = STATS.tstat; % outgroup lower vs. consistent 

    Results_TMP(nBin).inH_M = mean (thisBinNeuralRDM(:,1));
    Results_TMP(nBin).inH_SE = std (thisBinNeuralRDM(:,1)) / length(thisBinNeuralRDM(:,1));
    Results_TMP(nBin).inL_M = mean (thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inL_SE = std (thisBinNeuralRDM(:,2)) / length(thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inC_M = mean (thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inC_SE = std (thisBinNeuralRDM(:,3)) / length(thisBinNeuralRDM(:,3));
    Results_TMP(nBin).outH_M = mean (thisBinNeuralRDM(:,4));
    Results_TMP(nBin).outH_SE = std (thisBinNeuralRDM(:,4)) / length(thisBinNeuralRDM(:,4));
    Results_TMP(nBin).outL_M = mean (thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outL_SE = std (thisBinNeuralRDM(:,5)) / length(thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outC_M = mean (thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outC_SE = std (thisBinNeuralRDM(:,6)) / length(thisBinNeuralRDM(:,6));

    
    
end
Results_Pre = Results_TMP;

    %% Permutation Test
    if ( strcmp (perm_test, 'yes') )
        
    

    end

% post-learning
Results_TMP = struct; 
for nBin = 1:n_bins

    thisBinNeuralRDM = allsubject_neuralrdm_allbin_Post_AVG(:,:,nBin);
    
    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inHL_H = H;
    Results_TMP(nBin).inHL_P = P;
    Results_TMP(nBin).inHL_T = STATS.tstat; % ingroup higher vs. lower

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inHC_H = H;
    Results_TMP(nBin).inHC_P = P;
    Results_TMP(nBin).inHC_T = STATS.tstat; % ingroup higher vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,2), thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inLC_H = H;
    Results_TMP(nBin).inLC_P = P;
    Results_TMP(nBin).inLC_T = STATS.tstat; % ingroup lower vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outHL_H = H;
    Results_TMP(nBin).outHL_P = P;
    Results_TMP(nBin).outHL_T = STATS.tstat; % outgroup higher vs. lower

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outHC_H = H;
    Results_TMP(nBin).outHC_P = P;
    Results_TMP(nBin).outHC_T = STATS.tstat; % outgroup higher vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,5), thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outLC_H = H;
    Results_TMP(nBin).outLC_P = P;
    Results_TMP(nBin).outLC_T = STATS.tstat; % outgroup lower vs. consistent 

    Results_TMP(nBin).inH_M = mean (thisBinNeuralRDM(:,1));
    Results_TMP(nBin).inH_SE = std (thisBinNeuralRDM(:,1)) / length(thisBinNeuralRDM(:,1));
    Results_TMP(nBin).inL_M = mean (thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inL_SE = std (thisBinNeuralRDM(:,2)) / length(thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inC_M = mean (thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inC_SE = std (thisBinNeuralRDM(:,3)) / length(thisBinNeuralRDM(:,3));
    Results_TMP(nBin).outH_M = mean (thisBinNeuralRDM(:,4));
    Results_TMP(nBin).outH_SE = std (thisBinNeuralRDM(:,4)) / length(thisBinNeuralRDM(:,4));
    Results_TMP(nBin).outL_M = mean (thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outL_SE = std (thisBinNeuralRDM(:,5)) / length(thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outC_M = mean (thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outC_SE = std (thisBinNeuralRDM(:,6)) / length(thisBinNeuralRDM(:,6));
   
end
Results_Post = Results_TMP;


% difference
Results_TMP = struct; 
for nBin = 1:n_bins

    thisBinNeuralRDM = allsubject_neuralrdm_allbin_Diff_AVG(:,:,nBin);
    
    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inHL_H = H;
    Results_TMP(nBin).inHL_P = P;
    Results_TMP(nBin).inHL_T = STATS.tstat; % ingroup higher vs. lower

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,1), thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inHC_H = H;
    Results_TMP(nBin).inHC_P = P;
    Results_TMP(nBin).inHC_T = STATS.tstat; % ingroup higher vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,2), thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inLC_H = H;
    Results_TMP(nBin).inLC_P = P;
    Results_TMP(nBin).inLC_T = STATS.tstat; % ingroup lower vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outHL_H = H;
    Results_TMP(nBin).outHL_P = P;
    Results_TMP(nBin).outHL_T = STATS.tstat; % outgroup higher vs. lower

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,4), thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outHC_H = H;
    Results_TMP(nBin).outHC_P = P;
    Results_TMP(nBin).outHC_T = STATS.tstat; % outgroup higher vs. consistent 

    [H,P,CI,STATS] = ttest (thisBinNeuralRDM(:,5), thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outLC_H = H;
    Results_TMP(nBin).outLC_P = P;
    Results_TMP(nBin).outLC_T = STATS.tstat; % outgroup lower vs. consistent 

    Results_TMP(nBin).inH_M = mean (thisBinNeuralRDM(:,1));
    Results_TMP(nBin).inH_SE = std (thisBinNeuralRDM(:,1)) / length(thisBinNeuralRDM(:,1));
    Results_TMP(nBin).inL_M = mean (thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inL_SE = std (thisBinNeuralRDM(:,2)) / length(thisBinNeuralRDM(:,2));
    Results_TMP(nBin).inC_M = mean (thisBinNeuralRDM(:,3));
    Results_TMP(nBin).inC_SE = std (thisBinNeuralRDM(:,3)) / length(thisBinNeuralRDM(:,3));
    Results_TMP(nBin).outH_M = mean (thisBinNeuralRDM(:,4));
    Results_TMP(nBin).outH_SE = std (thisBinNeuralRDM(:,4)) / length(thisBinNeuralRDM(:,4));
    Results_TMP(nBin).outL_M = mean (thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outL_SE = std (thisBinNeuralRDM(:,5)) / length(thisBinNeuralRDM(:,5));
    Results_TMP(nBin).outC_M = mean (thisBinNeuralRDM(:,6));
    Results_TMP(nBin).outC_SE = std (thisBinNeuralRDM(:,6)) / length(thisBinNeuralRDM(:,6));
   
end
Results_Diff = Results_TMP;

cd (path_diff);
save(strcat(savName, '_Dissimilarity.mat'),...
    'Results_Pre', 'Results_Post', 'Results_Diff',...
    'allsubject_neuralrdm_allbin_Pre_AVG','allsubject_neuralrdm_allbin_Post_AVG','allsubject_neuralrdm_allbin_Diff_AVG');
writetable(struct2table(Results_Pre), strcat(savName,'_Dissimilarity_Pre.xlsx'));
writetable(struct2table(Results_Post), strcat(savName,'_Dissimilarity_Post.xlsx'));
writetable(struct2table(Results_Diff), strcat(savName,'_Dissimilarity_Diff.xlsx'));



