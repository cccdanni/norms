% This is a script for conducting Statistical analysis with RSA results
% Author: Danni Chen, Ziqing Yao & Xiaoqing Hu (The University of Hong
% Kong) 
% Reference: Fellner, M. C., Waldhauser, G. T., & Axmacher, N.
% (2020). Tracking selective rehearsal and active inhibition of memory
% traces in directed forgetting. Current Biology, 30(13), 2638-2644. 
% Update date: 12/01/2021

%% define basic information 
subs = [1205:1242,1244:1253];
badsubs   = [1201:1204, 1220, 1229, 1238, 1239, 1249]; % 1201 - 1204: Local, 1229, 1238, 1239,1249 (reject epochs > 100)
subs = setdiff (subs, badsubs);
nSubs = length(subs); %number of subjects

work_dir = '/home/chendanni/Documents/Norms/analysis/';
cd (work_dir);

erp  = "full_range";

stat_test = "yes"; 
perm_test = "yes";
ztrans = "yes"; 
average_erp = "no";

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
win = 0.15;
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

%% define path_in & path_out

path_in  = data_folder;
savName = strcat ('Results_RSA_concepImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 

path_pre  = fullfile(saveLocation, 'Results',strcat('RSA_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
path_post = fullfile(saveLocation, 'Results',strcat('RSA_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
path_diff = fullfile(saveLocation, 'Stats',strcat('RSA_concepImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
if ~exist (path_diff)
    mkdir(path_diff)
end
path_out = path_diff;

%% load neural RDM of post-learning and pre-learning

RDM_pre = load ( strcat (path_pre, '/', 'Results_RSA_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_trans.mat') );
RDM_post = load ( strcat (path_post, '/', 'Results_RSA_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_trans.mat') );
% allsubject_behardm_allbins: nItems x nItems x nBins x nSubject
% allsubject_conceprdm_conflict_allbins: nItems x nItems x nBins x nSubject
% allsubject_conceprdm_inconflict_allbins: nItems x nItems x nBins x nSubject
% allsubject_conceprdm_outconflict_allbins: nItems x nItems x nBins x nSubject
% allsubject_neuralrdm_allbin: nItems x nItems x nBins x nSubject

concep_rdm_both = squeeze (RDM_pre.allsubject_coceprdm_conflict_allbins (:,:,1,1));
concep_rdm_both_vector = matrix2vector (concep_rdm_both, "lower");
concep_rdm_in   = squeeze (RDM_pre.allsubject_coceprdm_inconflict_allbins (:,:,1,1));
concep_rdm_in_vector   = matrix2vector (concep_rdm_in, "lower");
concep_rdm_out  = squeeze (RDM_pre.allsubject_coceprdm_outconflict_allbins (:,:,1,1));
concep_rdm_out_vector  = matrix2vector (concep_rdm_out, "lower");

Results_NeuralConpAll = struct;
Results_NeuralConpIn  = struct;
Results_NeuralConpOut = struct;
Results_NeuralConpDiff= struct;

%% correlate neural RDM and conceptual RDM
allsubject_neuralbeha_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
allsubject_neuralconp_all_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);
allsubject_neuralconp_in_R = nan (size (allsubject_neuralrdm_allbins, 4),  n_bins);

for nSub = 1: size (allsubject_neuralrdm_allbins, 4)

    for nBin = 1: size (allsubject_neuralrdm_allbins, 3)

        behardm     = allsubject_behardm_allbins (:, :, nBin, nSub);
        conprdm_all = allsubject_coceprdm_conflict_allbins (:, :, nBin, nSub);
        conprdm_in  = allsubject_coceprdm_inconflict_allbins (:, :, nBin, nSub);
        conprdm_out = allsubject_coceprdm_outconflict_allbins (:, :, nBin, nSub);
        neuralrdm   = allsubject_neuralrdm_allbins (:, :, nBin, nSub);

        %% Correlation between neural rdm and conceptual rdm (in and out conflict)
        conprdm_VECTOR = matrix2vector (conprdm_all, "lower");

        %% Correlation between neural rdm and conceptual rdm (in conflict only)
        conprdm_VECTOR = matrix2vector (conprdm_in, "lower");

        %% Correlation between neural rdm and conceptual rdm (out conflict only)
        conprdm_VECTOR = matrix2vector (conprdm_out, "lower");




    end

end


[neuralbeha_corr, neuralbeha_corr_p]  = corr(neuralrdm_VECTOR', behardm', 'type','Spearman');


%% significant Analysis 














for nBin = 1:n_bins

    %% Correlation between neural rdm and conceptual rdm (in and out conflict)
    [Results_NeuralConpAll(nBin).H, Results_NeuralConpAll(nBin).P] = ttest (allsubject_neuralconp_all_R(:,nBin), 0, 'tail', 'right');
    Results_NeuralConpAll(nBin).R = mean (allsubject_neuralconp_all_R(:,nBin));
    Results_NeuralConpAll(nBin).SE = std (allsubject_neuralconp_all_R(:,nBin)) / sqrt(length (allsubject_neuralconp_all_R(:,nBin)));

    %% Correlation between neural rdm and conceptual rdm (in conflict only)
    [Results_NeuralConpIn(nBin).H, Results_NeuralConpIn(nBin).P] = ttest (allsubject_neuralconp_in_R(:,nBin), 0, 'tail', 'right');
    Results_NeuralConpIn(nBin).R = mean (allsubject_neuralconp_in_R(:,nBin));
    Results_NeuralConpIn(nBin).SE = std (allsubject_neuralconp_in_R(:,nBin)) / sqrt(length (allsubject_neuralconp_in_R(:,nBin)));

    %% Correlation between neural rdm and conceptual rdm (out conflict only)
    [Results_NeuralConpOut(nBin).H, Results_NeuralConpOut(nBin).P] = ttest (allsubject_neuralconp_out_R(:,nBin), 0, 'tail', 'right');
    Results_NeuralConpOut(nBin).R = mean (allsubject_neuralconp_out_R(:,nBin));
    Results_NeuralConpOut(nBin).SE = std (allsubject_neuralconp_out_R(:,nBin)) / sqrt(length (allsubject_neuralconp_out_R(:,nBin)));

    %% Compare the model between in vs. out 
    [Results_NeuralConpDiff(nBin).H, Results_NeuralConpDiff(nBin).P] = ttest (allsubject_neuralconp_in_R(:,nBin), allsubject_neuralconp_out_R(:,nBin), 'tail', 'both');
    Results_NeuralConpDiff(nBin).R_In   = mean (allsubject_neuralconp_in_R(:,nBin));
    Results_NeuralConpDiff(nBin).SE_In  = std (allsubject_neuralconp_in_R(:,nBin)) / sqrt(length (allsubject_neuralconp_in_R(:,nBin)));
    Results_NeuralConpDiff(nBin).R_Out  = mean (allsubject_neuralconp_out_R(:,nBin));
    Results_NeuralConpDiff(nBin).SE_Out = std (allsubject_neuralconp_out_R(:,nBin)) / sqrt(length (allsubject_neuralconp_out_R(:,nBin)));
end

%% Significant Analysis 

contrast_pairs = [1 2; 1 3; 3 2; 4 5; 4 6; 6 1];
contrast_names = {'in_HL', 'in_HC', 'in_CL', 'out_HL', 'out_HC', 'out_CL'};
contrast_cond_names = {'Higher','Lower'; 'Higher','Consistent'; 'Consistent','Lower';...
    'Higher','Lower'; 'Higher','Consistent'; 'Consistent','Lower'};

%% pre-learning
% Statistics Testing
Results_Stats = struct; 
for iContrast = 1:6
    
    DiffAVG = squeeze (allsubject_neuralrdm_allbin_Pre_AVG(:,contrast_pairs(iContrast, 1),:) - allsubject_neuralrdm_allbin_Pre_AVG(:,contrast_pairs(iContrast, 2),:));
    [H, P, ~, STATS] = ttest (DiffAVG, 0, 'alpha', p_first);
    Results_Stats(iContrast).H = H;
    Results_Stats(iContrast).P = P;
    Results_Stats(iContrast).STATS = STATS;
    Results_Stats(iContrast).conditions = contrast_names (iContrast);
    Results_Stats(iContrast).session    = "prelearning";
    
    % Find negative cluster 
    data_neg=squeeze(H.*(STATS.tstat<0));
    data_neg(isnan(data_neg)) = 0; % find significant and t<0 time-window 
    [data_L_neg,data_num_neg] = bwlabel(data_neg); % bwlabel: Label connected components in 2-D binary image.
    for neg=1:data_num_neg
        m=find(data_L_neg==neg);
        data_negt(neg)=sum(STATS.tstat(m));
    end
    if isempty(neg)
        data_negt=0;
    end
    [data_negt,ind_negt]=sort(data_negt,'ascend'); % abs(t) biggest
    
    % Find positive cluster 
    data_pos=squeeze(H.*(STATS.tstat>0));
    data_pos(isnan(data_pos)) = 0; % find significant and t>0 time-window 
    [data_L_pos,data_num_pos] = bwlabel(data_pos); % find cluster, data_L_pos: label components in the cluster; data_num_pos: how many clusters in the data
    for pos=1:data_num_pos
        m=find(data_L_pos==pos);
        data_post(pos)=sum(STATS.tstat(m));
    end % calculate the cluster_T (sum of T within each cluster)
    if isempty(pos)
        data_post=0;
    end
    [data_post,ind_post]=sort(data_post,'descend');
    
    Results_Stats(iContrast).cluster_neg.L       = data_L_neg;
    Results_Stats(iContrast).cluster_neg.num     = data_num_neg;
    Results_Stats(iContrast).cluster_neg.T_value = data_negt;
    Results_Stats(iContrast).cluster_neg.T_ind   = ind_negt;
    Results_Stats(iContrast).cluster_pos.L       = data_L_pos;
    Results_Stats(iContrast).cluster_pos.num     = data_num_pos;
    Results_Stats(iContrast).cluster_pos.T_value = data_post;
    Results_Stats(iContrast).cluster_pos.T_ind   = ind_post; % save all data in the Results_Stats
        
    % Permutation test
    if strcmp (perm_test, 'yes')
        
        dat1 = squeeze (allsubject_neuralrdm_allbin_Pre_AVG(:,contrast_pairs(iContrast, 1),:));
        dat2 = squeeze (allsubject_neuralrdm_allbin_Pre_AVG(:,contrast_pairs(iContrast, 2),:));
        dat(:,:,1) = dat1;
        dat(:,:,2) = dat2;
        
        pos_tsum = {}; 
        neg_tsum = {};
        
        for r = 1:rand
            
            % randomize conditions 
            randsel = randi([1,2], nSubs, n_bins);  % if x (in randsel) = 1 then dat1-dat2; = 2 then dat2-dat1.
            randsel1 = (randsel == 1);
            randsel2 = (randsel == 2);
            randdat1 = (dat1.*randsel1) + (dat2.*randsel2); 
            randdat2 = (dat1.*(~randsel1)) + (dat2.*(~randsel2)); % here we randomize label (cond1 vs cond2)
            
            % conduct statistics analysis
            randAVG = randdat1 - randdat2;
            [H, P, ~, STATS] = ttest (randAVG, 0, 'alpha', p_first);

            % Find negative cluster 
            rand_neg=squeeze(H.*(STATS.tstat<0));
            rand_neg(isnan(rand_neg)) = 0; 
            [rand_L_neg,rand_num_neg] = bwlabel(rand_neg); 
            for neg=1:rand_num_neg
                m=find(rand_L_neg==neg);
                rand_negt(neg)=sum(STATS.tstat(m));
            end
            if isempty(neg)
                rand_negt=0;
            end
    
            % Find positive cluster 
            rand_pos=squeeze(H.*(STATS.tstat>0));
            rand_pos(isnan(rand_pos)) = 0; 
            [rand_L_pos,rand_num_pos] = bwlabel(rand_pos); 
            for pos=1:rand_num_pos
                m=find(rand_L_pos==pos);
                rand_post(pos)=sum(STATS.tstat(m));
            end 
            if isempty(pos)
                rand_post=0;
            end
            
            pos_tsum {r} = sort (rand_post, 'descend');
            neg_tsum {r} = sort (rand_negt, 'ascend'); 
            
            clear randsel randsel1 randsel2 randdat1 randdat2 randAVG H P;
            clear STATS rand_neg rand_L_neg rand_num_neg rand_negt rand_pos rand_L_pos rand_num_pos m rand_post; 
            
        end % randomize RAND times, and found tsum in each data 
        
        % randomized cluster 
        max_pos=max(cellfun(@numel,pos_tsum)); % max positive cluster number 
        max_neg=max(cellfun(@numel,neg_tsum)); % max negative cluster number 
        
        for x=1:rand
            pos_tsum{x}= [pos_tsum{x},zeros(1,max_pos-numel(pos_tsum{x}))];
            neg_tsum{x}= [neg_tsum{x},zeros(1,max_neg-numel(neg_tsum{x}))];
        end % complete the matrix for the next step 
        pos_dist=reshape([pos_tsum{:}],max_pos,rand); % maximum cluster number * rand
        pos_dist=sort(pos_dist,2,'descend');
        neg_dist=reshape([ neg_tsum{:}],max_neg,rand); % maximum cluster number * rand
        neg_dist=sort(neg_dist,2,'ascend');
        
        % significant?
        % check positive clusters
        ind=1;
        p_threshold=p_cluster;
        p_threshold=p_threshold/2; % one sided test; DC: I edited here, so that we could compare the p_threshold
        sig_pos=zeros(size(data_pos));
        pos_cluster(ind)=nearest(pos_dist(ind,:),data_post(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(pos_dist(ind,:),data_post(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            sig_pos=sig_pos+(data_L_pos==ind_post(ind)); % locate cluster
            pos_cluster(ind)=nearest(pos_dist(ind,:),data_post(ind))./rand; % calculate p value
            ind=ind+1;
        if ind > numel(data_post) | ind > max_pos
            break % if all cluster for both the data or the shuffled data has been tested, then break 
        end
        end 
        
        % check negative clusters
        ind=1;
        sig_neg=zeros(size(data_neg));
        neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(neg_dist(ind,:),data_negt(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            sig_neg=sig_neg+(data_L_neg==ind_negt(ind)); % locate cluster
            neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand; % calculate p value
            ind=ind+1;
        if ind > numel(data_negt) | ind > max_neg
            break % if all cluster for both the data or the shuffled data has been tested, then break 
        end
        end 

        
        Results_Stats(iContrast).cluster_pos.perm_dist   = pos_dist;
        Results_Stats(iContrast).cluster_pos.perm_p      = pos_cluster;
        Results_Stats(iContrast).cluster_neg.perm_dist   = neg_dist;
        Results_Stats(iContrast).cluster_neg.perm_p      = neg_cluster;
        
        clear ind p_threshold
        
        alpha=sig_pos+sig_neg; % combining p_higher and p_lower 
        
        cd (path_out)
        tm = 1: n_bins; 
        M1 = mean (dat1); 
        M2 = mean (dat2);
        MAll = [M1, M2];
        SE1 = std(dat1)/sqrt(nSubs);
        SE2 = std(dat2)/sqrt(nSubs);
        hAVG = alpha;
        figure(2)
        cl=colormap(parula(50));
        boundedline(tm,M1,SE1,'cmap',cl(42,:),'alpha','transparency',0.35);
        hold on;
        boundedline(tm,M2,SE2,'cmap',cl(21,:),'alpha','transparency',0.35);
        legend (contrast_cond_names{iContrast,1},'',contrast_cond_names{iContrast,2},'');
        ylabel ('dissimilarity');
        xlabel ('timepoint');
        liney = min (MAll, [], 'all');
        [L, num] = bwlabel (hAVG);
        for n = 1:num
            line (find(L == n), liney*ones(size(find(L == n))), 'LineWidth', 5, 'color','r')
        end
        title(strcat(Results_Stats(iContrast).session, string(Results_Stats(iContrast).conditions)));
        saveas(gcf,strcat(savName,strcat(Results_Stats(iContrast).session, string(Results_Stats(iContrast).conditions), '.png')));
        close all;
        
        clear dat1 dat2 M1 M2 MAll SE1 SE2 hAVG alpha sig_pos sig_neg liney;
        
    end 
    
end
Results_Pre_Stats = Results_Stats; 

cd (path_diff);
save(strcat(savName, 'Dissimilarity.mat'),...
    'Results_Pre_Stats', 'Results_Post_Stats', 'Results_Diff_Stats', ...
    'allsubject_neuralrdm_allbin_Pre_AVG','allsubject_neuralrdm_allbin_Post_AVG','allsubject_neuralrdm_allbin_Diff_AVG');
% writetable(struct2table(Results_Pre), strcat(savName,'Dissimilarity_Pre.xlsx'));
% writetable(struct2table(Results_Post), strcat(savName,'Dissimilarity_Post.xlsx'));
% writetable(struct2table(Results_Diff), strcat(savName,'_Dissimilarity_Diff.xlsx'));




