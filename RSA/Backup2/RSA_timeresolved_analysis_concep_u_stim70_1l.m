% This is a script for conducting Statistical analysis with RSA results
% Author: Danni Chen, Ziqing Yao & Xiaoqing Hu (The University of Hong
% Kong) 
% Reference: Fellner, M. C., Waldhauser, G. T., & Axmacher, N.
% (2020). Tracking selective rehearsal and active inhibition of memory
% traces in directed forgetting. Current Biology, 30(13), 2638-2644. 
% Update date: 12/01/2021

%% define basic information 

clear;clc;

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
win = 0.20;
slide = 0.01;
nCond = 7;

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
rand_set = 2000;
p_cluster = 0.05;
p_first   = 0.05;

%% define path_in & path_out

path_in  = data_folder;
savName = strcat ('Results_RSA1l_concepImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_'); 
path_pre  = fullfile(saveLocation, 'Results',strcat('RSA1l_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
path_post = fullfile(saveLocation, 'Results',strcat('RSA1l_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
path_out = fullfile(saveLocation, 'Stats',  strcat('RSA1l_concepImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide),'_', ztrans));
if ~exist (path_out)
    mkdir(path_out)
end
cd (path_out);

%% load neural RDM of post-learning and pre-learning

RDM_pre = load ( strcat (path_pre, '/', 'Results_RSA1l_preLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_trans.mat') );
RDM_post = load ( strcat (path_post, '/', 'Results_RSA1l_postLearningImplicit_SR',num2str(sr),'WIN',num2str(win),'SLIDE',num2str(slide), '_', ztrans, '_trans.mat') );
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
allsubject_neuralbeha_R = nan (size (RDM_pre.allsubject_neuralrdm_allbins, 4),  n_bins);
allsubject_neuralconp_all_R = nan (size (RDM_pre.allsubject_neuralrdm_allbins, 4),  n_bins);
allsubject_neuralconp_in_R = nan (size (RDM_pre.allsubject_neuralrdm_allbins, 4),  n_bins);
allsubject_neuralconp_out_R = nan (size (RDM_pre.allsubject_neuralrdm_allbins, 4),  n_bins);
nBins = size (RDM_pre.allsubject_neuralrdm_allbins, 3);


%% Pre-learning ------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% construct the correlation matrix between conceptual RDM and neural RDM

fprintf ("-----------START OF PRELEARNING ANALYSIS------------\n")

fprintf ("construct the correlation matrix between conceptual RDM and neural RDM\n");
tStart = tic;

% Prepare dataset
behardm = RDM_pre.allsubject_behardm_allbins;
neurdm  = RDM_pre.allsubject_neuralrdm_allbins;
inrdm   = RDM_pre.allsubject_coceprdm_inconflict_allbins;
outrdm  = RDM_pre.allsubject_coceprdm_outconflict_allbins; 
allrdm  = RDM_pre.allsubject_coceprdm_conflict_allbins;
this_session = "prelearning";
this_task    = "implicit_perception";

% Correlation between neural rdm and conceptual rdm (in and out conflict)
con_vec = allrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_all_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

% Correlation between neural rdm and conceptual rdm (in conflict only)
con_vec = inrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_in_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;


% Correlation between neural rdm and conceptual rdm (out conflict only)
con_vec = outrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_out_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

% Correlation between neural rdm and behavioral rdm 
con_vec = behardm;
neu_vec = neurdm;
neu_vec(isnan(allrdm)) = nan;
con_vec(isnan(allrdm)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralbeha_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

toc(tStart);

%% Ststistical Testing
Results_Stats = struct; 
contrast_names = {'neural-conp-all', 'neural-conp-in', 'neural-conp-out', 'neural-beha', 'in-vs-out'};
all_R = {allsubject_neuralconp_all_R,...
    allsubject_neuralconp_in_R,... 
    allsubject_neuralconp_out_R,... 
    allsubject_neuralbeha_R,...
    (allsubject_neuralconp_in_R-allsubject_neuralconp_out_R)};
conp_RDM = {allrdm, inrdm, outrdm, behardm};  

for iContrast = 1 : length (contrast_names)
    
    % Conduct statistical test 
    [H, P, ~, STATS] = ttest (all_R{iContrast}, 0, 'alpha', p_first);% test whether the R is higher than zero
    Results_Stats(iContrast).H = H;
    Results_Stats(iContrast).P = P;
    Results_Stats(iContrast).STATS = STATS;
    Results_Stats(iContrast).conditions = contrast_names (iContrast);
    Results_Stats(iContrast).session    = this_session;
    Results_Stats(iContrast).task       = this_task;

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
        
    %% Permutation test

    fprintf (strcat("running permutation test on: ", contrast_names (iContrast),"\n"));
    tStart = tic;
    
    if isequal(data_post, 0) && isequal (data_negt, 0) % none of the cluster reach significance before correction
        rand = 1;
    else
        rand = rand_set;
    end
    
    if strcmp (perm_test, 'yes')
        
        pos_tsum = {};
        neg_tsum = {};
        
        gate_cmp = 0;
        if strcmp (contrast_names (iContrast), 'neural-beha') 
            gate_cmp = 1;
        end
            
        if ~strcmp (contrast_names (iContrast), 'in-vs-out')
            
            thisconprdm = conp_RDM{iContrast};
            
            for r = 1:rand
                
                r_perm = nan (nSubs, nBins); 
            
                %  shuffle (face group based shuffle)
                if gate_cmp 
                    randseq = randperm (nCond*10); % if neural-behavioral correlation: randomize the sequence of all stim
                else
                    randseq = randperm (nCond);
                    randseq_mat = repmat(randseq, 10,1);
                    a = (1:10)';
                    for j = 1:length(randseq)
                        randseq_mat(:,j) = (randseq_mat(:,j) - 1)*10 + a;
                    end
                    randseq = reshape (randseq_mat, nCond*10, 1);
                    randseq = randseq';
                end % if neural-conceptual correlation: shuffle group label
                
                % make conceptual rdm nan

                shuffledneurdm = neurdm(randseq, :, :, :); % only shuffle the x-axis, not y-axis
                shuffledneurdm( isnan (thisconprdm)) = nan;

                % reshape neural rdm & conceptual rdm
                shuffledneurdm_reshape = reshape (shuffledneurdm, size (shuffledneurdm, 1)*size (shuffledneurdm, 2), size (shuffledneurdm, 3), size (shuffledneurdm, 4) );
                shuffledconrdm_reshape = reshape (thisconprdm, size(thisconprdm,1)*size(thisconprdm,2), size(thisconprdm,3), size(thisconprdm,4));

                % correlation 
                for iwBin = 1:nBins
                    for iwSub = 1:nSubs
                        r_perm (iwSub, iwBin) = corr(shuffledneurdm_reshape(:, iwBin, iwSub), shuffledconrdm_reshape(:, iwBin, iwSub), 'rows', 'complete', 'type','Spearman');
                    end
                end % now we have a nSubs * nBins correlation coefficient matrix
                
                % conduct statistics analysis
                [H, P, ~, STATS] = ttest (r_perm, 0, 'alpha', p_first);

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
                
                clear shuffledneurdm_reshape shuffledconrdm_reshape r_perm randseq;
                clear STATS rand_neg rand_L_neg rand_num_neg rand_negt rand_pos rand_L_pos rand_num_pos m rand_post;  
                
            end
            
        elseif strcmp (contrast_names (iContrast), 'in-vs-out')
            
            dat = [];
            dat1 = [];
            dat2 = [];
            
            dat1 = squeeze (all_R{2});
            dat2 = squeeze (all_R{3});
            dat(:,:,1) = dat1;
            dat(:,:,2) = dat2;
            
            for r = 1:rand
                            
                % randomize conditions 
                randsel = randi([0,1], nSubs, nBins);  % if x (in randsel) = 1 then dat1-dat2; = 2 then dat2-dat1.
                randdat1 = (dat1.*randsel) + (dat2.*(~randsel)); 
                randdat2 = (dat1.*(~randsel)) + (dat2.*(randsel)); % here we randomize label (cond1 vs cond2)

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
                
                clear shuffledneurdm_reshape shuffledconrdm_reshape r_perm randseq randsel randsel1 randsel2s randAVG;
                clear STATS rand_neg rand_L_neg rand_num_neg rand_negt rand_pos rand_L_pos rand_num_pos m rand_post;  
                
            end

        end
        
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

        toc (tStart); 

        % significant?
        % check positive clusters
        ind=1;
        p_threshold=p_cluster;
        sig_pos=zeros(size(data_pos));
        pos_cluster(ind)=nearest(pos_dist(ind,:),data_post(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(pos_dist(ind,:),data_post(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            p_threshold=p_threshold/2;
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
        p_threshold=p_cluster;
        neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(neg_dist(ind,:),data_negt(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            p_threshold=p_threshold/2;
            sig_neg=sig_neg+(data_L_neg==ind_negt(ind)); % locate cluster
            neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand; % calculate p value
            ind=ind+1;
        if ind > numel(data_negt) | ind > max_neg
            break % if all cluster for both the data or the shuffled data has been tested, then break 
        end
        end 

        Results_Stats(iContrast).cluster_pos.perm_dist   = pos_dist;
        Results_Stats(iContrast).cluster_pos.perm_p      = pos_cluster;
        Results_Stats(iContrast).cluster_pos.sig_p       = sig_pos;
        Results_Stats(iContrast).cluster_neg.perm_dist   = neg_dist;
        Results_Stats(iContrast).cluster_neg.perm_p      = neg_cluster;
        Results_Stats(iContrast).cluster_neg.sig_p       = sig_neg;

        clear ind p_threshold

        alpha=sig_pos+sig_neg; % combining p_higher and p_lower 
        alpha_raw = Results_Stats(iContrast).H; % combining p_higher and p_lower (not corrected)

        cd (path_out);
        
        cmp = customcolormap_preset('red-white-blue');
    
        fig =  figure
        
        subplot(2,2,1:2)
        tm  = 1: n_bins; 
        M1  = mean (all_R{iContrast}); 
        SE1 = std(all_R{iContrast})/sqrt(nSubs);
        hAVG = alpha;
        boundedline(tm,M1,SE1,'cmap',cmp,'alpha','transparency',0.35);
        ylim ([-1.0000e-02,1.0000e-02])
        title(strcat(Results_Stats(iContrast).session, ' ', Results_Stats(iContrast).conditions));
        ylabel ('correlation');
        xlabel ('timepoint');
        hold on;
        line ([1,nBins], [0,0], 'linestyle', '--', 'LineWidth', 1, 'color','b')
        liney = 0;
        [L, num] = bwlabel (hAVG);
        for n = 1:num
            line (find(L == n), liney*ones(size(find(L == n))), 'LineWidth', 3, 'color','r')
        end
        hold on;
        [L, num] = bwlabel (alpha_raw);
        for n = 1:num
            line (find(L == n), (-0.005)*ones(size(find(L == n))), 'LineWidth', 3, 'color','k')
        end

        subplot(2,2,3)
        hist(pos_dist(1,:))
        hold on
        plot([data_post(1) data_post(1)],[0 rand],'r')
        title('distribution positive cluster')
        text(data_post(1),10,strcat('pcorr=',num2str(pos_cluster)))

        subplot(2,2,4)
        hist(neg_dist(1,:))
        hold on
        plot([data_negt(1) data_negt(1)],[0 rand],'r')
        title('distribution negative cluster')
        text(data_negt(1), 10,strcat('pcorr=',num2str(neg_cluster)))
        
        
        savefig(fig,fullfile(path_out,strcat('clustered_correlation_', Results_Stats(iContrast).session, '_', Results_Stats(iContrast).conditions)))
        
        close all;
        
        clear dat1 dat2 M1 M2 MAll SE1 SE2 hAVG alpha sig_pos sig_neg;

    end

end

Results_Pre_Stats = Results_Stats;

cd (path_out);
save(strcat(savName, strcat(Results_Pre_Stats(iContrast).session, '.mat')),...
    'Results_Pre_Stats');
clear behardm neurdm inrdm outrdm allrdm conprdm_all conprdm_in conprdm_out neuralrdm Results_Stats;



fprintf ("-----------END OF PRELEARNING ANALYSIS------------\n")


fprintf ("-----------START OF POSTLEARNING ANALYSIS------------\n")


%% Post-learning ------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% construct the correlation matrix between conceptual RDM and neural RDM

fprintf ("construct the correlation matrix between conceptual RDM and neural RDM\n");
tStart = tic;

% Prepare dataset
clear behardm neurdm inrdm outrdm allrdm this_session this_task
behardm = RDM_post.allsubject_behardm_allbins;
neurdm  = RDM_post.allsubject_neuralrdm_allbins;
inrdm   = RDM_post.allsubject_coceprdm_inconflict_allbins;
outrdm  = RDM_post.allsubject_coceprdm_outconflict_allbins; 
allrdm  = RDM_post.allsubject_coceprdm_conflict_allbins;
this_session = "postlearning";
this_task    = "implicit_perception";

% Correlation between neural rdm and conceptual rdm (in and out conflict)
con_vec = allrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_all_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

% Correlation between neural rdm and conceptual rdm (in conflict only)
con_vec = inrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_in_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;


% Correlation between neural rdm and conceptual rdm (out conflict only)
con_vec = outrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_out_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

% Correlation between neural rdm and behavioral rdm 
con_vec = behardm;
neu_vec = neurdm;
neu_vec(isnan(allrdm)) = nan;
con_vec(isnan(allrdm)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralbeha_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

toc(tStart);

%% Ststistical Testing
Results_Stats = struct; 
contrast_names = {'neural-conp-all', 'neural-conp-in', 'neural-conp-out', 'neural-beha', 'in-vs-out'};
all_R = {allsubject_neuralconp_all_R,...
    allsubject_neuralconp_in_R,... 
    allsubject_neuralconp_out_R,... 
    allsubject_neuralbeha_R,...
    (allsubject_neuralconp_in_R-allsubject_neuralconp_out_R)};
conp_RDM = {allrdm, inrdm, outrdm, behardm};  

for iContrast = 1 : length (contrast_names)
    
    % Conduct statistical test 
    [H, P, ~, STATS] = ttest (all_R{iContrast}, 0, 'alpha', p_first);% test whether the R is higher than zero
    Results_Stats(iContrast).H = H;
    Results_Stats(iContrast).P = P;
    Results_Stats(iContrast).STATS = STATS;
    Results_Stats(iContrast).conditions = contrast_names (iContrast);
    Results_Stats(iContrast).session    = this_session;
    Results_Stats(iContrast).task       = this_task;

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
        
    %% Permutation test

    fprintf (strcat("running permutation test on: ", contrast_names (iContrast),"\n"));
    tStart = tic;
    
    if isequal(data_post, 0) && isequal (data_negt, 0) % none of the cluster reach significance before correction
        rand = 1;
    else
        rand = rand_set;
    end
    
    if strcmp (perm_test, 'yes')
        
        pos_tsum = {};
        neg_tsum = {};
        
        gate_cmp = 0;
        if strcmp (contrast_names (iContrast), 'neural-beha') 
            gate_cmp = 1;
        end
            
        if ~strcmp (contrast_names (iContrast), 'in-vs-out')
            
            thisconprdm = conp_RDM{iContrast};
            
            for r = 1:rand
                
                r_perm = nan (nSubs, nBins); 
            
                %  shuffle (face group based shuffle)
                if gate_cmp 
                    randseq = randperm (nCond*10); % if neural-behavioral correlation: randomize the sequence of all stim
                else
                    randseq = randperm (nCond);
                    randseq_mat = repmat(randseq, 10,1);
                    a = (1:10)';
                    for j = 1:length(randseq)
                        randseq_mat(:,j) = (randseq_mat(:,j) - 1)*10 + a;
                    end
                    randseq = reshape (randseq_mat, nCond*10, 1);
                    randseq = randseq';
                end % if neural-conceptual correlation: shuffle group label
                
                % make conceptual rdm nan

                shuffledneurdm = neurdm(randseq, :, :, :); % only shuffle the x-axis, not y-axis
                shuffledneurdm( isnan (thisconprdm)) = nan;

                % reshape neural rdm & conceptual rdm
                shuffledneurdm_reshape = reshape (shuffledneurdm, size (shuffledneurdm, 1)*size (shuffledneurdm, 2), size (shuffledneurdm, 3), size (shuffledneurdm, 4) );
                shuffledconrdm_reshape = reshape (thisconprdm, size(thisconprdm,1)*size(thisconprdm,2), size(thisconprdm,3), size(thisconprdm,4));

                % correlation 
                for iwBin = 1:nBins
                    for iwSub = 1:nSubs
                        r_perm (iwSub, iwBin) = corr(shuffledneurdm_reshape(:, iwBin, iwSub), shuffledconrdm_reshape(:, iwBin, iwSub), 'rows', 'complete', 'type','Spearman');
                    end
                end % now we have a nSubs * nBins correlation coefficient matrix
                
                % conduct statistics analysis
                [H, P, ~, STATS] = ttest (r_perm, 0, 'alpha', p_first);

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
                
                clear shuffledneurdm_reshape shuffledconrdm_reshape r_perm randseq;
                clear STATS rand_neg rand_L_neg rand_num_neg rand_negt rand_pos rand_L_pos rand_num_pos m rand_post;  
                
            end
            
        elseif strcmp (contrast_names (iContrast), 'in-vs-out')
            
            dat = [];
            dat1 = [];
            dat2 = [];
            
            dat1 = squeeze (all_R{2});
            dat2 = squeeze (all_R{3});
            dat(:,:,1) = dat1;
            dat(:,:,2) = dat2;
            
            for r = 1:rand
                            
                % randomize conditions 
                randsel = randi([0,1], nSubs, nBins);  % if x (in randsel) = 1 then dat1-dat2; = 2 then dat2-dat1.
                randdat1 = (dat1.*randsel) + (dat2.*(~randsel)); 
                randdat2 = (dat1.*(~randsel)) + (dat2.*(randsel)); % here we randomize label (cond1 vs cond2)

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
                
                clear shuffledneurdm_reshape shuffledconrdm_reshape r_perm randseq randsel randsel1 randsel2s randAVG;
                clear STATS rand_neg rand_L_neg rand_num_neg rand_negt rand_pos rand_L_pos rand_num_pos m rand_post;  
                
            end

        end
        
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

        toc (tStart); 

        % significant?
        % check positive clusters
        ind=1;
        p_threshold=p_cluster;
        sig_pos=zeros(size(data_pos));
        pos_cluster(ind)=nearest(pos_dist(ind,:),data_post(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(pos_dist(ind,:),data_post(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            p_threshold=p_threshold/2;
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
        p_threshold=p_cluster;
        neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(neg_dist(ind,:),data_negt(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            p_threshold=p_threshold/2;
            sig_neg=sig_neg+(data_L_neg==ind_negt(ind)); % locate cluster
            neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand; % calculate p value
            ind=ind+1;
        if ind > numel(data_negt) | ind > max_neg
            break % if all cluster for both the data or the shuffled data has been tested, then break 
        end
        end 

        Results_Stats(iContrast).cluster_pos.perm_dist   = pos_dist;
        Results_Stats(iContrast).cluster_pos.perm_p      = pos_cluster;
        Results_Stats(iContrast).cluster_pos.sig_p       = sig_pos;
        Results_Stats(iContrast).cluster_neg.perm_dist   = neg_dist;
        Results_Stats(iContrast).cluster_neg.perm_p      = neg_cluster;
        Results_Stats(iContrast).cluster_neg.sig_p       = sig_neg;

        clear ind p_threshold

        alpha=sig_pos+sig_neg; % combining p_higher and p_lower 
        alpha_raw = Results_Stats(iContrast).H; % combining p_higher and p_lower (not corrected)

        cd (path_out);
        
        cmp = customcolormap_preset('red-white-blue');
    
        fig =  figure
        
        subplot(2,2,1:2)
        tm  = 1: n_bins; 
        M1  = mean (all_R{iContrast}); 
        SE1 = std(all_R{iContrast})/sqrt(nSubs);
        hAVG = alpha;
        boundedline(tm,M1,SE1,'cmap',cmp,'alpha','transparency',0.35);
        ylim ([-1.0000e-02,1.0000e-02])
        title(strcat(Results_Stats(iContrast).session, ' ', Results_Stats(iContrast).conditions));
        ylabel ('correlation');
        xlabel ('timepoint');
        hold on;
        line ([1,nBins], [0,0], 'linestyle', '--', 'LineWidth', 1, 'color','b')
        liney = 0;
        [L, num] = bwlabel (hAVG);
        for n = 1:num
            line (find(L == n), liney*ones(size(find(L == n))), 'LineWidth', 3, 'color','r')
        end
        hold on;
        [L, num] = bwlabel (alpha_raw);
        for n = 1:num
            line (find(L == n), (-0.005)*ones(size(find(L == n))), 'LineWidth', 3, 'color','k')
        end

        subplot(2,2,3)
        hist(pos_dist(1,:))
        hold on
        plot([data_post(1) data_post(1)],[0 rand],'r')
        title('distribution positive cluster')
        text(data_post(1),10,strcat('pcorr=',num2str(pos_cluster)))

        subplot(2,2,4)
        hist(neg_dist(1,:))
        hold on
        plot([data_negt(1) data_negt(1)],[0 rand],'r')
        title('distribution negative cluster')
        text(data_negt(1), 10,strcat('pcorr=',num2str(neg_cluster)))
        
        
        savefig(fig,fullfile(path_out,strcat('clustered_correlation_', Results_Stats(iContrast).session, '_', Results_Stats(iContrast).conditions)))
        
        close all;
        
        clear dat1 dat2 M1 M2 MAll SE1 SE2 hAVG alpha sig_pos sig_neg;

    end

end

Results_Post_Stats = Results_Stats;

cd (path_out);
save(strcat(savName, strcat(Results_Stats(iContrast).session, '.mat')), 'Results_Post_Stats');
clear behardm neurdm inrdm outrdm allrdm conprdm_all conprdm_in conprdm_out neuralrdm Results_Stats;

fprintf ("-----------END OF POSTLEARNING ANALYSIS------------\n")





fprintf ("-----------START OF DIFFERENCE ANALYSIS------------\n")


%% Difference: Post-Pre ------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% construct the correlation matrix between conceptual RDM and neural RDM

fprintf ("construct the correlation matrix between conceptual RDM and neural RDM\n");
tStart = tic;

% Prepare dataset
clear behardm neurdm inrdm outrdm allrdm this_session this_task
behardm = RDM_post.allsubject_behardm_allbins - RDM_pre.allsubject_behardm_allbins;
neurdm  = RDM_post.allsubject_neuralrdm_allbins - RDM_pre.allsubject_neuralrdm_allbins;
inrdm   = RDM_post.allsubject_coceprdm_inconflict_allbins;
outrdm  = RDM_post.allsubject_coceprdm_outconflict_allbins; 
allrdm  = RDM_post.allsubject_coceprdm_conflict_allbins;
this_session = "difference";
this_task    = "implicit_perception";

% Correlation between neural rdm and conceptual rdm (in and out conflict)
con_vec = allrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_all_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

% Correlation between neural rdm and conceptual rdm (in conflict only)
con_vec = inrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_in_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;


% Correlation between neural rdm and conceptual rdm (out conflict only)
con_vec = outrdm;
neu_vec = neurdm;
neu_vec(isnan(con_vec)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralconp_out_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

% Correlation between neural rdm and behavioral rdm 
con_vec = behardm;
neu_vec = neurdm;
neu_vec(isnan(allrdm)) = nan;
con_vec(isnan(allrdm)) = nan;
neu_vec2 = reshape (neu_vec, [ size(neu_vec,1) * size(neu_vec,2), size(neu_vec,3)*size(neu_vec,4) ]); 
con_vec2 = reshape (con_vec, [ size(con_vec,1) * size(con_vec,2), size(con_vec,3)*size(con_vec,4) ]); 
corr_tmp = corr (neu_vec2, con_vec2, 'rows', 'complete', 'type', 'Spearman');
corr_tmp = reshape (corr_tmp, [ size(neu_vec,3), size(neu_vec,4), size(neu_vec,3), size(neu_vec,4) ] );
corr_tmp2 = nan (size(neu_vec,3), size(neu_vec,4));
for is = 1:nSubs
    for ib = 1:nBins
        corr_tmp2 (ib, is) = corr_tmp(ib, is, ib, is);
    end
end
allsubject_neuralbeha_R = corr_tmp2';
clear con_vec neu_vec con_vec2 neu_vec2 corr_tmp corr_tmp2;

toc(tStart);

%% Ststistical Testing
Results_Stats = struct; 
contrast_names = {'neural-conp-all', 'neural-conp-in', 'neural-conp-out', 'neural-beha', 'in-vs-out'};
all_R = {allsubject_neuralconp_all_R,...
    allsubject_neuralconp_in_R,... 
    allsubject_neuralconp_out_R,... 
    allsubject_neuralbeha_R,...
    (allsubject_neuralconp_in_R-allsubject_neuralconp_out_R)};
conp_RDM = {allrdm, inrdm, outrdm, behardm};  

for iContrast = 1 : length (contrast_names)
    
    % Conduct statistical test 
    [H, P, ~, STATS] = ttest (all_R{iContrast}, 0, 'alpha', p_first);% test whether the R is higher than zero
    Results_Stats(iContrast).H = H;
    Results_Stats(iContrast).P = P;
    Results_Stats(iContrast).STATS = STATS;
    Results_Stats(iContrast).conditions = contrast_names (iContrast);
    Results_Stats(iContrast).session    = this_session;
    Results_Stats(iContrast).task       = this_task;

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
        
    %% Permutation test

    fprintf (strcat("running permutation test on: ", contrast_names (iContrast),"\n"));
    tStart = tic;
    
    if isequal(data_post, 0) && isequal (data_negt, 0) % none of the cluster reach significance before correction
        rand = 1;
    else
        rand = rand_set;
    end
    
    if strcmp (perm_test, 'yes')
        
        pos_tsum = {};
        neg_tsum = {};
        
        gate_cmp = 0;
        if strcmp (contrast_names (iContrast), 'neural-beha') 
            gate_cmp = 1;
        end
            
        if ~strcmp (contrast_names (iContrast), 'in-vs-out')
            
            thisconprdm = conp_RDM{iContrast};
            
            for r = 1:rand
                
                r_perm = nan (nSubs, nBins); 
            
                %  shuffle (face group based shuffle)
                if gate_cmp 
                    randseq = randperm (nCond*10); % if neural-behavioral correlation: randomize the sequence of all stim
                else
                    randseq = randperm (nCond);
                    randseq_mat = repmat(randseq, 10,1);
                    a = (1:10)';
                    for j = 1:length(randseq)
                        randseq_mat(:,j) = (randseq_mat(:,j) - 1)*10 + a;
                    end
                    randseq = reshape (randseq_mat, nCond*10, 1);
                    randseq = randseq';
                end % if neural-conceptual correlation: shuffle group label
                
                % make conceptual rdm nan

                shuffledneurdm = neurdm(randseq, :, :, :); % only shuffle the x-axis, not y-axis
                shuffledneurdm( isnan (thisconprdm)) = nan;

                % reshape neural rdm & conceptual rdm
                shuffledneurdm_reshape = reshape (shuffledneurdm, size (shuffledneurdm, 1)*size (shuffledneurdm, 2), size (shuffledneurdm, 3), size (shuffledneurdm, 4) );
                shuffledconrdm_reshape = reshape (thisconprdm, size(thisconprdm,1)*size(thisconprdm,2), size(thisconprdm,3), size(thisconprdm,4));

                % correlation 
                for iwBin = 1:nBins
                    for iwSub = 1:nSubs
                        r_perm (iwSub, iwBin) = corr(shuffledneurdm_reshape(:, iwBin, iwSub), shuffledconrdm_reshape(:, iwBin, iwSub), 'rows', 'complete', 'type','Spearman');
                    end
                end % now we have a nSubs * nBins correlation coefficient matrix
                
                % conduct statistics analysis
                [H, P, ~, STATS] = ttest (r_perm, 0, 'alpha', p_first);

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
                
                clear shuffledneurdm_reshape shuffledconrdm_reshape r_perm randseq;
                clear STATS rand_neg rand_L_neg rand_num_neg rand_negt rand_pos rand_L_pos rand_num_pos m rand_post;  
                
            end
            
        elseif strcmp (contrast_names (iContrast), 'in-vs-out')
            
            dat = [];
            dat1 = [];
            dat2 = [];
            
            dat1 = squeeze (all_R{2});
            dat2 = squeeze (all_R{3});
            dat(:,:,1) = dat1;
            dat(:,:,2) = dat2;
            
            for r = 1:rand
                            
                % randomize conditions 
                randsel = randi([0,1], nSubs, nBins);  % if x (in randsel) = 1 then dat1-dat2; = 2 then dat2-dat1.
                randdat1 = (dat1.*randsel) + (dat2.*(~randsel)); 
                randdat2 = (dat1.*(~randsel)) + (dat2.*(randsel)); % here we randomize label (cond1 vs cond2)

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
                
                clear shuffledneurdm_reshape shuffledconrdm_reshape r_perm randseq randsel randsel1 randsel2s randAVG;
                clear STATS rand_neg rand_L_neg rand_num_neg rand_negt rand_pos rand_L_pos rand_num_pos m rand_post;  
                
            end

        end
        
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

        toc (tStart); 

        % significant?
        % check positive clusters
        ind=1;
        p_threshold=p_cluster;
        sig_pos=zeros(size(data_pos));
        pos_cluster(ind)=nearest(pos_dist(ind,:),data_post(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(pos_dist(ind,:),data_post(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            p_threshold=p_threshold/2;
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
        p_threshold=p_cluster;
        neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand;  % nearest: which randomized cluster T is closest to the current data
        while nearest(neg_dist(ind,:),data_negt(1))<=round(p_threshold.*rand) % when the rank of the closest number is smaller than critical rank (i.e., if p_threshold = 0.05, rand = 100, for the cluster is significant, the cluster T should be higher than the top 5th randomized T, hence the rank of nearest number should be smaller than 5)
            p_threshold=p_threshold/2;
            sig_neg=sig_neg+(data_L_neg==ind_negt(ind)); % locate cluster
            neg_cluster(ind)=nearest(neg_dist(ind,:),data_negt(ind))./rand; % calculate p value
            ind=ind+1;
        if ind > numel(data_negt) | ind > max_neg
            break % if all cluster for both the data or the shuffled data has been tested, then break 
        end
        end 

        Results_Stats(iContrast).cluster_pos.perm_dist   = pos_dist;
        Results_Stats(iContrast).cluster_pos.perm_p      = pos_cluster;
        Results_Stats(iContrast).cluster_pos.sig_p       = sig_pos;
        Results_Stats(iContrast).cluster_neg.perm_dist   = neg_dist;
        Results_Stats(iContrast).cluster_neg.perm_p      = neg_cluster;
        Results_Stats(iContrast).cluster_neg.sig_p       = sig_neg;

        clear ind p_threshold

        alpha=sig_pos+sig_neg; % combining p_higher and p_lower 
        alpha_raw = Results_Stats(iContrast).H; % combining p_higher and p_lower (not corrected)

        cd (path_out);
        
        cmp = customcolormap_preset('red-white-blue');
    
        fig =  figure
        
        subplot(2,2,1:2)
        tm  = 1: n_bins; 
        M1  = mean (all_R{iContrast}); 
        SE1 = std(all_R{iContrast})/sqrt(nSubs);
        hAVG = alpha;
        boundedline(tm,M1,SE1,'cmap',cmp,'alpha','transparency',0.35);
        ylim ([-1.0000e-02,1.0000e-02])
        title(strcat(Results_Stats(iContrast).session, ' ', Results_Stats(iContrast).conditions));
        ylabel ('correlation');
        xlabel ('timepoint');
        hold on;
        line ([1,nBins], [0,0], 'linestyle', '--', 'LineWidth', 1, 'color','b')
        liney = 0;
        [L, num] = bwlabel (hAVG);
        for n = 1:num
            line (find(L == n), liney*ones(size(find(L == n))), 'LineWidth', 3, 'color','r')
        end
        hold on;
        [L, num] = bwlabel (alpha_raw);
        for n = 1:num
            line (find(L == n), (-0.005)*ones(size(find(L == n))), 'LineWidth', 3, 'color','k')
        end

        subplot(2,2,3)
        hist(pos_dist(1,:))
        hold on
        plot([data_post(1) data_post(1)],[0 rand],'r')
        title('distribution positive cluster')
        text(data_post(1),10,strcat('pcorr=',num2str(pos_cluster)))

        subplot(2,2,4)
        hist(neg_dist(1,:))
        hold on
        plot([data_negt(1) data_negt(1)],[0 rand],'r')
        title('distribution negative cluster')
        text(data_negt(1), 10,strcat('pcorr=',num2str(neg_cluster)))
        
        
        savefig(fig,fullfile(path_out,strcat('clustered_correlation_', Results_Stats(iContrast).session, '_', Results_Stats(iContrast).conditions)))
        
        close all;
        
        clear dat1 dat2 M1 M2 MAll SE1 SE2 hAVG alpha sig_pos sig_neg;

    end

end

Results_Diff_Stats = Results_Stats;

cd (path_out);
save(strcat(savName, strcat(Results_Stats(iContrast).session, '.mat')), 'Results_Diff_Stats');
clear behardm neurdm inrdm outrdm allrdm conprdm_all conprdm_in conprdm_out neuralrdm Results_Stats;

fprintf ("-----------END OF DIFFERENCE ANALYSIS------------\n")

