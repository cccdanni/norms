%% define work folder 

project_folder = "/home/chendanni/Documents/Norms/data";
rsa_folder = fullfile(project_folder, "EEGAnalysisResults_PreprocessingVersion1", "RSA", "data");
erp_folder = fullfile(project_folder, "EEGAnalysisResults_PreprocessingVersion1");


cd ( rsa_folder );

pre_rsa_folder = fullfile(rsa_folder, "all_trials_pre_sr100window200slide10");
post_rsa_folder = fullfile(rsa_folder, "all_trials_post_sr100window200slide10");

%% define parameters

subs = [1204:1242, 1244:1253];
badsubs = [1228, 1237, 1239, 1229];
subs = setdiff(subs, badsubs);
all_subs = {};
for i = 1:length(subs)
    all_subs (i) =  { num2str(subs(i)) } ;
end

stim_nr = 80;

% define sliding windo
win = 0.02; % in sec 
slide = 0.01;

% define new Sampling rate
sr=100; % in Hz
step=win/(1/sr)+1;

channels='all';

%define start and end of item window
t_start= 0;
t_end=0.996;
tois=t_start:1/sr:t_end;
t1=t_start:slide:(t_end-win);
t2=t1+win;
ind_t1=1:slide/(1/sr):((numel(tois)-win/(1/sr)));
ind_t2=ind_t1+win/(1/sr);
n_bins=numel(t1);


%% prepare behavioral rdm and conpcetual rdm for each ppts
% 
% pre_alltrials_full_trans  = readmatrix ("/home/chendanni/Documents/Norms/data/BehaData/PreLearningExplicitBehaData2021-02-23.csv");
% post_alltrials_full_trans = readmatrix ("/home/chendanni/Documents/Norms/data/BehaData/PostLearningExplicitBehaData2021-02-23.csv");
% 
% 
% for sub_nr = 1:numel(subs)
%     
%     thissub = num2str(subs(sub_nr));
%     load (['all_trials_post_sr100window200slide10/' thissub 'post_alltrials_timeresolved.mat']);
%     
%     thisERPfolder = fullfile (erp_folder, thissub);
%     cd(thisERPfolder);
%     EEG = pop_loadset ([thissub '_postLearningImplicit_EpochArtRej_FalseRespRej.set']);
%     
%     thisbini = [];
%     thiscode = [];
%     for i = 1: size(EEG.event, 2)
%         thisbini(i) = EEG.event(i).bini;
%         thislabel = EEG.event(i).codelabel;
%         thiscode(i) = str2double(thislabel(2:end));
%     end
%     uniquecode = unique(thiscode, 'first');
%     thisbini = thisbini (uniquecode);
%     thiscode = thiscode (uniquecode);
%     thisall  = [thisbini;thiscode]';
%     % these code select unique bin 
%     
%     thisbinlist = zeros(8, 10);
%     for i = 1:8 % 8 bins
%         x = ceil((min(thiscode(find(thisbini == i))))/10);
%         x1 = ((x-1)*10+1):x*10;
%         thisbinlist(i,:) = x1;
%     end % check each bin includes which figures
%     newseq = reshape(thisbinlist',1,80);
%     
%     conceptual_rdm_allconf = ones (size(corr_res,1), size(corr_res,3));
%     for i = 1: size(corr_res,1)
%         for j = 1:size(corr_res,3)
%             if j>=i 
%                 conceptual_rdm_allconf (i, j ) = nan;
%             elseif (ismember(j, [1:10, 31:40])) & (ismember(i, [71:80]))
%                     conceptual_rdm_allconf (i, j ) = 0;
%             end
%         end
%     end % conceptual rdm 1: both ingroup higher and outgroup higher set as zeros (dissimilarity = 0)
%     allsubject_coceprdm_conf (:,:, sub_nr) = conceptual_rdm_allconf;
%     
%     conceptual_rdm_in_conf = ones (size(corr_res,1), size(corr_res,3));
%     for i = 1: size(corr_res,1)
%         for j = 1:size(corr_res,3)
%             if j>=i 
%                 conceptual_rdm_in_conf (i, j ) = nan;
%             elseif (ismember(j, [1:10])) & (ismember(i, [71:80]))
%                     conceptual_rdm_in_conf (i, j ) = 0;
%             end
%         end
%     end % conceptual rdm 2: only ingroup higher set as zeros (dissimilarity = 0)
%     allsubject_coceprdm_in_conf (:,:, sub_nr) = conceptual_rdm_in_conf;
%     
%     conceptual_rdm_outconf = ones (size(corr_res,1), size(corr_res,3));
%     for i = 1: size(corr_res,1)
%         for j = 1:size(corr_res,3)
%             if j>=i 
%                 conceptual_rdm_outconf (i, j ) = nan;
%             elseif (ismember(j, [31:40])) & (ismember(i, [71:80]))
%                     conceptual_rdm_outconf (i, j ) = 0;
%             end
%         end
%     end % conceptual rdm 3: only outgroup higher set as zeros (dissimilarity = 0)
%     allsubject_coceprdm_conf (:,:, sub_nr) = conceptual_rdm_outconf;
%     
%     behadata = pre_alltrials_full_trans(:,[8,30,81]);
%     behavior_rdm = zeros (size(corr_res,1), size(corr_res,3));
%     thisbeha = behadata(find(behadata(:,3) == subs(sub_nr)),:);
%     for i = 1:size(corr_res,1)
%         for j = 1:size(corr_res,3)
%             rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
%             rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
%             behavior_rdm(i,j) = abs(rating_i - rating_j);
%         end
%     end % behavioral rdm
%     pre_learning_beha_rdm = behavior_rdm;
%     allsubject_pre_behardm (:,:,sub_nr) = behavior_rdm; % behavioral rdm (prelearning)
%     
%     behadata = post_alltrials_full_trans(:,[8,30,81]);
%     behavior_rdm = zeros (size(corr_res,1), size(corr_res,3));
%     thisbeha = behadata(find(behadata(:,3) == subs(sub_nr)),:);
%     for i = 1:size(corr_res,1)
%         for j = 1:size(corr_res,3)
%             rating_i = thisbeha(find(thisbeha(:,1)==newseq(i)), 2);
%             rating_j = thisbeha(find(thisbeha(:,1)==newseq(j)), 2);
%             behavior_rdm(i,j) = abs(rating_i - rating_j);
%         end
%     end % behavioral rdm
%     post_learning_beha_rdm = behavior_rdm;
%     allsubject_post_behardm (:,:,sub_nr) = behavior_rdm; % behavioral rdm (postlearning)
%     
% 
%     cd( rsa_folder );
%     save([thissub 'concept_beha_rdm_trans.mat'], 'thisbinlist',...
%         'conceptual_rdm_allconf', 'conceptual_rdm_in_conf', 'pre_learning_beha_rdm','post_learning_beha_rdm','conceptual_rdm_outconf');
%     clear thissub corr_res corr_res_update thisbinlist newseq thisbinlist conceptual_rdm_allconf conceptual_rdm_in_conf pre_learning_beha_rdm post_learning_beha_rdm;
% end

%% --- Statistical analysis 1 
% statistical analysis 1: correlate neural RDM with two different conceptual RDM

cd ( rsa_folder );

conds = 2;
models = 4;
rdm_correlation = nan (numel(subs), n_bins, conds, models);

for sub_nr = 1:numel(subs)
    
    thissub = num2str(subs(sub_nr));
    
    thisrsa_pre  = load (strcat("all_trials_pre_sr100window200slide10/", thissub, "pre_alltrials_timeresolved.mat"));
    thisrsa_post = load (strcat("all_trials_post_sr100window200slide10/", thissub, "post_alltrials_timeresolved.mat"));
    concept_rsa  = load (strcat(thissub, "concept_beha_rdm_trans.mat"));
    
    thisrsa_pre  = thisrsa_pre.corr_res_trial;
    thisrsa_post = thisrsa_post.corr_res_trial;
    
    thisneural_rdm_pre  = 1-thisrsa_pre;
    thisneural_rdm_post = 1-thisrsa_post;
    thiscon_rdm_allconf = concept_rsa.conceptual_rdm_allconf;
    thiscon_rdm_inconf  = concept_rsa.conceptual_rdm_in_conf;
    thiscon_rdm_outconf  = concept_rsa.conceptual_rdm_outconf;
    thisbeha_rdm_pre    = concept_rsa.pre_learning_beha_rdm;
    thisbeha_rdm_post   = concept_rsa.post_learning_beha_rdm;
    
    for bin_nr = 1:n_bins
   
        for cond_nr = 1:conds
            
            switch cond_nr
                case 1
                    thiscond = "pre";
                    corr_neural_rdm = squeeze( thisneural_rdm_pre(bin_nr,:,:) );
                    corr_beha_rdm = thisbeha_rdm_pre;
                case 2
                    thiscond = "post";
                    corr_neural_rdm = squeeze( thisneural_rdm_post(bin_nr,:,:) );
                    corr_beha_rdm = thisbeha_rdm_post;
            end
            
            for model_nr = 1:models
                
                switch model_nr 
                    case 1
                        thismodel = "conceptual_rdm_all_conflict";
                        corr_concept_rdm = thiscon_rdm_allconf;
                        rdm1 = matrix2vector(corr_neural_rdm,'lower')';
                        rdm2 = matrix2vector(corr_concept_rdm,'lower')';
                        rdm_correlation (sub_nr, bin_nr, cond_nr, model_nr) = corr(rdm1, rdm2, 'type', 'Spearman');
                    case 2
                        thismodel = "conceptual_rdm_in_conflict";
                        corr_concept_rdm = thiscon_rdm_inconf;
                        rdm1 = matrix2vector(corr_neural_rdm,'lower')';
                        rdm2 = matrix2vector(corr_concept_rdm,'lower')';
                        rdm_correlation (sub_nr, bin_nr, cond_nr, model_nr) = corr(rdm1, rdm2, 'type', 'Spearman');
                    case 3 
                        thismodel = "conceptual_rdm_out_conflict";
                        corr_concept_rdm = thiscon_rdm_outconf;
                        rdm1 = matrix2vector(corr_neural_rdm,'lower')';
                        rdm2 = matrix2vector(corr_concept_rdm,'lower')';
                        rdm_correlation (sub_nr, bin_nr, cond_nr, model_nr) = corr(rdm1, rdm2, 'type', 'Spearman');
                    case 4
                        thismodel = "behavioral_model";
                        rdm1 = matrix2vector(corr_neural_rdm,'lower')';
                        rdm2 = matrix2vector(corr_beha_rdm,'lower')';
                        rdm_correlation (sub_nr, bin_nr, cond_nr, model_nr) = corr(rdm1, rdm2, 'type', 'Spearman');
                end
                
            end
            
        end
        
    end
    
end

save(strcat("RDM_Correlation",date(),'.mat'),'rdm_correlation');

for model_nr = 1:models
    
    switch model_nr 
        case 1
            thismodel = "ConceptualAllConflict";
        case 2
            thismodel = "ConceptualIngroupConflict";
        case 3
            thismodel = "ConceptualOutgroupConflict";
        case 4
            thismodel = "Behavioral"
    end    
    
    thiscorr = squeeze(mean(rdm_correlation(:,:,:,model_nr),1));
    
    figure;
    plot (thiscorr);
    lgd = legend ('pre', 'post');
    title (thismodel);
    xticks (1:5:size(thiscorr,1));
    ylabel ('r');
    xticklabels (string(t1(1:5:numel(t1))));
    saveas(gcf, strcat(thismodel,'_Correlation', date(),'.png'));
    close all;
    
end


%% --- Statistical analysis 2 
% statistical analysis 2: directly compare dissimilarity

cd ( rsa_folder );

conds = 6;
parts = 2;
neural_rdm_all = nan (numel(subs), n_bins, conds, parts);

for sub_nr = 1:numel(subs)
    
    thissub = num2str(subs(sub_nr));
    
    thisrsa_pre  = load (strcat("all_trials_pre_sr100window200slide10/", thissub, "pre_alltrials_timeresolved.mat"));
    thisrsa_post = load (strcat("all_trials_post_sr100window200slide10/", thissub, "post_alltrials_timeresolved.mat"));
    thisrsa_pre  = thisrsa_pre.corr_res_trial;
    thisrsa_post = thisrsa_post.corr_res_trial;
    thisneural_rdm_pre  = 1-thisrsa_pre;
    thisneural_rdm_post = 1-thisrsa_post;

    
    for bin_nr = 1:n_bins
   
        for cond_nr = 1:conds
            
            for part_nr = 1:parts
                
                switch part_nr 
                    case 1 
                        rdm = squeeze ( thisneural_rdm_pre(bin_nr,:,:));
                        neural_rdm_all(sub_nr, bin_nr, cond_nr, part_nr) = mean( mean(rdm(  (10*(cond_nr-1)+1) :(10*cond_nr), 71:80 )));
                        
                    case 2
                        rdm = squeeze ( thisneural_rdm_post(bin_nr, :, :));
                        neural_rdm_all(sub_nr, bin_nr, cond_nr, part_nr) = mean( mean(rdm(  (10*(cond_nr-1)+1) :(10*cond_nr), 71:80 )));
                end

                
            end
            
        end
        
    end
    
end

save(strcat("RDM_",date(),'.mat'),'neural_rdm_all');

%% Plot 
data_pre = neural_rdm_all(:,:,1:3,1);
data_post = neural_rdm_all(:,:,1:3,2);
data = data_post - data_pre;
data = squeeze (mean(data,1));
size(data);

figure;
plot (data);
lgd = legend ('higher', 'lower', 'consistent');
xticks (1:5:size(data,1));
ylabel ('Dissimilarity Update');
xticklabels (string(t1(1:5:numel(t1))));
saveas(gcf, strcat('Dissimilarity_Ingroup', date(),'.png'));
close all;



data_pre = neural_rdm_all(:,:,4:6,1);
data_post = neural_rdm_all(:,:,4:6,2);
data = data_post - data_pre;
data = squeeze (mean(data,1));
size(data);
figure;
plot (data);
lgd = legend ('higher', 'lower', 'consistent');
xticks (1:5:size(data,1));
ylabel ('Dissimilarity Update');
xticklabels (string(t1(1:5:numel(t1))));
saveas(gcf, strcat('Dissimilarity_Outgroup', date(),'.png'));
close all;