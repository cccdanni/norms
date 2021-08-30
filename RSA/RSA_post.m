project_folder = 'C:/Users/Psychology/Desktop/Research/social-norms-learning/Exp1_Norms_Learning/rawData'
rsa_folder = fullfile(project_folder, "RSA", "data");
erp_folder = fullfile(project_folder, "EEGAnalysisResults");
beha_folder = fullfile(project_folder, 'Norm_BehaData');

cd( beha_folder);
behadata = readmatrix ('PostLearningExplicitBehaData2021-02-23.csv');
behadata = behadata (:,[8,30,81]);

cd( rsa_folder );

%% post-learning 

this_rsa_folder = fullfile(rsa_folder, 'all_trials_post_sr100windowfullslidenone');
cd( this_rsa_folder );

subs = [1204:1242, 1244:1252];
badsubs = [1228,1229,1237,1239];
subs = setdiff(subs, badsubs);

allsubject_neuralrdm = zeros (80, 80, numel(subs));
allsubject_behardm = zeros (80, 80, numel(subs));
allsubject_coceprdm_conflict = zeros (80, 80, numel(subs));
allsubject_coceprdm_inconflict = zeros (80, 80, numel(subs));

allsubject = struct;

% load EEGlab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

cnt = 1;

for sub_nr = 1:numel(subs)
    thissub = num2str(subs(sub_nr));
    
    cd(this_rsa_folder);
    load ([thissub 'post_alltrials_full.mat']); % corr_res, erpdata, erpdata3
    
    corr_res = 1-corr_res; % simialrity to RDM
    
    thisERPfolder = fullfile (erp_folder, thissub);
    cd(thisERPfolder);
    EEG = pop_loadset ([thissub '_postLearningImplicit_EpochArtRej_FalseRespRej.set']);
    
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
    thisall  = [thisbini;thiscode]';
    
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
    end % conceptual rdm 1 
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
    end % conceptual rdm 2
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
    save([thissub 'post_alltrials_full_trans.mat'], 'neural_rdm', 'thisbinlist',...
        'conceptual_rdm_conflict', 'conceptual_rdm_inconflict', 'behavior_rdm');
    clear thissub corr_res corr_res_update thisbinlist newseq;
end

save('post_alltrials_full_trans.mat', 'allsubject_behardm', 'allsubject_neuralrdm',....
    'allsubject_coceprdm_inconflict', 'allsubject_coceprdm_conflict', 'allsubject');