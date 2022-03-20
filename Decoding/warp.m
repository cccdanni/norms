work_dir = '/home/chendanni/Documents/Norms/analysis/EEGAnalysisResults';

subs = [1205:1242,1244:1253];
badsubs   = [1201:1204, 1214, 1215, 1238]; % 1201 - 1204: HK Local; 1214, 1215 & 1238 - visually detected bad subjects
subs = setdiff (subs, badsubs);

task_list = {'preLearningImplicit', 'postLearningImplicit'};
condition_list = {'ingroup-higher-vs-morph', 'ingroup-lower-vs-morph','ingroup-consistent-vs-morph',...
    'outgroup-higher-vs-morph', 'outgroup-lower-vs-morph','outgroup-consistent-vs-morph'};

for i = 1:2
    task = task_list{i};
    for j = 1:6
        condition = condition_list{j};
        SVM_ECOC_ERP_Decoding_morphvsexp(subs, work_dir, condition, task);
    end
end
        

