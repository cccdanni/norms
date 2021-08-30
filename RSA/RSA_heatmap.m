%% RSA figure 
load ('post_alltrials_full_trans.mat');

% Behavioral RDM
allsubject_behardm_mean = mean(allsubject_behardm, 3);

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
imagesc (allsubject_behardm_mean)
colorbar
title ('Post-learning Attractiveness Rating')
set(gca,'xtick',[5:10:80],'xticklabel',CustomXLabels,...
    'ytick',[5:10:80],'yticklabel',CustomXLabels)
xtickangle(45)
hold on 
for i = [0:10:80]
    for j = [0:10:80]
        rectangle ('Position', [i+.5 j+.5  10 10], 'LineWidth', 2);
    end
end

% Neural RDM
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
title ('Post-learning Neural Dissimilarity')
set(gca,'xtick',[5:10:80],'xticklabel',CustomXLabels,...
    'ytick',[5:10:80],'yticklabel',CustomXLabels)
xtickangle(45)
hold on 
for i = [0:10:80]
    for j = [0:10:80]
        rectangle ('Position', [i+.5 j+.5  10 10], 'LineWidth', 2);
    end
end

%% RSA figure 
load ('pre_alltrials_full_trans.mat');

% Behavioral RDM
allsubject_behardm_mean = mean(allsubject_behardm, 3);

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
imagesc (allsubject_behardm_mean)
colorbar
title ('Pre-learning Attractiveness Rating')
set(gca,'xtick',[5:10:80],'xticklabel',CustomXLabels,...
    'ytick',[5:10:80],'yticklabel',CustomXLabels)
xtickangle(45)
hold on 
for i = [0:10:80]
    for j = [0:10:80]
        rectangle ('Position', [i+.5 j+.5  10 10], 'LineWidth', 2);
    end
end

% Neural RDM
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
title ('Pre-learning Neural Dissimilarity')
set(gca,'xtick',[5:10:80],'xticklabel',CustomXLabels,...
    'ytick',[5:10:80],'yticklabel',CustomXLabels)
xtickangle(45)
hold on 
for i = [0:10:80]
    for j = [0:10:80]
        rectangle ('Position', [i+.5 j+.5  10 10], 'LineWidth', 2);
    end
end