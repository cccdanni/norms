%% first step:
project_folder='\data_share\';
toolbox_folder='\matlab_tools\';
%% add toolboxes
addpath (fullfile(toolbox_folder,'fieldtrip-20190611'))
% add path with additional functions
addpath (fullfile(project_folder,'scripts','additional_functions'));
%% draw trials for alpha, levene test
% test for stable variance in subsamples
% draw n times a certain number of trials
path_in=fullfile(project_folder,'data');

all_subs={'01';'02';'03';'04';'05';'06';'08';'09';'12';'13';'14';'15';'16';'17';'18';'19';'22';'23'};

nrand=20;
num2draw=5:25;
toi=[2 4];

for n=1:numel(all_subs)
    sel_sub=all_subs{n};
    load(fullfile(path_in,strcat(sel_sub,'_data')));
    
    % get alpha power for these trials
    cfg=[];
    cfg.method='wavelet';
    cfg.foi=8:0.5:13;
    cfg.toi=-0.5:0.05:4;
    cfg.width=5;
    cfg.output='pow';
    cfg.keeptrials='yes';
    cfg.channel='POz',
    data=ft_freqanalysis(cfg, data);
    cfg=[];
    cfg.time=[-0.5 4.5];
    data=z_trans_TF_seltime(cfg,data);
    t(1)=nearest(data.time,toi(1));
    t(2)=nearest(data.time,toi(2));
    
    freq=squeeze(nanmean(nanmean(nanmean(data.powspctrm(:,:,:,t(1):t(2)),2),3),4));
    true_var(n)=var(freq);
    true_mean(n)=mean(freq);
    % trialnumbers
    for draw=1:numel(num2draw)
        for r=1:nrand
            sel_draw=num2draw(draw);
            randvec=randperm(numel(sel_trials),sel_draw);
            avg_pow(n,draw,r)=nanmean(freq(randvec));
            
        end
    end
end


%calculate variance
var_all = var(avg_pow,0,3);

% compute levene test for every subject/ number vs 30
for n=1:numel(all_subs)
    for draw=1:(numel(num2draw)-1)
        pow_vec=[squeeze(avg_pow(n,draw,:));squeeze(avg_pow(n,end,:))];
        group_vec=[ones(nrand,1).*draw;ones(nrand,1).*num2draw(end)];
        p(n,draw) = vartestn(pow_vec,group_vec,'TestType','LeveneAbsolute','Display', 'off');
    end
end

for draw=1:(numel(num2draw)-1)
    group_pval(draw) = fisher_pvalue_meta_analysis( p(:,draw));
end

min_trials =[10 18 14 17 23 30 10 22 15 15 27 19 16 35 13 11 15 22];

figure
subplot(1,4,1)
mf_lineplusdots(ones(1,numel(num2draw)).*true_mean(5),squeeze(avg_pow(5,:,:))',num2draw,'mean power')
title('means for every draw for one subject')

subplot(1,4,2)
imagesc(num2draw(1:end-1),1:numel(vp),var_all,[0 0.1])
colormap(gca,'parula')
ylabel('single subjects')
xlabel('number of trials')
title(['variance across random draws'])

subplot(1,4,3)
imagesc(num2draw(1:end-1),1:numel(vp),p,[0 1])
colormap(gca,'bone')
ylabel('single subjects')
xlabel('number of trials')
title(['p-value of levene test:','variance in x number trials vs 20 trials'])
hold on
scatter(min_trials,1:numel(vp),'o','r','filled')

subplot(1,4,4)
yyaxis left
plot(num2draw,mean(var_all))
ylabel('mean variance')
hold on
yyaxis right
plot(num2draw(1:end-1),sum(p>=0.05)/numel(vp))

ylabel('relative number of subjects with non sig levene test')
ylim([0 1])
