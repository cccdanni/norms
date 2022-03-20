% edited by Danni and edited from below source:
% for norms project.


% This is an updated script that computes average decoding accuracy 
% and performes cluster-based permutation analysis. 

% In order to perform cluster-based permutation analysis, users must have
% already simulated null-distribution to test against. See
% "SVM_ECOC_ERP_PermutationTesting_revised.". Plot function requires matlab
% function "boundedline" & "inpaint_nans" (Kearny 2020) in working directory. 
% see here: https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m

% Edited by Aaron Simmons (UC Davis)
% Original Author: Gi-Yeul Bae (Arizona State University)

% Note: Randomization routine included in this analysis (e.g., random integer generator)
% can produce stats slightly different from those reported in the paper.
% Gi-Yeul Bae  2017-10-3

%% Parameters

% Subject List 
subs = [1205:1242,1244:1253];
badsubs   = [1201:1204, 1214, 1215, 1238]; 
subList = setdiff (subs, badsubs);
Nsub = length(subList); % N subjects

% Decoding Parameters (must match svmECOC object)  
Nblock = 3; % cross-validation
Nitr = 10; % iteration
Ntp = 60; % # of time points; 
NBins = 2; % # of stimulus bins 

% Stats (see "Plot Significant Clusters" (line 204) for resulting plot formatting)
ClusterMassAnalysis = 0; % 0 = do not perform / 1 = perform
NPermutations = 100; % N permutations of null distribution
chancelvl = 1/NBins; %chance lvl of avg. decoding
tm = -200:20:996; % resampled timepoints from decoding (20ms/tp)
%define time point for statistical analysis
releventTime = 1:60; % corresponds to [220, 1480] ms
% find(tm==220); % obtain index of relevant time points in sampling units

workdir = "/home/chendanni/Documents/Norms/analysis/EEGAnalysisResults/Decoding";
% Individual subject filenames
condition = "ingroup-higher-vs-morph";
task = "preLearningImplicit";




task_list = {'preLearningImplicit', 'postLearningImplicit'};
condition_list = {'ingroup-higher-vs-morph', 'ingroup-lower-vs-morph','ingroup-consistent-vs-morph',...
    'outgroup-higher-vs-morph', 'outgroup-lower-vs-morph','outgroup-consistent-vs-morph'};


for i = 1:2
    task = task_list{i};
    for j = 1:6
        condition = condition_list{j};
        
        fileLocation = fullfile (workdir, 'Results', task, condition); 
        fName = strcat ('Results_ERPbased_', task , '_', condition, '_'); 


%% Plotting Decoding Results (no stats)

% Create/preallocate empty matrix
AverageAccuracy = nan(Nsub,Ntp); 

% First we will obtain average decoding accuracy within subjects at each
% timepoint. This results in size(AverageAccuracy). 

for sub = 1:Nsub 

    DecodingAccuracy = nan(Ntp,Nblock,Nitr);
    % We will compute decoding accuracy per subject in DecodingAccuracy,
    % enter DecodingAccuracy into AverageAccuray, then overwrite next subj.
    
    %% load SVM_ECOC output files
    readThis =strcat(fileLocation,filesep,fName,num2str(subList(sub)),'.mat');
    load(readThis)
     
    % Obtain predictions from SVM-ECOC model
    svmPrediction = squeeze(svmECOC.modelPredict);
    tstTargets = squeeze(svmECOC.targets);
    %clear svmECOC (is now a large unnecessary objet)
    clear svmECOC
    
    %% Step 5: Compute decoding accuracy of each decoding trial
    for block = 1:Nblock
        for itr = 1:Nitr
            for tp = 1:Ntp  

                prediction = squeeze(svmPrediction(itr,tp,block,:)); % this is predictions from models
                TrueAnswer = squeeze(tstTargets(itr,tp,block,:)); % this is predictions from models
                Err = TrueAnswer - prediction; %compute error. No error = 0
                ACC = mean(Err==0); %Correct hit = 0 (avg propotion of vector of 1s and 0s)
                DecodingAccuracy(tp,block,itr) = ACC; % average decoding accuracy at tp & block

            end
        end
    end
      
     % Average across block and iterations
     grandAvg = squeeze(mean(mean(DecodingAccuracy,2),3));
    
     % Perform temporal smoothing (5 point moving avg) 
     smoothed = nan(1,Ntp);
     for tAvg = 1:Ntp
         if tAvg ==1
           smoothed(tAvg) = mean(grandAvg((tAvg):(tAvg+2)));
         elseif tAvg ==2
           smoothed(tAvg) = mean(grandAvg((tAvg-1):(tAvg+2)));
         elseif tAvg == (Ntp-1)
           smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg+1)));
         elseif tAvg == Ntp
           smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg)));
         else
           smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg+2)));  
         end

     end
     
     % Save smoothe data
     AverageAccuracy(sub,:) =smoothed; % average across iteration and block
     
end %End of subject

%compute average accuracy across participants and SE of the mean.
cd (fileLocation);

subAverage = squeeze(mean(AverageAccuracy,1)); 
seAverage = squeeze(std(AverageAccuracy,1))/sqrt(Nsub); 
save (strcat(fName, 'all.mat'), 'subAverage', 'seAverage');

%% Visualization:Plotting the decoding accuracy 
% across all subjects at each timepoint
% Not publication quality
figure(1)
cl=colormap(parula(50));
plot(tm, subAverage); %plot
boundedline(tm,subAverage,seAverage,'cmap',cl(42,:),'alpha','transparency',0.35)
line([tm(1),tm(length(tm))],[chancelvl,chancelvl]); %chance line
saveas(figure(1),strcat(fName, 'all'),'pdf');
close all;

%% Perform cluster mass analyses
if ClusterMassAnalysis == 1
    
    % t-test at each relevent time point
    Ps = nan(2,length(releventTime));
    for i = 1:length(releventTime)
        tp = releventTime(i);
        
        [H,P,CI,STATS] =  ttest(AverageAccuracy(:,tp), chancelvl, 'tail', 'right'); % one sample t-test
        
        Ps(1,i) = STATS.tstat;
        Ps(2,i) = P;
    end
    
    % find significant time points
    candid = Ps(2,:) <= .05;
    
    %remove orphan time points
    candid_woOrphan = candid;
    candid_woOrphan(1,1) = candid(1,1);
    for i = 2:(length(releventTime)-1)
        
        if candid(1,i-1) == 0 && candid(1,i) ==1 && candid(1,i+1) ==0
            candid_woOrphan(1,i) = 0;
        else
            candid_woOrphan(1,i) = candid(1,i);
        end
        
    end
    
    % combine whole time range with relevent time & significant information
    clusters = zeros(length(tm),1);
    clusterT = zeros(length(tm),1);
    clusters(releventTime,1) = candid_woOrphan;
    clusterT(releventTime,1) = Ps(1,:);
    clusterTsum = sum(Ps(1,logical(candid_woOrphan)));
    
    %%find how many clusters are there, and compute summed T of each cluster
    tmp = zeros(10,300);
    cl = 0;
    member = 0;
    for i = 2:length(clusters)-1
        
        if clusters(i-1) ==0 && clusters(i) == 1 && clusters(i+1) == 1
            cl = cl+1;
            member = member +1;
            tmp(cl,member) = i;
            
        elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 0
            member = member +1;
            tmp(cl,member) = i;
            member = 0;
        elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 1
            member = member +1;
            tmp(cl,member) = i;
            
        else
            
        end
    end
    
    
    HowManyClusters = cl;
    a = tmp(1:cl,:); % subset significant clusters
    eachCluster = a(:,logical(sum(a,1)~=0)); % cut the size at the maximum cluster
    
    %now, compute summed T of each cluster
    dat_clusterSumT = nan(HowManyClusters,1);
    for c = 1:HowManyClusters
        dat_clusterSumT(c,1) = sum(clusterT(eachCluster(c,eachCluster(c,:) ~=0)));
    end
    
    
    
    %% note: Load Simulation tests (simulation takes very long time)    
    SimclusterTvalue = nan(1,NPermutations);

    for itr = 1:NPermutations
        
        Ps = nan(2,length(releventTime));

        % generate fake data. Size should be the same as actual data

          simacc = nan(Nitr,Nblock);
          simaccSub = nan(Nsub,Ntp);  
          for sub = 1:Nsub
              for t = 1:Ntp 
                for it = 1:Nitr
                    for blk = 1:Nblock
                    simaccTest = (1:NBins) - randi(NBins,1,NBins); % sample with replacement
                    simacc(it,blk) = sum(simaccTest==0)/NBins;
                    end
                 end
            simaccSub(sub,t) = mean(mean(simacc(:,:),1),2);
              end
          end

         % apply temporal smoothing to the fake data
          smtFake = nan(Nsub,Ntp);
          for tAvg = 1:Ntp
             if tAvg ==1
               smtFake(:,tAvg) = mean(simaccSub(:,(tAvg):(tAvg+2)),2);
             elseif tAvg ==2
               smtFake(:,tAvg) = mean(simaccSub(:,(tAvg-1):(tAvg+2)),2);
             elseif tAvg == (Ntp-1)
               smtFake(:,tAvg) = mean(simaccSub(:,(tAvg-2):(tAvg+1)),2);
             elseif tAvg == Ntp
               smtFake(:,tAvg) = mean(simaccSub(:,(tAvg-2):(tAvg)),2);
             else
               smtFake(:,tAvg) = mean(simaccSub(:,(tAvg-2):(tAvg+2)),2);
             end

          end
         % compute t value at relevent tiem points     
        for i = 1:length(releventTime) 
            tp = releventTime(i);

            [H,P,CI,STATS] =  ttest(smtFake(:,i),chancelvl, 'tail','right' );

            Ps(1,i) = STATS.tstat;
            Ps(2,i) = P;
        end

        % find significant points
        candid = Ps(2,:) <= .05;

        candid_woOrphan = zeros(1,length(candid));
        candid_woOrphan(1,1) = candid(1,1);
        candid_woOrphan(1,length(candid)) = candid(1,length(candid));
        %remove orphan time points
        for i = 2:(length(releventTime)-1)

            if candid(1,i-1) == 0 && candid(1,i) ==1 && candid(1,i+1) ==0
            candid_woOrphan(1,i) = 0; 
            else
            candid_woOrphan(1,i) = candid(1,i);     
            end

        end

        % combine whole time range with relevent time & significant information
        clusters = zeros(length(tm),1);
        clusterT = zeros(length(tm),1);
        clusters(releventTime,1) = candid_woOrphan;
        clusterT(releventTime,1) = Ps(1,:);


        %find how many clusters are there, and compute summed T of each cluster
        tmp = zeros(50,300);
        cl = 0;
        member = 0;
        for i = 2:length(clusters)-1

                if clusters(i-1) ==0 && clusters(i) == 1 && clusters(i+1) == 1 
                cl = cl+1;
                member = member +1;
                tmp(cl,member) = i;    

                elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 0 
                cl = cl+1;
                member = member +1;  
                tmp(cl,member) = i;    
                member = 0;  

                elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 1             
                member = member +1;  
                tmp(cl,member) = i;    

                else

                end
        end

        HowManyClusters = cl;
        if HowManyClusters >0
        a = tmp(1:cl,:);
        sim_eachCluster = a(:,logical(sum(a,1)~=0));

        %now, compute summed T of each cluster 
        sim_clusterSumT = zeros(HowManyClusters,1);
        for c = 1:HowManyClusters
           sim_clusterSumT(c,1) = sum(clusterT(sim_eachCluster(c,sim_eachCluster(c,:) ~=0)));
        end

        %find the most extreme cluster-level t value
        record = abs(sim_clusterSumT) == max(abs(sim_clusterSumT));
        SimclusterTvalue(1,itr) = sim_clusterSumT(record);
        else
        SimclusterTvalue(1,itr) = 0;    
        end

    end % End of simulation

    % sort simulated summed t values
    sortedTvalues = sort(SimclusterTvalue,2);

    %% find critical t-value
    permutedT = sortedTvalues;
    cutOff = NPermutations - NPermutations * 0.05; %one tailed
    critT = permutedT(cutOff); % t-mass of top 95%
    sigCluster = dat_clusterSumT > critT;
    
    
    %% Plot significant clusters
    figure(2)
    cl=colormap(parula(50));
    
    % draw average accuracy
    accEst = squeeze(subAverage);
    % draw clusters
    draw = eachCluster(sigCluster,:);
    draw = sort(reshape(draw,1,size(draw,1)*size(draw,2)));
    draw = draw(draw>0);
    
    w = zeros(Ntp,1);
    w(draw)=1;
    a = area(1:length(tm), accEst.*w');
    a.EdgeColor = 'none';
    a.FaceColor = [0.8,0.8,0.8];
    child = get(a,'Children');
    set(child,'FaceAlpha',0.9)
    % draw mean and SE of average decoding accuracy
    hold on
    mEI = boundedline(1:length(tm),subAverage,seAverage, 'cmap',cl(42,:),'alpha','transparency',0.35);
    
    %% Plot Formatting 
    xlabel('Time (ms)');ylabel('Decoding Accuracy')
    ax = gca;
    ax.YLim = [0.05, 0.15];
    ax.YTick = [0,0.02,0.04,0.06,0.08,0.10,0.12,0.14];
    ax.XTick = [1 26 51 76 100];
    ax.XTickLabel = {'-500','0','500','1000','1500'};
    h = line(1:length(tm),chancelvl* ones(1,Ntp));
    h.LineStyle = '--';
    h.Color = [0.1,0.1,0.1];
    hold off
    
    saveas(figure(2),'plot_decoding_accuracy_with_stat','pdf')
    
    close all;

end

    end
end
