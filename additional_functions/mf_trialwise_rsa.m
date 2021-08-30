% rsa function
%data1 & data2: fieldtrip trial-based data (dimord: 'rpt_chan_time'

% cfg.trial_combi= 'in_trial' only on diagonal, important: order of trials in data1 and data2 need to match
%                 'across_trial' correlation between all trial combinations
% cfg_rsa.win=  window size of sliding window in sec
% cfg_rsa.slide= slide i.e. steps in which window is moved, in sec
% cfg_rsa.sr= sampling rate in Hz
% cfg_rsa.step=cfg_rsa.win/(1/cfg_rsa.sr)+1;
% cfg.corr_type='Spearman' 'Pearson'
%cfg definiton

function corr_trials=mf_trialwise_rsa(cfg,data1,data2)

%enco: define start and end of item window
t_start1=data1.time(1);
t_end1=data1.time(end-1);

tois1=t_start1:1/cfg.sr:t_end1;
t1_1=t_start1:cfg.slide:(t_end1-cfg.win);
t2_1=t1_1+cfg.win;

ind_t1_1=1:cfg.slide/(1/cfg.sr):((numel(tois1)-cfg.win/(1/cfg.sr)));
ind_t2_1=ind_t1_1+cfg.win/(1/cfg.sr);

n_bins1=numel(t1_1);


%cue define start and end of item window
t_start2=data2.time(1);
t_end2=data2.time(end-1);

tois2=t_start2:1/cfg.sr:t_end2;
t1_2=t_start2:cfg.slide:(t_end2-cfg.win);
t2_2=t1_2+cfg.win;

ind_t1_2=1:cfg.slide/(1/cfg.sr):((numel(tois2)-cfg.win/(1/cfg.sr)));
ind_t2_2=ind_t1_2+cfg.win/(1/cfg.sr);

n_bins2=numel(t1_2);


      data1_vec=zeros(size(data1.trial,1),n_bins1,numel(data1.label)*cfg.step);
    for bin=1:n_bins1
       % vectorize sel_bins: data_vec(trials, nbins, features)
       data_vec_tmp=data1.trial(:,:,ind_t1_1(bin):ind_t2_1(bin));
       data1_vec(:,bin,:)=reshape(data_vec_tmp,size(data1.trial,1),[]);
       % repmat  
     %  data_vec_tmp=dataenc.trial(:,:,ind_t1_e(bin):ind_t2_e(bin));
     %  data_enc_vec_nonorm(:,bin,:)=reshape(data_vec_tmp,num_1+num_2,[]);
    end
    clear data_vec_tmp 
   

      data2_vec=zeros(size(data2.trialinfo,1),n_bins2,numel(data2.label)*cfg.step);
    for bin=1:n_bins2
       % vectorize sel_bins: data_vec(trials, nbins, features)
       data_vec_tmp=data2.trial(:,:,ind_t1_2(bin):ind_t2_2(bin));
       data2_vec(:,bin,:)=reshape(data_vec_tmp,size(data2.trial,1),[]);
       % repmat
     %  data_vec_tmp=datacue.trial(:,:,ind_t1_c(bin):ind_t2_c(bin));
      % data_cue_vec_nonorm(:,bin,:)=reshape(data_vec_tmp,num_1+num_2,[]);       
    end
    clear data_vec_tmp  
    
    
    data1_vec2=reshape(data1_vec, size(data1_vec,1)*size(data1_vec,2),size(data1_vec,3));
    data2_vec2=reshape(data2_vec, size(data2_vec,1)*size(data2_vec,2),size(data2_vec,3));
    corr12=corr(data1_vec2', data2_vec2', 'type',cfg.corr_type);
 
    % size corr_cue_enc= trials enco x bins enco x trials cue x bins cue
    corr12=  reshape(corr12,size(data1_vec,1),size(data1_vec,2),size(data2_vec,1),size(data2_vec,2));
       
    % fisher z transform correlations
    corr12=  0.5.*log(((ones(size(corr12))+corr12)./(ones(size(corr12))-corr12)));

switch cfg.trial_combi
    case 'in_trial'
  % prepare data: get on diagonal i.e only in trial correlations
  for tr=1:size(data1.trialinfo,1)
     corr_mat(tr,:,:)= squeeze(corr12(tr,:,tr,:));
  end
      corr_trials.trialinfo=data1.trialinfo; % here trialinfos match
    case 'across_trials'
        % keep all possible combinations
     corr_mat(tr,:,:)= permute(corr12(1,3,2,4));  
     corr_trials.trialinfo_dim1=data1.trialinfo; % here trialinfos donot need to match
     corr_trials.trialinfo_dim2=data2.trialinfo; % here trialinfos do not need to match
   
end
                
    corr_trials.corr_mat=corr_mat;
    corr_trials.time1=[t1_1;t2_1];
    corr_trials.time2=[t1_2;t2_2];
    corr_trials.cfg=cfg;

    