% z_transformation for fieldtrip timefrequency data 

% important: use only on all conditions concatenated! (otherwise you will
% cancel out you effects of interest)
%cfg.time=[tmin tmax]


% script by Marie-Christin Fellner mariefellner@gmx.de

function [z_freq]=z_trans_TF_seltime(cfg,data)
            switch data.dimord
                case 'rpt_chan_freq_time'
                    trial_dim=1;
                otherwise
                    error('fix dimension order')               
            end
            
            
            t1=nearest(data.time,cfg.time(1));
            t2=nearest(data.time,cfg.time(2));
            

           powspctrm=data.powspctrm(:,:,:,t1:t2);
           
           % concatenate single trials
           pow=permute(powspctrm,[4,1,2,3]);
           pow=reshape(pow,size(pow,1)*size(pow,2),size(pow,3),size(pow,4));
           
           % compute mean/std for each chan/freq
           mean_pow=squeeze(nanmean(pow));
           std_pow=squeeze(nanstd(pow));
            clear pow
           
           % reshape mean/std for matrix based calculations
            powspctrm=data.powspctrm;
           mean_pow=repmat(reshape(mean_pow,[1,size(mean_pow),1]),[size(powspctrm,1),1,1,size(powspctrm,4)]);
           std_pow=repmat(reshape(std_pow,[1,size(std_pow),1]),[size(powspctrm,1),1,1,size(powspctrm,4)]);

           
           % z-trans pow (pow-mean/std)

           powspctrm=(powspctrm-mean_pow)./std_pow;
           
           z_freq=data;
           z_freq.powspctrm=powspctrm;
           z_freq.cfg.baseline='z_trans using z_trans_TF';