% smoothes fieldtrip time frequency GA data of a group of subjects with gaussian kernel

%settings: FWHM for time and frequency
%cfg.fwhm_t=time in sec i.e. 0.2
%cfg.fwhm_f=freq in Hz i.e 2

% matlab image processing toolbox needed!


% script by Marie-Christin Fellner mariefellner@gmx.de

function [data_smoothed]=smooth_TF_GAgroup(cfg,data)



if strmatch(data.dimord,'subj_chan_freq_time')==0
    error('data.dimord does not match')
end
    
% settings for smoothing kernel
fwhm_time=cfg.fwhm_t; % in sec
fwhm_freq=cfg.fwhm_f; % in Hz

% get bin sizes
resol_t=data.time(1,2)-data.time(1,1);
resol_f=data.freq(1,2)-data.freq(1,1);

%get fwhm in bins
fwhm_t=fwhm_time/resol_t;
fwhm_f=fwhm_freq/resol_f;


% calculate std based on fwhm (fwhm~=2.355*std)
std_t=fwhm_t/2.355;
std_f=fwhm_f/2.355;

% filt size should at least be 5 times std (if it should be near zero at the edges), needs to be odd
size_t=round(std_t*5);
size_t=size_t+(mod(size_t,1)==0);
size_f=round(std_f*5);
size_f=size_f+(mod(size_f,1)==0);

%  no nans (no smearing of Nans!), data zero padded
nan_ind=find(isnan(data.powspctrm));
data.powspctrm(nan_ind)=0;

% create filters with both std
h_t=fspecial('gaussian',[1,size_t],std_t);
h_f=fspecial('gaussian',[size_f,1],std_f);
h=h_f*h_t;

% to check smoothing kernel plot:
%figure
%surf(h)

powspctrm=data.powspctrm;

%loop over channels and trials
for sub=1:size(powspctrm,1)
for chan=1:size(powspctrm,2)
% use imfilter, conv2 smears data to different tf bins, and is slower!
%w = conv2(squeeze(freq.powspctrm(1,1,:,:)),h); 
powspctrm(sub,chan,:,:)=imfilter(squeeze(powspctrm(sub,chan,:,:)),h,'conv'); 
end

display(strcat ('smoothed sub ',num2str(sub)));
end

% replace with nan_inds with nan
powspctrm(nan_ind)=NaN;
data_smoothed=data;
data_smoothed.powspctrm=powspctrm;
data_smoothed.cfg.smoothingkernel=h;