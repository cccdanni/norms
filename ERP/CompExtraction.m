function CompVol = CompExtraction(EEG, tw, bltw, fs, ROIChans, type, admtw_half)
% Args:
%     Inputs:
%         data (matrix): ChannelNr X Timewindow X Trials
%         tw (vector): time window range of interest (second) - zero point:
%         stimuli onset
%         bttw (vector): baseline twime window range (second)
%         fs (int): sampling frequency (Hz)
%         ROIChans (cell): channels of interest
%         type: 'Mean', 'AdaptiveMean', 'Peak'
%         admtw (int): adpative time window (second)
%     Output: 
%         CompVol (vector): return extracted voltage
% Author: Chen Danni (dnchen@connect.hku.hk)
% Update date: 2022-02-28


%% transform time window range to data point range
bllen = max(abs(bltw));
% bltw = (bltw + bllen) * fs;
% tw = (tw + bllen) * fs;
admtw_half = round(admtw_half/2 * fs); 

%% select ROI channels
EEG = pop_select( EEG, 'channel', ROIChans );

%% baseline correction 
data = EEG.data;
data = squeeze(mean(data,1)); % average across ROI channels

times = EEG.times;
tw_seconds = tw * 1000;
idx_start  = nearest(times, tw_seconds(1));
idx_end    = nearest(times, tw_seconds(2));

twdata = data(idx_start:idx_end,:);

%% extract component based on different type
if strcmp (type, 'Mean')
    CompVolMean = mean(twdata,1);
    CompVolMean = squeeze(CompVolMean);
    CompVol = CompVolMean;
elseif strcmp (type, 'AdaptiveMean')
    tmptwdata = abs(twdata); % Window of Interest * nTrial
    [M,I] = max (tmptwdata); % Find the peak within the corresponding time window; M - maximum; I - the indices
    for j = 1:length(I) % loop for each trial
        peakPoint = I(j);
        startPoint = I(j) - admtw_half;
        endPoint = I(j) + admtw_half;
        if startPoint <= 0
            startPoint = 1;
        end
        if endPoint > size(tmptwdata,1)
            endPoint = size(tmptwdata,1);
        end
        thisTrial = twdata((startPoint:endPoint),j);
        CompVolAdpMean(j) = mean(thisTrial);
    end
    CompVol = CompVolAdpMean;
elseif strcmp (type, 'Peak')
    tmptwdata = abs(twdata);
    [M,I] = max (tmptwdata);
    for j = 1:length(I)
        CompVolPeak(j) = twdata(I(j), j);
    end
    CompVol = CompVolPeak;
end


end