function CompVol = CompExtraction(EEG, tw, bltw, fs, ROIChans, isBaselineCorrection, type, admtw)
%{
Args:
    Inputs:
        data (matrix): ChannelNr X Timewindow X Trials
        tw (vector): time window range of interest (second) - zero point:
        stimuli onset
        bttw (vector): baseline twime window range (second)
        fs (int): sampling frequency (Hz)
        ROIChans (cell): channels of interest
        isBaselineCorrection (boolean): conduct baseline correction or not
        type: 'Mean', 'AdaptiveMean', 'Peak'
        admtw (int): adpative time window (second)
    Output: 
        CompVol (vector): return extracted voltage
Author: Chen Danni (dnchen@connect.hku.hk)
Update date: 7/24/2021
%}

%% transform time window range to data point range
bllen = max(abs(bltw));
bltw = (bltw + bllen) * fs;
tw = (tw + bllen) * fs;
admtw = round(admtw/2 * fs);

%% select ROI channels
EEG = pop_select( EEG, 'channel', ROIChans );

%% baseline correction 
data = EEG.data;
data = squeeze(mean(data,1)); % average across ROI channels
bldata = data(int64((bltw(1)+1):bltw(2)),:);
twdata = data(int64((tw(1)+1):tw(2)),:);
blmean = mean(bldata,1); % average time points
if isBaselineCorrection
    twdata = twdata - blmean;
end

%% extract component based on different type
if strcmp (type, 'Mean')
    CompVolMean = mean(twdata,1);
    CompVolMean = squeeze(CompVolMean);
    CompVol = CompVolMean;
elseif strcmp (type, 'AdaptiveMean')
    tmptwdata = abs(twdata);
    [M,I] = max (tmptwdata);
    for j = 1:length(I)
        peakPoint = I(j);
        startPoint = I(j) - admtw;
        endPoint = I(j) + admtw;
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