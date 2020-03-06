function [fixationtimes,saccadetimes,saccadeIdx,info] = adaptive_VT_mad(eyedat,samplingFreq,useMAD,exciseIntersaccade,initVel,lambda,fixedThreshVal)
%function [fixationtimes,saccadetimes,saccadeIdx,info] = adaptive_VT_mad(eyedat,samplingFreq,useMAD,exciseIntersaccade,initVel,lambda,fixedThreshVal)
%
% Adaptive threshold algorithm. 
%
% adapted from:
% Nyström, M., & Holmqvist, K. (2010). An adaptive algorithm for fixation, saccade, and glissade de-tection in eyetracking data. Behavior Research Methods, 42(1), 188?204. 
%
% original code: http://dev.humlab.lu.se/www-transfer/people/marcus-nystrom/EventDetector1.1.zip
%
% SDK, 2017: re-written for looking at recording data since wasn't handling NaNs wells
% BV 2018: updated to use MAD to estimate threshold

info = [];

if isempty(fixedThreshVal)
    useFixedThresh = 0;
else
    useFixedThresh = 1;
end

%---Default Parameters---%
X = eyedat(1,:);
Y = eyedat(2,:);
minFixDur = 0.040; % in seconds
minSaccadeDur = 0.010; % in seconds

%---Filter Data and Calculate Velocity/ederation of Eye---%

% Lowpass filter window length
smoothInt = minSaccadeDur; % in seconds

% Span of filter
span = ceil(smoothInt*samplingFreq);

% Calculate unfiltered data
Xorg = X;
Yorg = Y;

velXorg = [0 diff(X)]/samplingFreq; %removed angleInPixelsV
velYorg = [0 diff(Y)]/samplingFreq; %remvoed angleInPixelsV
velOrg = sqrt(velXorg.^2 + velYorg.^2);

% Pixel values, velocities, and accelerations
N = 2;                 % Order of polynomial fit
F = 2*ceil(span)-1;    % Window length
[b,g] = sgolay(N,F);   % Calculate S-G coefficients
Nf = F;

% Calculate the velocity and acceleration
tempX = conv(X, g(:,1)', 'same');
tempX(1:(Nf-1)/2) = X(1:(Nf-1)/2);
tempX(end-(Nf-3)/2:end) = X(end-(Nf-3)/2:end);
tempY = conv(Y, g(:,1)', 'same');
tempY(1:(Nf-1)/2) = Y(1:(Nf-1)/2);
tempY(end-(Nf-3)/2:end) = Y(end-(Nf-3)/2:end);
X = tempX;
Y = tempY;

tempX = conv(X, -g(:,2)', 'same');
tempX(1:(Nf-1)/2) = 0;
tempX(end-(Nf-3)/2:end) = 0;
tempY = conv(Y, -g(:,2)', 'same');
tempY(1:(Nf-1)/2) = 0;
tempY(end-(Nf-3)/2:end) = 0;
velX = tempX;
velY = tempY;
vel = sqrt(velX.^2 + velY.^2)*samplingFreq;

tempX = conv(X, g(:,3)', 'same');
tempX(1:(Nf-1)/2) = 0;
tempX(end-(Nf-3)/2:end) = 0;
tempY = conv(Y, g(:,3)', 'same');
tempY(1:(Nf-1)/2) = 0;
tempY(end-(Nf-3)/2:end) = 0;
accX = tempX;
accY = tempY;
acc = sqrt(accX.^2 + accY.^2)*samplingFreq^2;

info.vel = vel;
info.acc = acc;

%---Remove Blinks---%
nanIdx = zeros(1,length(X));
nanIdx(isnan(X)) = 1;

%---Detect Fixations and Sacccades---%
oldPeak  = inf;
peakDetectionThreshold = initVel;
niter = 0;

if ~useFixedThresh
    threshHist = peakDetectionThreshold;
    while abs(peakDetectionThreshold -  oldPeak) > 1
        niter=niter+1;
        oldPeak = peakDetectionThreshold;

        tempcross = zeros(1,length(X));
        InitialVelPeakIdx = find(vel > peakDetectionThreshold);
        tempcross(InitialVelPeakIdx) = 1;

        possibleFixationIdx = ~tempcross;
        fixLabeled = bwlabel(possibleFixationIdx);

        % Process one inter-peak-saccadic periods (called fixations below,
        % although they are not identified as fixations yet).
        fixNoise = [];
        for k = 1:max(fixLabeled)

            % The samples related to the current fixation
            fixIdx = find(fixLabeled == k);

            % Check that the fixation duration exceeds the minimum duration criteria.
            if length(fixIdx)/samplingFreq < minFixDur
                continue
            end

            % Extract the samples from the center of the fixation
            if exciseIntersaccade
                centralFixSamples = minFixDur*samplingFreq/6;
            else
                centralFixSamples = 0;
            end
            fNoise = vel(floor(fixIdx(1)+centralFixSamples):ceil(fixIdx(end)-centralFixSamples));
            fixNoise = [fixNoise fNoise];
        end

        % calcualte average and std with median/MAD?
        if useMAD
            avgNoise = nanmedian(fixNoise);
            stdNoise = 1.4826 * mad(fixNoise,1); 
        else
            avgNoise = nanmean(fixNoise);
            stdNoise = nanstd(fixNoise);
        end

        % Base the peak velocity threshold on the noise level
        peakDetectionThreshold =  avgNoise + lambda*stdNoise;
        saccadeVelocityTreshold = avgNoise + 3*stdNoise;
        velPeakIdx  = vel > peakDetectionThreshold;

        threshHist(numel(threshHist)+1,1) = peakDetectionThreshold;
    end
else
    peakDetectionThreshold = fixedThreshVal;
    saccadeVelocityTreshold = fixedThreshVal - 10;
    velPeakIdx  = vel > peakDetectionThreshold;
    threshHist = peakDetectionThreshold;
end


info.threshHist = threshHist;
info.peakDetectionThreshold = peakDetectionThreshold;
info.saccadeVelocityTreshold = saccadeVelocityTreshold;
info.niter = niter;

%%
%---Detect Saccades---%
len = length(vel);

% Preallocate memory
velLabeled = bwlabel(velPeakIdx);
saccadeIdx = zeros(1,len);    % Saccade index

% If no saccades are detected, return
if isempty(velLabeled)
    return
end

% Process one velocity peak at the time
kk = 1;

for k = 1:max(velLabeled) 

    %----------------------------------------------------------------------
    % Check the saccade peak samples
    %----------------------------------------------------------------------
    % The samples related to the current saccade
    peakIdx = find(velLabeled == k);
    
    
    % If the peak consists of =< minPeakSamples consequtive samples, it it probably
    % noise (1/6 or the min saccade duration)
    minPeakSamples = ceil(minSaccadeDur/6*samplingFreq);
    if length(peakIdx) <= minPeakSamples, continue, end
    
    % Check whether this peak is already included in the previous saccade
    % (can be like this for glissades)
    if kk > 1
        if ~isempty(intersect(peakIdx,[find(saccadeIdx) find(glissadeIdx)]))
            continue
        end
    end
    
    %----------------------------------------------------------------------
    % DETECT SACCADE
    %----------------------------------------------------------------------
    
    % Detect saccade start.  AND acc <= 0
    saccadeStartIdx = find(vel(peakIdx(1):-1:1) <= saccadeVelocityTreshold &...% vel <= global vel threshold
        [diff(vel(peakIdx(1):-1:1)) 0] >= 0);          % acc <= 0
    if isempty(saccadeStartIdx), continue, end
    saccadeStartIdx = peakIdx(1) - saccadeStartIdx(1) + 1;
    
    % Calculate local fixation noise (the adaptive part)
    localVelNoise = vel(saccadeStartIdx:-1: max(1,ceil(saccadeStartIdx - minFixDur*samplingFreq)));
    localVelNoise = mean(localVelNoise) + 3*std(localVelNoise);
    localsaccadeVelocityTreshold = localVelNoise*0.3 + saccadeVelocityTreshold*0.7; % 30% local + 70% global
    
    % Check whether the local vel. noise exceeds the peak vel. threshold.
    if localVelNoise > peakDetectionThreshold, continue, end
    
    % Detect end of saccade (without glissade)
    saccadeEndIdx = find(vel(peakIdx(end):end) <= localsaccadeVelocityTreshold &...             % vel <= adaptive vel threshold
        [diff(vel(peakIdx(end):end)) 0] >= 0);        % acc <= 0
    
    if isempty(saccadeEndIdx), continue, end
    saccadeEndIdx = peakIdx(end) + saccadeEndIdx(1) - 1;
    
    % If the saccade contains NaN samples, continue
    if any(nanIdx(saccadeStartIdx:saccadeEndIdx)), continue, end
    
    % Make sure the saccade duration exceeds the minimum duration.
    saccadeLen = saccadeEndIdx - saccadeStartIdx;
    if saccadeLen/samplingFreq < minSaccadeDur
        continue
    end
    
    % If all the above criteria are fulfilled, label it as a saccade.
    saccadeIdx(saccadeStartIdx:saccadeEndIdx) = 1;
    localSaccadeVelocityTreshold = localsaccadeVelocityTreshold;
    
    % Collect information about the saccade
    saccadeInfo.start = saccadeStartIdx/samplingFreq; % in ms
    saccadeInfo.end = saccadeEndIdx/samplingFreq; % in ms
    saccadeInfo.duration = saccadeInfo.end - saccadeInfo.start;
    saccadeInfo.amplitude = sqrt(((X(saccadeEndIdx)-...
        (X(saccadeStartIdx))))^2 + ...
        ((Y(saccadeEndIdx)-...
        (Y(saccadeStartIdx))))^2   );
    saccadeInfo.peakVelocity = max(vel(saccadeStartIdx:saccadeEndIdx));
    saccadeInfo.peakAcceleration = max(acc(saccadeStartIdx:saccadeEndIdx));
    
    
end


%use format from other code
[saccadetimes] = BehavioralIndex(find(saccadeIdx));

%---Detect Fixations---%
%going to do it in a comparable manner to Cluster Fix so that this part
%doesn't change the methodology

fixationindexes = ones(1,len);
fixationindexes(nanIdx == 1) = 0;
fixationindexes(saccadeIdx ==1) = 0;
[fixationtimes] = BehavioralIndex(find(fixationindexes));
fixdur = diff(fixationtimes,1)+1;
too_short = find(fixdur < minFixDur*samplingFreq);
fixationtimes(:,too_short) = [];

end

function [behaviortime] = BehavioralIndex(behavind)
[behaveind]=findgaps(behavind);
if isempty(behavind)
    behaviortime = [];
    return
end
for index=1:size(behaveind,1)
    rowfixind = behaveind(index,:);
    rowfixind(rowfixind == 0) = [];
    behaviortime(:,index) = [rowfixind(1);rowfixind(end)];
end
end

function [broken_ind]=findgaps(input_ind)
% finds gaps (greater than 0) in between indeces in a vector
%
% rechecked for bugs by SDK on 1/5/2017

gaps =find(abs(diff(input_ind)) > 1);
broken_ind = zeros(length(gaps),50);
if ~isempty(gaps)
    for gapind = 1:length(gaps)+1;
        if gapind == 1;
            temp = input_ind(1:gaps(gapind));
        elseif gapind == length(gaps)+1
            temp = input_ind(gaps(gapind-1)+1:end);
        else
            temp = input_ind(gaps(gapind-1)+1:gaps(gapind));
        end
        broken_ind(gapind,1:length(temp)) = temp;
    end
else
    broken_ind = input_ind;
end

end
