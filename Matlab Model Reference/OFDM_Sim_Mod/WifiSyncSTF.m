function [EstVals, syncMetric, diffSyncM] = WifiSyncSTF(RxSignal, delayFactor)
%% Parameters
L = 16;
PLATEAU_THRESHOLD = 0.30;
DIFF_THRESHOLD = -0.20;
MIN_POWER_THRESHOLD = 0.0004;

%% Return values
EstVals = struct;
EstVals.index = 0; % Index of LTF start
EstVals.syncSucceeded = false;
EstVals.cfo = 0; % rad/sample

%% Metrics
if length(RxSignal) < 2*L
    return;
end

delayedRxSignal = delayseq(RxSignal, L); % RxSignal_{m-16}
autocorrSignal = RxSignal .* conj(delayedRxSignal); 
powerSignal = abs(RxSignal).^2;

corrFilter = ones(L,1);
autocorrMetric = filter(corrFilter, 1, autocorrSignal);
powerMetric = filter(corrFilter, 1, powerSignal);
powerMetricDelayed = delayseq(powerMetric, L);

syncMetric = abs(autocorrMetric).^2 ./ ...
    max(powerMetric .* powerMetricDelayed, MIN_POWER_THRESHOLD);

%% Differences
delayedSyncM = delayseq(syncMetric, delayFactor); % M{D-delayFactor}
diffSyncM = syncMetric - delayedSyncM;

%% Detect STF
plateauMask = syncMetric > PLATEAU_THRESHOLD;
% Implementar filtro. Al menos una ventana grande antes de superar el
% umbral.

diffMask = diffSyncM <= DIFF_THRESHOLD;

windowLen = 2*L;
longPlateauWindowMask = zeros(1, length(syncMetric));
plateauCount = 0;
for i = 1:length(syncMetric)
    if plateauMask(i) && plateauCount<windowLen 
        plateauCount = plateauCount+1; 
    elseif ~plateauMask(i) && plateauCount>0
        plateauCount = plateauCount-1;
    end
    longPlateauWindowMask(i) = plateauCount >= windowLen-3;
end

candidatesMask = longPlateauWindowMask & diffMask;

if ~any(candidatesMask)
    return;
end

% Best index for FFT
candidatesIndexes = find(candidatesMask);
EstVals.index = candidatesIndexes(1)-(delayFactor-1);

% Ensure enough samples remain for CFO estimation
if EstVals.index - L + 1 < 1
    return;
end

%% Coarse CFO estimation
cfoWindow = EstVals.index - L + 1: EstVals.index;
meanAutocorr = mean(autocorrSignal(cfoWindow));

EstVals.cfo = angle(meanAutocorr) / L;
EstVals.syncSucceeded = true;

end