function [EstVals] = OfdmSync(RxSignal,L,CP)
%% Parameters
PLATEAU_THRESHOLD = 0.01;
%RxSignalSinc = 0;
DELAY_FACTOR = 8;
MAJORITY_WINDOW_SIZE = 7;
DETECTION_THRESHOLD = 0.67;
MIN_POWER_THRESHOLD = 0.004;

%% Return values
EstVals = struct;
EstVals.index = L*2 + CP + 4;
EstVals.syncSucceeded = false;
EstVals.cfo = 0;

%% Sincronía
delayedRxSignal = delayseq(RxSignal,L);
autocorrProducts = conj(RxSignal) .* delayedRxSignal;
signalPower = conj(RxSignal) .* RxSignal;

movingAverage = dsp.MovingAverage(L);
autocorrMetric = movingAverage(autocorrProducts);
movingAverage.reset;
powerMetric = movingAverage(signalPower);

syncMetric = (abs(autocorrMetric).^2) ./ max(powerMetric.^2, MIN_POWER_THRESHOLD);
%figure(10); plot(syncMetric)

detectionMask = syncMetric > DETECTION_THRESHOLD;
if(sum(detectionMask))
    [~,initialDetectionIndex] = max(detectionMask ~= 0,[],1);
    if (initialDetectionIndex < 2*L + CP)
        % Ignore detections before the earliest possible preamble location
        detectionMask(1:2*L + CP) = 0;
        [~,initialDetectionIndex] = max(detectionMask ~= 0,[],1);
    end
    if (initialDetectionIndex + L > length(RxSignal))
        return
    end

    % Detect the flat region (plateau) of the S&C metric
    delayedSyncMetric = delayseq(syncMetric,DELAY_FACTOR);
    syncMetric = syncMetric(initialDetectionIndex:initialDetectionIndex + L);
    delayedSyncMetric = delayedSyncMetric(initialDetectionIndex:initialDetectionIndex + L);
    autocorrMetric = autocorrMetric(initialDetectionIndex:initialDetectionIndex + L);

    metricDifference = syncMetric - delayedSyncMetric;

    plateauMask = (metricDifference < PLATEAU_THRESHOLD) & (metricDifference > -(2*PLATEAU_THRESHOLD));

    plateauStartIndex = binmodeFilter(plateauMask,MAJORITY_WINDOW_SIZE);
    if ((plateauStartIndex == 0) || ...
        ((plateauStartIndex + 1 + floor(0.55*CP) + 2) > L))
        return
    end
    %[~,plateauStartIndex] = max(plateauMask);
    estimatedSyncIndex = plateauStartIndex - DELAY_FACTOR + initialDetectionIndex + floor(0.7*CP);

    EstVals.index = estimatedSyncIndex;

    estimatedCFO = -angle(autocorrMetric(plateauStartIndex + 1 + floor(0.55*CP))) / L;
    EstVals.cfo = estimatedCFO;
    EstVals.syncSucceeded = true;
end
end

function index = binmodeFilter(x,GSize)
% Majority filter: checks whether most samples are true in a local window.
index = int32(0);

for i = 1:numel(x)
    % Define the local window around the current sample.
    localWindow = x(max(i-(GSize+1)/2,1):min(i+(GSize+1)/2,numel(x)));
    majorityDetected = sum(localWindow) > (GSize-1)/2;

    % Return the first position where the majority condition is satisfied.
    if (index < 1) && majorityDetected
        index = i;
        break;
    end
end
end