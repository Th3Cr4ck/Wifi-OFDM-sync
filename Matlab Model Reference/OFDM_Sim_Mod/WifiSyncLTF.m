function [EstVals] = WifiSyncLTF(RxSignal)
%% Parameters
DETECTION_THRESHOLD = 0.72;
MIN_POWER_THRESHOLD = 0.004;
L = 64; % Period for correlation
GI2 = 32; % Guard interval LTF

%% Return values
EstVals = struct;
EstVals.index = 0; % Fin del preambulo e inicio del simbolo OFDM (con CP)
EstVals.syncSucceeded = false;
EstVals.cfo = 0;

if length(RxSignal) < 2*L
    return;
end

delayedRxSignal = delayseq(RxSignal, L); % RxSignal_{m-64}
autocorrSignal = RxSignal .* conj(delayedRxSignal); 

%% Fine CFO estimation
cfoWindow = 160 - L + 1 : 160;
meanAutocorr = mean(autocorrSignal(cfoWindow));

EstVals.cfo = angle(meanAutocorr) / L;
EstVals.syncSucceeded = true;
end