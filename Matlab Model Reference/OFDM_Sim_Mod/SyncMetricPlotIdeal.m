clc;
clear;
close all;

%% Parameters
Fs = 10e6;
Fc = 2.4e9;
Mary = 4;   % QAM order
M = 20;     % OFDM Symbols by frame
Nfft = 64;
GI = 16;
GI2 = 32;
initialDelay = 200;  % Initial zeros
SNR_dB = 20;% 0:5:50;
it = 100;

% syncFailedSTF = false(length(SNR_dB), it);
% syncFailedLTF = false(length(SNR_dB), it);
syncMetricSTF = [];
syncMetricLTF = [];
diffSyncMetricSTF  = [];
diffSyncMetricLTF  = [];

delayFactor = [1 2 3 4 5 6 7 8 9 10 11 12 13];

ACTIVE_TDL_CHANNEL = false;

%% Channel Model
if ACTIVE_TDL_CHANNEL              
    v_kmh = 100.0;                % UE velocity in km/h
    c = physconst('lightspeed');  % speed of light in m/s
    v_ms = v_kmh * 1000/3600;     % convert to m/s
    Fd_max = (v_ms/c) * Fc;       % maximum Doppler frequency in Hz
    
    tdl = nrTDLChannel;
    tdl.DelayProfile = 'TDL-A';
    tdl.DelaySpread = 400e-9;
    tdl.MaximumDopplerShift = Fd_max;
    tdl.SampleRate = Fs;
    tdl.NumReceiveAntennas = 1;
end


%% Generate preamble
preamShort = generatePreambleSTF();
preamLong = generatePreambleLTF();

% Normalize power
preamPower = 1;
scaleFactorShort = sqrt(preamPower) / sqrt(mean(abs(preamShort).^2));
scaleFactorLong = sqrt(preamPower) / sqrt(mean(abs(preamLong).^2));

preamShort = preamShort * scaleFactorShort;
preamLong = preamLong * scaleFactorLong;

%for SNR = SNR_dB   
for delay = delayFactor
    for i = 1:it
        %% Generate payload
        % 52 active subcarriers, excluding DC
        activeSubcarriers = [-26:-1 1:26];
        
        transmittedBits = round(rand(52*log2(Mary), 1));
        payloadSymbols = qammod(transmittedBits, Mary, ...
            'InputType', 'bit', 'UnitAveragePower', true).';
           
        X = zeros(1,Nfft);
        activeBins = mod(activeSubcarriers, Nfft) + 1;
        X(activeBins) = payloadSymbols;
        
        payloadUseful = ifft(X);               
                   
        % Normalize power
        payloadScale = 1 / sqrt(mean(abs(payloadUseful).^2));
        payloadUseful = payloadUseful * payloadScale;        

        % Add CP
        payload = [payloadUseful(end-GI+1:end), payloadUseful];                             

        %% Assemble transmitted frame
        txFrame = [zeros(1,initialDelay), ...
           preamShort, ...
           preamLong, ...
           payload, ...
           zeros(1, initialDelay/2)];
        
       %% TDL Channel
        if ACTIVE_TDL_CHANNEL                          
            [rxFrame, pathGains] = tdl(transpose(txFrame));
            offsetChannel = nrPerfectTimingEstimate(pathGains, getPathFilters(tdl));
            trueSTFStart = initialDelay + 1 + offsetChannel;
            trueLTFStart = initialDelay + length(preamShort) + 1 + offsetChannel;
            rxFrame = transpose(rxFrame);
        else
            offsetChannel = 0;
            trueSTFStart = initialDelay + 1;
            trueLTFStart = initialDelay + length(preamShort) + 1;
            rxFrame = txFrame;
        end
     

        %% Add noise
        % activeSignal = rxFrame(initialDelay+1:end-initialDelay/2);
        % signalPower = mean(abs(activeSignal).^2);
        % noisePower = signalPower / 10^(SNR_dB/10);
        % 
        % noise = sqrt(noisePower/2) * ...
        %     (randn(1,length(rxFrame)) + 1j*randn(1, length(rxFrame)));
        % 
        % rxFrame = rxFrame + noise;
        % 
        %% Synchronization
        [EstSTF, syncMSTF, diffSyncMSTF] = WifiSyncSTF(rxFrame, delay);
        
        if isempty(syncMetricSTF)            
            syncMetricLenSTF = length(syncMSTF);
            diffSyncMLenSTF = length(diffSyncMSTF);
            syncMetricSTF = zeros(length(delayFactor), it, syncMetricLenSTF);
            diffSyncMetricSTF = zeros(length(delayFactor), it, diffSyncMLenSTF);
        end
        syncMetricSTF(delay==delayFactor, i, :) = syncMSTF;
        diffSyncMetricSTF(delay==delayFactor, i, :) = diffSyncMSTF;


        [EstLTF, syncMLTF, diffSyncMLTF] = ...
            WifiSyncLTF(rxFrame(trueLTFStart:trueLTFStart+260), delay);

        if isempty(syncMetricLTF)            
            syncMetricLenLTF = length(syncMLTF);
            diffSyncMLenLTF = length(diffSyncMLTF);
            % syncMetricLTF = zeros(length(SNR_dB), it, syncMetricLenLTF);
            % diffSyncMetricLTF = zeros(length(SNR_dB), it, diffSyncMLenLTF);
            syncMetricLTF = zeros(length(delayFactor), it, syncMetricLenLTF);
            diffSyncMetricLTF = zeros(length(delayFactor), it, diffSyncMLenLTF);
        end
        % syncMetricLTF(SNR==SNR_dB, i, :) = syncMLTF;
        % diffSyncMetricLTF(SNR==SNR_dB, i, :) = diffSyncMLTF;
        syncMetricLTF(delay==delayFactor, i, :) = syncMLTF;
        diffSyncMetricLTF(delay==delayFactor, i, :) = diffSyncMLTF;

   
        truePayloadStart = trueLTFStart + length(preamLong);
        % payloadStart = EstSTF.index + EstLTF.index - 1;
        % 
        % % Sync index must be a valid point for the fft
        % if payloadStart < truePayloadStart - GI + 1 ...
        %         || payloadStart > truePayloadStart + 1
        %     syncFailedLTF(SNR == SNR_dB, i) = true;
        %     continue;
        % end
            
    end
end

%% Results 

meanSyncMetric = squeeze(mean(syncMetricSTF, 2));
meanDiffSyncM = squeeze(mean(diffSyncMetricSTF, 2));
colors = lines(length(SNR_dB));
colors = lines(length(delayFactor));

figure;
hold on; grid on;
for k = 1:length(delayFactor)
    % name = sprintf('SNR=%2d dB', SNR_dB(k));
    name = sprintf('Delay Factor=%d', delayFactor(k));

    plot(meanSyncMetric(k,:), '-', ...
        'Color', colors(k,:), ...
        'DisplayName', name, ...
        'LineWidth', 3);

    plot(meanDiffSyncM(k,:), '--', ...
        'Color', colors(k,:), ...
        'HandleVisibility', 'off', ...
        'LineWidth', 1);    
end

xline(trueLTFStart, 'k--', 'LTF Start', ...
    'HandleVisibility', 'off', 'LineWidth', 1);

% legend('Location', 'northeast');
xlabel('Samples');
ylabel('M(d)');
title('STF Sync Metric');

meanSyncMetric = squeeze(mean(syncMetricLTF, 2));
meanDiffSyncM = squeeze(mean(diffSyncMetricLTF, 2));

figure; 
hold on; grid on;
for k = 1:length(delayFactor)
    % name = sprintf('SNR=%2d dB', SNR_dB(k));
    name = sprintf('Delay Factor=%d', delayFactor(k));

    plot(meanSyncMetric(k,:), '-', ...
        'Color', colors(k,:), ...
        'DisplayName', name, ...
        'LineWidth', 3);

    plot(meanDiffSyncM(k,:), '--', ...
        'Color', colors(k,:), ...
        'HandleVisibility', 'off', ...
        'LineWidth', 1);
end

xline(truePayloadStart-trueLTFStart, 'k--', 'Payload Start', ...
    'HandleVisibility', 'off', 'LineWidth', 1);

% legend('Location', 'northeast');
xlabel('Samples');
ylabel('M(d)');
title('LTF Sync Metric');


%% Only diff figures
% One figure with multiple delay factors
colors = lines(length(delayFactor));
figure; hold on;
name = sprintf('STF Diff Sync Metric');
for k = 1:length(delayFactor)
    leg = sprintf('K = %d', k);
    plot(squeeze(diffSyncMetricSTF(k,1,:)), '-', ...
    'Color', colors(k,:), ...
    'LineWidth', 2, ...
    'DisplayName', leg);
end
xline(trueLTFStart, '--', 'LTF Start', ...
'HandleVisibility', 'off', 'LineWidth', 1);
xlim([trueLTFStart-2 trueLTFStart+15]);
xlabel('Samples');
ylabel('M(d)');
grid on;
title(name);
legend;

%% Figure per delay factor
% colors = lines(length(delayFactor));
% for k = 1:length(delayFactor)
%     figure; hold on;
%     name = sprintf('STF Diff Sync Metric, Delay Factor = %d', delayFactor(k));
% 
%     for l = 1:7    
%     plot(squeeze(diffSyncMetricSTF(k,l,:)), '-', ...
%         'Color', colors(l,:), ...
%         'LineWidth', 2);  
%     end
% 
%     xline(trueLTFStart, '--', 'LTF Start', ...
%     'HandleVisibility', 'off', 'LineWidth', 1);
% 
%     xlim([trueLTFStart-2 trueLTFStart+15]);
%     xlabel('Samples');
%     ylabel('M(d)');
%     grid on;
%     title(name);
% end

% % LTF
% for k = 1:length(delayFactor)
%     figure; hold on;
%     name = sprintf('LTF Diff Sync Metric, Delay Factor = %d', delayFactor(k));
% 
%     for l = 10:17    
%     plot(squeeze(diffSyncMetricLTF(k,l,:)), '-', ...
%         'Color', colors(l-9,:), ...
%         'LineWidth', 2);  
%     end
% 
%     xline(truePayloadStart-trueLTFStart, 'k--', 'Payload Start', ...
%     'HandleVisibility', 'off', 'LineWidth', 1);
% 
%     xlim([158 170]);
%     xlabel('Samples');
%     ylabel('M(d)');
%     grid on;
%     title(name);
% end

%% Variance
v = squeeze(diffSyncMetricSTF(3, :, trueLTFStart-2:trueLTFStart+15));
varDiffSTF = mean(var(v, 0, 1));


shortPayloadStart = truePayloadStart-trueLTFStart;
v = squeeze(diffSyncMetricLTF(3, :, shortPayloadStart-2:shortPayloadStart+15));
varDiffLTF = mean(var(v, 0, 1));

fprintf('Variance Diff STF = %0.5f\n', sqrt(varDiffSTF));
fprintf('Variance Diff LTF = %0.5f\n', sqrt(varDiffLTF));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
function preamLong = generatePreambleLTF()
    k = -26:26;
    L = [1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 ...
         0 ...
         1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1];
    
    Lfft = zeros(1,64);
    Lfft(mod(k,64)+1) = L;          % k negativo -> bins 64+k
    
    l64  = ifft(Lfft);
    GI2 = 32;
    preamLong = [l64(GI2+1:64) l64 l64];    % GI2 + T1 + T2 = 160 muestras
end

function preamShort = generatePreambleSTF()
    k    = [-24 -20 -16 -12 -8 -4 4 8 12 16 20 24];
    vals = (1+1j) * [1 -1 1 -1 -1 1 -1 -1 1 1 1 1];
    
    S = zeros(1,64);
    S(mod(k,64)+1) = sqrt(13/6) * vals;   % mapeo a bins FFT (k negativo -> 64+k)
    
    s64   = ifft(S);                       % 64 muestras, periodo de 16
    preamShort = repmat(s64(1:16), 1, 10);      % 160 muestras = 8 µs
end

function [EstVals, syncMetric, diffSyncM] = WifiSyncLTF(RxSignal, DELAY_FACTOR)
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
powerSignal = abs(RxSignal).^2;

corrFilter = ones(L,1);
autocorrMetric = filter(corrFilter, 1, autocorrSignal);
powerMetric = filter(corrFilter, 1, powerSignal);
delayedPowerMetric = delayseq(powerMetric, L);

syncMetric = abs(autocorrMetric).^2 ./ ...
    max(powerMetric .* delayedPowerMetric, MIN_POWER_THRESHOLD);

%% Detect LTF
% detectionMask = syncMetric > DETECTION_THRESHOLD;
% plateauIndexes = find(detectionMask);
% 
% if isempty(plateauIndexes)
%     return;
% end
% 
% initialDetectionIndex = plateauIndexes(1);
% lastDetectionIndex = plateauIndexes(end);
% 
% % Ensure enough samples remain for CFO estimation
% if initialDetectionIndex + L - 1 > length(RxSignal)
%     return;
% end

%% Differences

delayedSyncM = delayseq(syncMetric, DELAY_FACTOR); % M{D-6}
diffSyncM = syncMetric - delayedSyncM;
%mask = detectionMask & 
% index = min(delayedSyncM);
% negDiffVals = diffSyncM(mask);
% negDiffIndexes = find(mask);
return;


%% Fine CFO estimation
cfoWindow = initialDetectionIndex : initialDetectionIndex + L - 1;
meanAutocorr = mean(autocorrMetric(cfoWindow));

EstVals.index = floor((initialDetectionIndex + lastDetectionIndex)/2);
EstVals.index = lastDetectionIndex - round(0.22 * GI2);
EstVals.cfo = angle(meanAutocorr) / L;
EstVals.syncSucceeded = true;
end

function [EstVals, syncMetric, diffSyncM] = WifiSyncSTF(RxSignal, DELAY_FACTOR)
%% Parameters
L = 16;
DETECTION_THRESHOLD = 0.72;
MIN_POWER_THRESHOLD = 0.0004;

%% Return values
EstVals = struct;
EstVals.index = 0; % Index of LTF start
EstVals.syncSucceeded = false;
EstVals.cfo = 0; % rad/sample

%% STF detection
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
% plot(syncMetric);


%% Detect STF
% detectionMask = syncMetric > DETECTION_THRESHOLD;
% plateauIndexes = find(detectionMask);
% 
% if isempty(plateauIndexes)
%     return;
% end
% 
% initialDetectionIndex = plateauIndexes(1);
% lastDetectionIndex = plateauIndexes(end);
% 
% % Ensure enough samples remain for CFO estimation
% if lastDetectionIndex - L + 1 < 0
%     return;
% end

%% Differences
delayedSyncM = delayseq(syncMetric, DELAY_FACTOR); % M{D-6}
diffSyncM = syncMetric - delayedSyncM;
return;

%% Best index for FFT
EstVals.index = ceil((lastDetectionIndex + initialDetectionIndex)/2);

%% Coarse CFO estimation
cfoWindow = lastDetectionIndex - L + 1: lastDetectionIndex;
meanAutocorr = mean(autocorrMetric(cfoWindow));

EstVals.cfo = angle(meanAutocorr) / L;
EstVals.syncSucceeded = true;

end