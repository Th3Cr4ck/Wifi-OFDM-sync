clear;
close all;
clc;

%% Parameters
Fs = 10e6;
Fc = 2.4e9;
Nfft = 64;
Mary = 64;
GI = 16;
SNR_dB = 5:5:50;
iter = 1000;
initialDelay = 200;  % Initial zeros

CFO_norm_values = [0.1 0.2 0.5 0.8 1.2 1.6 1.9]; % CFO normzalized by subcarrier spacing
%CFO_norm_values = [0.0064 0.0032 0.0016 0.00064 0.00032 0.00016 0.000064];
trueCFO_Hz = CFO_norm_values*Fs/Nfft;
trueCFO = 2*pi*trueCFO_Hz/Fs;  % rad/sample

ACTIVE_TDL_CHANNEL = true;

fullSize = size(zeros(length(CFO_norm_values), length(SNR_dB), iter));
coarseCFO = zeros(fullSize);
totalCFO = zeros(fullSize);
validSynch = false(fullSize);
BERFullSync = zeros(fullSize);
BERHalfSync = zeros(fullSize);

%% Generate preamble
preamShort = generatePreambleSTF();
preamLong = generatePreambleLTF();

% Normalize power
preamShort = preamShort / sqrt(mean(abs(preamShort).^2));
preamLong = preamLong / sqrt(mean(abs(preamLong).^2));

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


for m = 1:length(CFO_norm_values)
    for l = 1:length(SNR_dB)
        lSNR = SNR_dB(l);
        for i = 1:iter
            %% Generate payload
            % 52 active subcarriers, excluding DC
            activeSubcarriers = [-26:-1 1:26];
            
            txBits = round(rand(52*log2(Mary), 1));
            payloadSymbols = qammod(txBits, Mary, ...
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
            
            %% Apply CFO
            n = 0:length(rxFrame)-1;
            rxFrame = rxFrame .* exp(1j*trueCFO(m)*n);
            
            %% AGWN Noise
            activeSignal = rxFrame(initialDelay+1:end-initialDelay/2);
            signalPower = mean(abs(activeSignal).^2);
            noisePower = signalPower / 10^(lSNR/10);
            
            noise = sqrt(noisePower/2) * ...
                (randn(size(rxFrame)) + 1j*randn(size(rxFrame)));
            
            rxFrame = rxFrame + noise;
            
            %% Run STF synchronization            
            EstSTF = WifiSyncSTF(rxFrame, 4);
                        
            if ~EstSTF.syncSucceeded
                continue;
            end
            
            % Correct coarse CFO
            coarseCFO(m, l, i) = EstSTF.cfo*Fs/(2*pi);            
            rxCoarseCorrected = rxFrame .* exp(-1j*EstSTF.cfo*n);            
            
            %% Run LTF synchronization               
            EstLTF = WifiSyncLTF(transpose(rxCoarseCorrected(EstSTF.index:end)));

            if ~EstLTF.syncSucceeded
                continue;                             
            end
    
            validSynch(m, l, i) = true;

            % intCFO_norm = round((EstSTF.cfo - EstLTF.cfo)/(2*pi/Nfft));
            % estCFO = EstLTF.cfo + intCFO_norm * (2*pi/Nfft);
            % totalCFO(m, l, i) = estCFO*Fs/(2*pi);

            % Correct fine CFO
            totalCFO(m,l,i) = (EstLTF.cfo + EstSTF.cfo) * Fs / (2*pi);
            rxFineCorrected = rxCoarseCorrected .* exp(-1j*EstLTF.cfo*n);
    
    
            %% Equalization and Demodulation
                           
            BERHalfSync(m,l,i) = equalizeAndDemod(rxCoarseCorrected, ...
                EstSTF.index, activeBins, payloadScale, txBits, Mary, ...
                preamLong);

            BERFullSync(m,l,i) = equalizeAndDemod(rxFineCorrected, ...
                EstSTF.index, activeBins, payloadScale, txBits, Mary, ...
                preamLong);
        end
    end
end


%% CFO Estimation Performance
numCFO = length(CFO_norm_values);

% Residual CFO for each CFO, SNR and iteration
trueCFO_Hz_3D = reshape(trueCFO_Hz, [], 1, 1);

subcarrierSpacing = Fs/Nfft;
residualCoarse = (trueCFO_Hz_3D - coarseCFO) ./ subcarrierSpacing;
residualTotal = (trueCFO_Hz_3D - totalCFO) ./ subcarrierSpacing;

% RMSE across Monte Carlo iterations
residualCoarse(~validSynch) = NaN;
residualTotal(~validSynch) = NaN;

rmseCoarseOnly = sqrt(mean(residualCoarse.^2, 3, 'omitnan'));
rmseCoarseFine = sqrt(mean(residualTotal.^2, 3, 'omitnan'));

% Improvement in RMSE
improvementCFO_dB = 20*log10( ...
    rmseCoarseOnly ./ max(rmseCoarseFine, eps));

% Plot 1: Residual CFO RMSE vs SNR
figure('Name', 'Residual CFO RMSE');

semilogy(SNR_dB, rmseCoarseOnly(CFO_norm_values==1.2,:), '-o', ...
    'LineWidth', 1.4, 'MarkerSize', 5, ...
    'DisplayName', 'Coarse Only');
hold on;
semilogy(SNR_dB, rmseCoarseFine(CFO_norm_values==1.2,:), '-o', ...
    'LineWidth', 1.4, 'MarkerSize', 5, ...
    'DisplayName', 'Coarse + Fine');
grid on; box on;
xlabel('SNR (dB)');
ylabel('Normalized Residual CFO RMSE (\epsilon_r)');
legend('Location', 'best');
title('Residual CFO Estimation Performance With Applied CFO \epsilon=1.2');


% Plot 2: Improvement due to fine CFO estimation
figure('Name', 'Fine CFO Improvement');

imagesc(SNR_dB, CFO_norm_values, improvementCFO_dB);
axis xy;
colorbar;
hold on;

xlabel('SNR (dB)');
ylabel('Normalized CFO (subcarrier spacings)');
title('RMSE Improvement from Fine CFO Estimation (dB)');

set(gca, 'YTick', CFO_norm_values);
grid on;
box on;

% Plot 3: RMSE improvement vs SNR
figure('Name', 'RMSE Improvement');
hold on;

for m = 1:numCFO
    plot(SNR_dB, improvementCFO_dB(m,:), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 5, ...
        'DisplayName', sprintf('CFO = %.1f', CFO_norm_values(m)));
end

grid on;
box on;
xlabel('SNR (dB)');
ylabel('RMSE Improvement (dB)');
title('Benefit of Fine CFO Estimation');
legend('Location', 'best');

% Plot 4: BER vs SNR
figure('Name', 'BER vs SNR');

BERHalfSync(~validSynch) = NaN;
BERFullSync(~validSynch) = NaN;

meanBERHalf = mean(BERHalfSync, 3, 'omitnan');
meanBERFull = mean(BERFullSync, 3, 'omitnan');

hold on;
semilogy(SNR_dB, meanBERHalf(CFO_norm_values==1.2,:), '-o', ...
    'LineWidth', 1.4, 'MarkerSize', 5, ...
    'DisplayName', 'Coarse Only');
semilogy(SNR_dB, meanBERFull(CFO_norm_values==1.2,:), '-o', ...
    'LineWidth', 1.4, 'MarkerSize', 5, ...
    'DisplayName', 'Coarse + Fine');

grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER With Applied CFO \epsilon=1.2');
legend;

% Plot 5: BER Gain
figure;

gainBER = 10*log10(meanBERHalf ./ meanBERFull);
gainBER_plot = gainBER;

% Ganancia muy alta cuando BER Coarse + Fine es cero
gainBER_plot(isinf(gainBER_plot) & gainBER_plot > 0) = 50;

% Sin mejora cuantificable cuando ambas BER son cero
gainBER_plot(isnan(gainBER_plot)) = 0;

imagesc(SNR_dB, CFO_norm_values, gainBER);
axis xy;
colorbar;
hold on;

xlabel('SNR (dB)');
ylabel('BER Gain');
title('BER Improvement from Fine CFO Estimation (dB)');

set(gca, 'YTick', CFO_norm_values);
grid on;
box on;

% Plot 6
figure;
hold on;

for k = 1:length(CFO_norm_values)

    plot(SNR_dB, gainBER_plot(k,:), '-o', ...
        'LineWidth', 1.4, ...
        'MarkerSize', 5, ...
        'DisplayName', ...
        sprintf('\\epsilon = %.1f', CFO_norm_values(k)));

end

grid on;
xlabel('SNR (dB)');
ylabel('BER Gain');
title('BER Improvement vs SNR');
legend('Location', 'best');

%% Display results
% fprintf('--- STF synchronization test ---\n');
% fprintf('Actual LTF start index (MATLAB): %d\n', initialDelay+1);
% fprintf('Detected index:                  %d\n', EstSTF.index);
% fprintf('Synchronization successful:     %d\n', EstSTF.syncSucceeded);

%% Plot received signal and detection index
% figure('Name', 'STO Estimation');
% plot(abs(rxFrame));
% hold on;
% 
% xline(initialDelay+1, '--g', 'STF start');
% 
% if EstSTF.syncSucceeded
%     xline(EstSTF.index, '--r', 'STF detected index');
% end
% 
% xlabel('Sample index');
% ylabel('Magnitude');
% title('802.11a frame: Synchronization');
% grid on;
% legend('Received signal', 'STF start', 'STF Detected index');
% 

function BER = equalizeAndDemod(rxFrame, ltfStart, activeBins, ...
    payloadScale, txBits, Mary, preamLong)
GI = 16;
GI2 = 32;
Nfft = 64;
payloadStart = ltfStart + 160;
%% Channel estimation using LTF    

% Extract the two LTF symbols
ltf1 = rxFrame(ltfStart+GI2 : ltfStart+GI2+Nfft-1);
ltf2 = rxFrame(ltfStart+GI2+Nfft : ltfStart+GI2+2*Nfft-1);

Yltf1 = fft(ltf1, Nfft);
Yltf2 = fft(ltf2, Nfft);

% Known transmitted LTF in frequency domain
Xltf = fft(preamLong(GI2+1:GI2+Nfft), Nfft);

% Estimate channel and average both LTF symbols
H_hat = 0.5 * (Yltf1 ./ Xltf + Yltf2 ./ Xltf);

%% Equalization                       
% Remove CP
rxPayload = rxFrame(payloadStart+GI : payloadStart+GI+Nfft-1);

Y = fft(rxPayload, Nfft);

% Zero-forcing equalization
equalizedSymbols = Y(activeBins) ./ H_hat(activeBins);     
equalizedSymbols = equalizedSymbols / payloadScale;

% Constellation
%scatter(real(equalizedSymbols), imag(equalizedSymbols));

%% BER calculation
detectedBits = qamdemod(equalizedSymbols(:), Mary, ...
    'OutputType', 'bit', 'UnitAveragePower', true);

BER = mean(detectedBits ~= txBits);
end