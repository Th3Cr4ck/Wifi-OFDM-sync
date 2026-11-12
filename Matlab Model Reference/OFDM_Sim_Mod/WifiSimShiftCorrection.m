%% TODO: 
% Requiere estimar el canal si se usa TDL para comprobar la correctitud del
% metodo por desplazamiento. En canal AGWN no hace falta.
clear;
close all;
clc;

%% Parameters
Fs = 10e6;
Nfft = 64;
Ng = 16;
Fc = 2.4e9;
SNR_dB = 10:5:50;
iter = 3000;
initialDelay = 200;  % Initial zeros

CFO_norm_values = [0.1 0.2 0.5 0.8 1.2 1.6 1.9]; % CFO normzalized by subcarrier spacing
trueCFO_Hz = CFO_norm_values*Fs/Nfft;
trueCFO = 2*pi*trueCFO_Hz/Fs;  % rad/sample

ACTIVE_TDL_CHANNEL = false;

coarseCFO = zeros(length(CFO_norm_values),length(SNR_dB), iter);
fineCFO = zeros(length(CFO_norm_values), length(SNR_dB), iter);
validSynch = false(length(CFO_norm_values), length(SNR_dB), iter);

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
            %% Generate one OFDM QPSK payload symbol
            % 52 active subcarriers, excluding DC
            activeSubcarriers = [-26:-1 1:26];
            
            payloadSymbols = (2*randi([0 1],1,52)-1) + ...
                1j*(2*randi([0 1],1,52)-1);
            payloadSymbols = payloadSymbols/sqrt(2);
            
            X = zeros(1,Nfft);
            X(mod(activeSubcarriers,Nfft)+1) = payloadSymbols;
            
            payloadUseful = ifft(X);
            payloadUseful = payloadUseful / sqrt(mean(abs(payloadUseful).^2));
            payload = [payloadUseful(end-Ng+1:end), payloadUseful];
            
            %% Assemble transmitted frame
            txFrame = [zeros(1,initialDelay), ...
                       preamShort, ...
                       preamLong, ...
                       payload, ...
                       zeros(1, initialDelay/2)];
   
            %% TDL Channel
            if ACTIVE_TDL_CHANNEL                          
                rxFrame = tdl(transpose(txFrame));
                rxFrame = transpose(rxFrame);
            else
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
            EstSTF = WifiSyncSTF(rxFrame);
            
            % Estimate Coarse CFO
            if ~EstSTF.syncSucceeded
                warning('STF Synchronization not succeeded');
                continue;
            end
            
            coarseCFO(m, l, i) = EstSTF.cfo*Fs/(2*pi);                        
            
            %% Run LTF synchronization          
            EstLTF = WifiSyncLTF(transpose(rxFrame(EstSTF.index:end)), 64);

            % Estimate fine CFO
            if ~EstLTF.syncSucceeded
                warning('LTF Synchronization not succeeded');
                continue;                             
            end
    
            validSynch(m, l, i) = true;
            fineCFO(m, l, i) = EstLTF.cfo*Fs/(2*pi);

            intCFO_norm = round((EstSTF.cfo - EstLTF.cfo)/(2*pi/Nfft));        

            % Correct Fine CFO
            rxFrameFineCorrected = rxFrame .* exp(-1j*EstLTF.cfo*n);
            
            % Correct Total CFO Reference signal
            totalCFORef = EstLTF.cfo + intCFO_norm*(2*pi/Nfft);
            rxFrameReference = rxFrame .* exp(-1j*totalCFORef*n);

            %% Correct Coarse CFO
            
            % Extract symbol and remove CP
            CP = 16;
            symbolStart = EstLTF.index+EstSTF.index-1;
            y = rxFrameFineCorrected(symbolStart+CP : symbolStart+CP+Nfft-1);
            yRef = rxFrameReference(symbolStart+CP : symbolStart+CP+Nfft-1);
            
            % FFT
            Y = fft(y, Nfft);
            YRef = fft(yRef, Nfft);

            YCorrected = circshift(Y, -intCFO_norm);

            n0 = symbolStart + CP - 1;  % Índice de la primera muestra, empezando en 0            
            phaseFactor = exp(-1j*2*pi*intCFO_norm*n0/Nfft);
            % Incluir la fase global para comparar con la referencia temporal
            YCorrected = phaseFactor * YCorrected;


            relativeError = norm(YRef - YCorrected)^2 / norm(YRef)^2;
        end
    end
end


%% CFO Estimation Performance

numCFO = length(CFO_norm_values);

% Residual CFO for each CFO, SNR and iteration
trueCFO_Hz_3D = reshape(trueCFO_Hz, [], 1, 1);

subcarrierSpacing = Fs/Nfft;
residualCoarse = (trueCFO_Hz_3D - coarseCFO) ./ subcarrierSpacing;
residualTotal = (trueCFO_Hz_3D - fineCFO) ./ subcarrierSpacing;

% RMSE across Monte Carlo iterations
residualCoarse(~validSynch) = NaN;
residualTotal(~validSynch) = NaN;

rmseCoarseOnly = sqrt(mean(residualCoarse.^2, 3, 'omitnan'));
rmseCoarseFine = sqrt(mean(residualTotal.^2, 3, 'omitnan'));

% Improvement in RMSE
improvement_dB = 20*log10( ...
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

imagesc(SNR_dB, CFO_norm_values, improvement_dB);
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
    plot(SNR_dB, improvement_dB(m,:), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 5, ...
        'DisplayName', sprintf('CFO = %.1f', CFO_norm_values(m)));
end

grid on;
box on;
xlabel('SNR (dB)');
ylabel('RMSE Improvement (dB)');
title('Benefit of Fine CFO Estimation');
legend('Location', 'best');


%% Display results
fprintf('--- STF synchronization test ---\n');
fprintf('Actual STF start index (MATLAB): %d\n', initialDelay+1);
fprintf('Detected index:                  %d\n', EstSTF.index);
fprintf('Synchronization successful:     %d\n', EstSTF.syncSucceeded);

%% Plot received signal and detection index
figure('Name', 'STO Estimation');
plot(abs(rxFrame));
hold on;

xline(initialDelay+1, '--g', 'STF start');

if EstSTF.syncSucceeded
    xline(EstSTF.index, '--r', 'STF detected index');
end

xlabel('Sample index');
ylabel('Magnitude');
title('802.11a frame: Synchronization');
grid on;
legend('Received signal', 'STF start', 'STF Detected index');

%% Generate Preamble Functions
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