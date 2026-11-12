clear;
close all;
clc;

%% Parameters
Fs = 10e6;
Fc = 2.4e9;
Mary = 4;   % QAM order
M = 20;     % OFDM Symbols by frame
Nfft = 64;
GI = 16;
GI2 = 32;
SNR_dB = 5:5:50;
initialDelay = 200;  % Initial zeros
iter = 4000;
delayFactor = 4;

CFO_norm = 1.2; % CFO normzalized by subcarrier spacing
trueCFO_Hz = CFO_norm*Fs/Nfft;
trueCFO = 2*pi*trueCFO_Hz/Fs;  % rad/sample

ACTIVE_TDL_CHANNEL = true;

coarseCFO = zeros(length(SNR_dB), iter);
totalCFO = zeros(length(SNR_dB), iter);
syncFailedSTF = false(length(SNR_dB), iter);
syncFailedLTF = false(length(SNR_dB), iter);
BER = zeros(length(SNR_dB), iter);

%% Generate preamble
preamShort = generatePreambleSTF();
preamLong = generatePreambleLTF();

% Normalize power
preamPower = 1;
scaleFactorShort = sqrt(preamPower) / sqrt(mean(abs(preamShort).^2));
scaleFactorLong = sqrt(preamPower) / sqrt(mean(abs(preamLong).^2));

preamShort = preamShort * scaleFactorShort;
preamLong = preamLong * scaleFactorLong;

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


for l = 1:length(SNR_dB)
    lSNR = SNR_dB(l);
    for i = 1:iter
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
     
        %% Apply CFO
        n = 0:length(rxFrame)-1;
        rxFrame = rxFrame .* exp(1j*trueCFO*n);

        %% Add noise
        activeSignal = rxFrame(initialDelay+1:end-initialDelay/2);
        signalPower = mean(abs(activeSignal).^2);
        noisePower = signalPower / 10^(lSNR/10);
        
        noise = sqrt(noisePower/2) * ...
            (randn(1,length(rxFrame)) + 1j*randn(1, length(rxFrame)));
        
        rxFrame = rxFrame + noise;
                      
        %% STF synchronization
        EstSTF = WifiSyncSTF(rxFrame, delayFactor);
                   
        if ~EstSTF.syncSucceeded
            syncFailedSTF(l, i) = true;
            continue;
        end
       
        % Sync index must be a valid point for the fft
        validIndexRange = trueLTFStart-GI+1:trueLTFStart;
        if ~ismember(EstSTF.index, validIndexRange)
            EstSTF.syncSucceeded = false;
            syncFailedSTF(l, i) = true;
            continue;
        end
        
        coarseCFO(l, i) = EstSTF.cfo*Fs/(2*pi); 

        % Correct coarse CFO
        %rxFrameCoarseCorrected = rxFrame .* exp(-1j*EstSTF.cfo*n);            
        
        %% LTF synchronization          
        EstLTF = WifiSyncLTF(transpose(rxFrame(trueLTFStart:end)));
        
        if ~EstLTF.syncSucceeded
            syncFailedLTF(l, i) = true;
            continue;                                             
        end

        % Correct coarse + fine CFO
        intCFO_norm = round((EstSTF.cfo - EstLTF.cfo)/(2*pi/Nfft));
        estCFO = EstLTF.cfo + intCFO_norm * (2*pi/Nfft);

        totalCFO(l, i) = estCFO*Fs/(2*pi);        
        
        rxFrame = rxFrame .* exp(-1j*estCFO*n);

        %% Channel estimation using LTF

        ltfStart = EstSTF.index;  % Start of LTF, including GI2
               
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

        truePayloadStart = trueLTFStart + length(preamLong);
        payloadStart = EstSTF.index+160;

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
 
        BER(l, i) = mean(detectedBits ~= transmittedBits);
    end
end


%% Plot received signal and detection index
figure;
plot(abs(rxFrame));
hold on;

xline(trueSTFStart, '--g', 'STF start');
xline(trueLTFStart, '--', 'LTF start')

if EstSTF.syncSucceeded
    xline(EstSTF.index, '--r', 'LTF start detected');
end

xlabel('Sample index');
ylabel('Magnitude');
title('802.11a frame: Synchronization');
grid on;
legend('Received signal', 'STF start', 'LTF start', 'LTF Detected index');

figure;
ok = ~(syncFailedSTF | syncFailedLTF);
BERmean = sum(BER .* ok, 2) ./ max(sum(ok, 2), 1);
semilogy(SNR_dB, BERmean);
grid on;
title('BER vs SNR');
xlabel('SNR (dB)');
ylabel('BER');

%% CFO Estimation Performance

numCFO = length(CFO_norm);

% Residual CFO for each CFO, SNR and iteration
trueCFO_Hz_3D = reshape(trueCFO_Hz, [], 1, 1);

subcarrierSpacing = Fs/Nfft;
residualCoarse = (trueCFO_Hz_3D - coarseCFO) ./ subcarrierSpacing;
residualTotal = (trueCFO_Hz_3D - totalCFO) ./ subcarrierSpacing;

% RMSE across Monte Carlo iterations
residualCoarse(~ok) = NaN;
residualTotal(~ok) = NaN;

rmseCoarseOnly = sqrt(mean(residualCoarse.^2, 2, 'omitnan'));
rmseCoarseFine = sqrt(mean(residualTotal.^2, 2, 'omitnan'));

% Improvement in RMSE
improvement_dB = 20*log10( ...
    rmseCoarseOnly ./ max(rmseCoarseFine, eps));

% Plot 1: Residual CFO RMSE vs SNR
figure('Name', 'Residual CFO RMSE');

semilogy(SNR_dB, rmseCoarseOnly, '-o', ...
    'LineWidth', 1.4, 'MarkerSize', 5, ...
    'DisplayName', 'Coarse Only');
hold on;
semilogy(SNR_dB, rmseCoarseFine, '-o', ...
    'LineWidth', 1.4, 'MarkerSize', 5, ...
    'DisplayName', 'Coarse + Fine');
grid on; box on;
xlabel('SNR (dB)');
ylabel('Normalized Residual CFO RMSE (\epsilon_r)');
legend;
t = sprintf('Residual CFO Estimation Performance With Applied CFO \e=%f', CFO_norm);
title(t);

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