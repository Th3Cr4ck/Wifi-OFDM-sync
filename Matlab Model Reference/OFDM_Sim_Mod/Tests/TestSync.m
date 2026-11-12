addpath('..');
clc;
clear;
close all;

%% Parameters
Fs = 10e6;
Fc = 2.4e9;
Mary = 16;   % QAM order
M = 20;     % OFDM Symbols by frame
Nfft = 64;
GI = 16;
GI2 = 32;
initialDelay = 200;  % Initial zeros
SNR_dB = 5:5:50;
it = 3000;
delayFactor = 4;

syncFailedSTF = false(length(SNR_dB), it);
syncFailedLTF = false(length(SNR_dB), it);
BER = zeros(length(SNR_dB), it);

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

for SNR = SNR_dB
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
        activeSignal = rxFrame(initialDelay+1:end-initialDelay/2);
        signalPower = mean(abs(activeSignal).^2);
        noisePower = signalPower / 10^(SNR/10);
        
        noise = sqrt(noisePower/2) * ...
            (randn(1,length(rxFrame)) + 1j*randn(1, length(rxFrame)));
        
        rxFrame = rxFrame + noise;
        
        %% Synchronization
        EstSTF = WifiSyncSTF(rxFrame, delayFactor);
    
        if ~EstSTF.syncSucceeded
            syncFailedSTF(SNR == SNR_dB, i) = true;
            continue;
        end
        
        % Sync index must be a valid point for the fft
        validIndexRange = trueLTFStart-GI+1:trueLTFStart;
        if ~ismember(EstSTF.index, validIndexRange)
            EstSTF.syncSucceeded = false;
            syncFailedSTF(SNR == SNR_dB, i) = true;
            continue;
        end

        % Correct CFO
        n = 0:length(rxFrame)-1;
        rxFrame = rxFrame .* exp(-1j*EstSTF.cfo*n);

        truePayloadStart = trueLTFStart + length(preamLong);
        payloadStart = EstSTF.index+160;

   
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
 
        BER(SNR == SNR_dB, i) = mean(detectedBits ~= transmittedBits);

    end
end

%% Results 
figure;
ok = ~(syncFailedSTF | syncFailedLTF);
BERmean = sum(BER .* ok, 2) ./ max(sum(ok, 2), 1);
nonzeroElements = length(find(BERmean));
semilogy(SNR_dB(1:nonzeroElements), BERmean(1:nonzeroElements));
xlabel('SNR (dB)');
ylabel('BER');
grid on;

figure;
failRatio = sum(syncFailedLTF | syncFailedSTF, 2) ./it;
plot(SNR_dB, failRatio);
xlabel('SNR (dB)');
ylabel('Sync Failed Ratio');
title('Synchronization Failure Ratio');
grid on;

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