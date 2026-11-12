addpath('..');

clear; close all; clc;
rng(1);

%% Parametros
N    = 64;
Ng   = 16;
Tu   = 4e-6;        % Duracion del simbolo OFDM
fc   = 2.4e9;
numSymbols = 20;
mary = 4;

% Offsets verdaderos (oscilador comun: CFO y SFO comparten epsilon)
epsTrue   = 20e-6;              % epsilon total
eps0Hat   = 5e-6;               % estimacion inicial de LT (eps_0)
epsR      = epsTrue - eps0Hat;  % CFO residual tras LT: eps_r = eps - eps_0
initialOffset = eps0Hat;

% Pilotos 802.11a
pilotIndex = [-21 -7 7 21];
pilotTx    = [1 1 1 -1].';
pilotPos   = pilotIndex + N/2 + 1; % Para matlab indexing

%% Trama OFDM
activeIndex = [-26:-1 1:26];
activePos   = activeIndex + N/2 + 1;
dataIndex   = [-26:-22 -20:-8 -6:-1 1:6 8:20 22:26];
dataPos     = dataIndex + N/2 + 1;

pilotSel      = ismember(activeIndex, pilotIndex).';   % columna, sobre activas
dataActivePos = ismember(activeIndex, dataIndex).';

txFreq = zeros(N, numSymbols);
dataBits    = randi([0 1], length(dataIndex)*log2(mary)*numSymbols, 1);
dataSymbols = qammod(dataBits, mary, 'InputType', 'bit', 'UnitAveragePower', true);
txFreq(dataPos, :)  = reshape(dataSymbols, [], numSymbols);
txFreq(pilotPos, :) = repmat(pilotTx, 1, numSymbols);

%% Modelo del paper
l = 1:numSymbols;
k = activeIndex;
phase = 2*pi*l .* (epsR*Tu*fc - (N+Ng)/N * epsTrue .* k.');

rxFreq = txFreq;
rxFreq(activePos, :) = rxFreq(activePos, :) .* exp(1j*phase);

%% Estimacion fina
betaEst    = zeros(numSymbols, 1);   % beta_l (fase acumulada)
dphiEst    = zeros(numSymbols, 1);   % angle(W)
fineCfoEst = zeros(numSymbols, 1);   % eps_r estimado
offsetEst  = zeros(numSymbols, 1);   % eps estimado (total)
rxCorrected = zeros(length(activeIndex), numSymbols);

prevOffsetEstimate = initialOffset;
prevPilotValues    = zeros(4, 1);
W = 0;

for l = 1:numSymbols
    rxSymbol = rxFreq(activePos, l);

    % 1. Corregir SFO (incluye pilotos), con eps_{l-1}
    sfoPhase = computeSfoPhase(l, N, Ng, prevOffsetEstimate, activeIndex);
    rxSymbol = rxSymbol .* exp(1j*sfoPhase(:));

    % Pilotos con polaridad conocida removida
    pilotValues = rxSymbol(pilotSel) .* pilotTx;

    % 2. Corregir CFO residual (solo datos; pilotos se dejan para el paso 3)
    channelGainAtPilots = ones(4,1);
    beta = estimateResidualCfoPhase(pilotValues, channelGainAtPilots);
    betaEst(l) = beta;
    rxSymbol(dataActivePos) = rxSymbol(dataActivePos) * exp(-1j*beta);

    % 3. Estimacion fina de epsilon para el siguiente simbolo
    if l > 1
        [fineCfo, estOffset, W] = estimateFineOffset( ...
            pilotValues, prevPilotValues, Tu, fc, initialOffset, W);
        dphiEst(l)    = angle(W);
        fineCfoEst(l) = fineCfo;
        offsetEst(l)  = estOffset;
        prevOffsetEstimate = estOffset;   % solo para el SFO del siguiente simbolo
    else
        offsetEst(l) = prevOffsetEstimate;
    end

    prevPilotValues = pilotValues;
    rxCorrected(:, l) = rxSymbol;
end

%% Referencias y metricas
wrap = @(a, b) angle(exp(1j*(a - b)));

l = 1:numSymbols;
betaTrue = angle(exp(1j*2*pi*l.' * epsR*Tu*fc));   % acumulada, envuelta
betaErr  = wrap(betaEst, betaTrue);

dphiTrue = 2*pi*epsR*Tu*fc;   % incremento fase por simbolo
idx = 2:numSymbols;           % W_l existe desde l = 2

% RMSE por simbolo sobre portadoras de datos
err = rxCorrected(dataActivePos, :) - txFreq(activePos(dataActivePos), :);
rmse = sqrt(mean(abs(err).^2, 1)).';

fprintf("\n--- beta_l (fase acumulada) ---\n");
fprintf("Error absoluto medio (envuelto): %.3e rad\n", mean(abs(betaErr)));

fprintf("\n--- Incremento por simbolo angle(W) ---\n");
fprintf("Esperado : %.6f rad\n", dphiTrue);
fprintf("Estimado : %.6f rad (l = %d)\n", dphiEst(end), numSymbols);

fprintf("\n--- CFO residual (eps_r) ---\n");
fprintf("Real     : %.8e\n", epsR);
fprintf("Estimado : %.8e\n", fineCfoEst(end));
fprintf("Error    : %.3e\n", fineCfoEst(end) - epsR);

fprintf("\n--- Offset total (eps) ---\n");
fprintf("Real     : %.8e\n", epsTrue);
fprintf("Estimado : %.8e\n", offsetEst(end));
fprintf("Error    : %.3e\n", offsetEst(end) - epsTrue);

fprintf("\n--- RMSE (datos) ---\n");
fprintf("l=1: %.3e | l=2: %.3e | l=3: %.3e | l=%d: %.3e\n", ...
    rmse(1), rmse(2), rmse(3), numSymbols, rmse(end));

%% Graficas
figure;
subplot(2,2,1);
plot(1:numSymbols, betaEst, 'o-'); hold on;
plot(1:numSymbols, betaTrue, 'x--'); grid on;
xlabel('Simbolo OFDM'); ylabel('\beta_l [rad]');
legend('Estimado', 'Real (envuelto)'); title('Fase acumulada de CFO residual');

subplot(2,2,2);
plot(idx, fineCfoEst(idx), 'o-'); hold on;
yline(epsR, '--'); grid on;
xlabel('Simbolo OFDM'); ylabel('\epsilon_r');
legend('Estimado', 'Real'); title('CFO residual (fino)');

subplot(2,2,3);
plot(1:numSymbols, offsetEst, 'o-'); hold on;
yline(epsTrue, '--'); grid on;
xlabel('Simbolo OFDM'); ylabel('\epsilon');
legend('Estimado', 'Real'); title('Offset total (se usa para el SFO)');

subplot(2,2,4);
semilogy(1:numSymbols, max(rmse, eps), 'o-'); grid on;
xlabel('Simbolo OFDM'); ylabel('RMSE');
title('RMSE de datos tras corregir SFO y CFO');



figure;
scatterplot(rxCorrected);