%% Config Pameters
%addpath('..');
clear;
close all;
it=10;
SNRvals=5:5:50;
OFDMConfig.Fs = 10e6; %Sampling Rate
OFDMConfig.Ts = 1/OFDMConfig.Fs;
OFDMConfig.N=1024;    %FFT Size 2^power
OFDMConfig.EffBW=0.85; % Efective ocupied bandwidth
OFDMConfig.MRU=8;    %Resource Elements Time, total OFDM symbols is M=MRU*RUStr.size(2)
OFDMConfig.CP=32;   %Cyclic prefix
OFDMConfig.Mary=16;
OFDMConfig.PilotPw=2;
OFDMConfig.PrePW=2;
OFDMConfig.VTmode=0;
OFDMConfig.DFTs=0;
OFDMConfig.DFTsFrame=1;

LMMSETv=1;
OFDMSymbolLen = OFDMConfig.N+OFDMConfig.CP;

OFDMConfig.WIFI_PREAM = false;
fails = 0;

%% OFDM RE Map Definition
%Resource Element Map Definition 1-Data, 2-Pilot, 0-Empty
% ------- LTE --------- %
% OFDMConfig.RUStrmap=uint8(...
%     [1, 1, 2, 1, 1, 1, 2;...   %Subcarrier-Rows   Symbol-Column 
%      2, 1, 1, 1, 1, 1, 1;...
%      1, 1, 1, 1, 2, 1, 1;...
%      1, 1, 1, 1, 1, 1, 1;...
%      1, 1, 1, 1, 1, 1, 1;...
%      1, 1, 1, 2, 1, 1, 1;...
%      1, 1, 1, 1, 1, 1, 1;...
%      1, 1, 1, 1, 1, 1, 2;...
%      1, 1, 1, 1, 2, 1, 1;...
%      1, 1, 1, 1, 1, 1, 1;...
%      2, 1, 1, 1, 1, 1, 1;...
%      1, 1, 1, 2, 1, 1, 1]);

OFDMConfig.RUStrmap=uint8(...
    [2, 1, 2, 1, 2, 1, 2;...   %Subcarrier-Rows   Symbol-Column 
     1, 1, 1, 1, 1, 1, 1;...
     1, 1, 1, 1, 1, 1, 1;...
     2, 1, 2, 1, 2, 1, 2;...
     1, 1, 1, 1, 1, 1, 1;...
     1, 1, 1, 1, 1, 1, 1;...
     2, 1, 2, 1, 2, 1, 2;...
     1, 1, 1, 1, 1, 1, 1;...
     1, 1, 1, 1, 1, 1, 1;...
     2, 1, 2, 1, 2, 1, 2;...
     1, 1, 1, 1, 1, 1, 1;...
     1, 1, 1, 1, 1, 1, 1]);

% OFDMConfig.RUStrmap=uint8(...
%     [2, 1, 2, 1, 2, 1, 2;...   %Subcarrier-Rows   Symbol-Column 
%      1, 1, 1, 1, 1, 1, 1;...
%      0, 1, 0, 1, 0, 1, 0;...
%      2, 1, 2, 1, 2, 1, 2;...
%      0, 1, 0, 1, 0, 1, 0;...
%      1, 1, 1, 1, 1, 1, 1;...
%      0, 1, 0, 1, 0, 1, 0;...
%      2, 1, 2, 1, 2, 1, 2;...
%      0, 1, 0, 1, 0, 1, 0;...
%      1, 1, 1, 1, 1, 1, 1;...
%      0, 1, 0, 1, 0, 1, 0;...
%      2, 1, 2, 1, 2, 1, 2]);

% ------- WIFI --------- %
% OFDMConfig.RUStrmap=uint8( ...
%     [1;1;1;1;1;1;2;...
%      1;1;1;1;1;1;1;...
%      1;1;1;1;1;1;2;...
%      1;1;1;1;1]);  %WiFi
% with VTmode=0

%% Channel Model
v = 100.0;                   % UE velocity in km/h
fc = 3e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
OFDMConfig.fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz

tdl = nrTDLChannel;
tdl.DelayProfile = 'TDL-A';
tdl.DelaySpread = 400e-9;
tdl.MaximumDopplerShift = OFDMConfig.fd;
tdl.SampleRate=OFDMConfig.Fs;
tdl.NumReceiveAntennas=1;


%% OFDM Parameters
[OFDMParams] = OFDMGetParams(OFDMConfig);

%% BEM 2D set up
[A_Pseu,PHI2D_cell] = BEM2D(OFDMParams,OFDMConfig);

%Path2TapMtx=transpose(chan.info.ChannelFilterCoefficients);
Path2TapMtx=tdl.getPathFilters;

%% Figures
if ~exist('f1','var')
    f1=figure;
    fConsRxZF=scatter([],[],'.b');hold on
    %fConsRxLMMSE=scatter([],[],'.r');
    f1.CurrentAxes.YLimMode="manual";
    f1.CurrentAxes.XLimMode="manual";
    f1.CurrentAxes.XLim = [-2 2];
    f1.CurrentAxes.YLim = [-2 2];

    % Appr Delay Doppler
    f5=figure;
    fDDCir = surf(zeros(size(Path2TapMtx,1),OFDMParams.M+1)); 
    shading interp; colormap colorcube
end

%% Simulation
BER=zeros(numel(SNRvals),1);

for SNR=SNRvals
    sprintf('SNR=%0.5g',SNR)
    errors=0;
    
    for k=1:it
        %% Tx
        [TxBits,TxDataSymbols] = OFDMGenPayLoad(OFDMParams);
        TxSignal = OFDMGenTxFrame(OFDMParams,TxDataSymbols);        
        
        %% Channel
        % Simulate frame start with leading zeros
        nZeros = OFDMConfig.N/2;
        RxSignalTx = [zeros(nZeros, 1); TxSignal];
        
        % Simulate channel
        [RxSignal, ChannCoeffs] = tdl(RxSignalTx);
        
        % Add AWGN
        NoisePW = 10^(-SNR/10)*OFDMParams.TxPW;
        awgn = sqrt(NoisePW/2) * ...
            (randn(numel(RxSignal),1) + 1j*randn(numel(RxSignal),1));
        
        RxSignal = RxSignal + awgn;
        
        % Channel impulse response
        ChannCoeffs = transpose(ChannCoeffs);
        % Remove the coefficients corresponding to the
        % artificial samples before the OFDM frame.
        ChannCoeffs = ChannCoeffs(:,nZeros+1:end);
        CIR = Path2TapMtx*ChannCoeffs;
        
        % Approximate Delay-Doppler
        fDDCir(1).ZData = abs(fftshift( ...
            fft(CIR(:,1:(OFDMSymbolLen):end),[],2),2));
        
        CIR = reshape(transpose(CIR), ...
            [OFDMSymbolLen, ...
             OFDMParams.M+1, size(Path2TapMtx,1)]);
        
        CIR = CIR(OFDMParams.CP+1:end,2:end,:);
        CIR = permute(CIR,[1,3,2]);
        CIRFr = squeeze(sum(CIR,1));
        
        %% RX
        
        % Simulate CFO
        L = OFDMParams.N/2;
        % CFO angular: rad/sample
        % Ambiguity-free region: |cfo| < pi/L
        cfo = pi/(L)/2;
        
        n = (0:length(RxSignal)-1).';
        RxSignal = RxSignal .* exp(1j*cfo*n);
        
        %% Synchornization
        % Frame syncronization        
        sync = OfdmSync(RxSignal, L, OFDMConfig.CP);
        %sync = WifiSyncLTF(RxSignal);
        if ~sync.syncSucceeded
            fails = fails + 1;
            continue
        end
    
        % Remove CFO
        RxSignalSynch = RxSignal .* exp(-1j*sync.cfo*n);
    
        % Crop frame to synchronized start index
        nOFDMSymbols = floor((length(RxSignal)-sync.index)/OFDMSymbolLen); % M  
        RxSignalSynch = RxSignalSynch(sync.index:sync.index+(nOFDMSymbols*OFDMSymbolLen)-1);

        RxSignalSynchRect = reshape(RxSignalSynch, ...
        OFDMSymbolLen, nOFDMSymbols);
                
        % Remove CP
        RxFrame = RxSignalSynchRect(OFDMParams.CP+1:end,:); 
        
        % Convert to freq
        RxFrecFrame = 1/sqrt(OFDMParams.N) * ...
            fftshift(fft(RxFrame,[],1),1);
                
       
        %% Channel Estimation        
        G = Channel_Est(OFDMParams,OFDMConfig,A_Pseu,PHI2D_cell,RxFrecFrame);        

        %% Data Recovery
        if(LMMSETv==1)            
            [RxBits, RxSymb] = OfdmLmmseEq(OFDMParams,RxFrecFrame,CIR,NoisePW+5e-3,G);
        else
            [RxBits, RxSymb] = OfdmZfoneTap(OFDMParams,RxFrecFrame,CIR,NoisePW,G);
        end        
        errors=errors+sum(xor(TxBits,RxBits));
        
        %Refresh Graphs
       %  if mod(k,1)==0
       %     RxSymbZF=RxSymbZF(1:1:OFDMParams.HmL);
       %     fConsRxZF.XData=(real(RxSymb));
       %     fConsRxZF.YData=(imag(RxSymb));
       %     drawnow
       %     figure(3); clf; 
       %     scatter(real(RxSymb), imag(RxSymb), 'red', 'DisplayName', 'Rx');           
       %     legend;
       %     titleStr = sprintf('RX Constelation SNR=%d', SNR);
       %     title(titleStr);
       %     drawnow
       % end
    end
        
    BER(SNR==SNRvals)=errors/(OFDMParams.SizeDataBits*it);
    
end

%% BER vs SNR
figure;
plot(SNRvals, BER, 'LineWidth', 2, 'Marker', 'x');
title('BER vs SNR');
xlabel('SNR');
ylabel('BER');

fprintf('Total Synchronization fails: %d/%d\n', fails, it*length(SNRvals));