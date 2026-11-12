

ConstSize=4;
NFFT=64; %FFT size
CP=16;
Fs=20e6; %Tasa se muestreo
Ts=1/Fs;
BandasG=NFFT/4-4; %En total
FrameSize=40; %datos



%% Carrier indexes
Carriers=1:NFFT;
%DataCarrierIn=Carriers;
GIndex=[1:BandasG/2 NFFT/2+1 NFFT-BandasG/2+2:NFFT];
UsedIndex=Carriers;
UsedIndex(GIndex)=[];
PilotIndex=UsedIndex(2:floor(numel(UsedIndex)/10):end);
DataCarrierIn=Carriers;
DataCarrierIn([GIndex PilotIndex])=[];

[~,LogicIndex]=sort([GIndex DataCarrierIn PilotIndex]);

%%
Tsymbol=Ts*(NFFT+CP);
FrameTime=Tsymbol*(FrameSize+2);
DataCarriers=numel(DataCarrierIn);
DataCarriersFrame=DataCarriers*FrameSize;
BitpSymbol=DataCarriers*log2(ConstSize);
BitRate=BitpSymbol/Tsymbol;
%% BEM
pilots=randsrc(numel(PilotIndex),1,[-1 1 1i 1i])*.71;
L=10;
B1=eye(L);
B1=[B1;zeros(NFFT-L,L)];

K=L;
dps_seq = dpss(K,K*(NFFT-BandasG+1)/NFFT/2,K);
% plot(dps_seq)
% size(dps_seq)
dps_seq=[dps_seq;zeros(NFFT-size(dps_seq,1),size(dps_seq,2))];
% size(dps_seq)
dps_seq=[dps_seq;zeros(NFFT-K,size(dps_seq,2))];
% size(dps_seq)


PHI=fftshift(fft(dps_seq,[],1),1);
PHIpilots=PHI(PilotIndex,:);
A=eye(numel(pilots))*PHIpilots;
Ainv=(A'*A)\A';
ChannEstData=PHI(DataCarrierIn,:)*Ainv;


%%Sincronia de trama
PreambuloSym=randsrc(NFFT/2,1,[-1 1]);
PreambuloFrec=upsample(PreambuloSym,2);
PreambuloFrec([1:12 end-11:end])=0; %guardas
PreambuloFrec=zeros(64,1);
PreambuloFrec(9:4:57)=4*(1+1i)*[1,-1,1,-1,-1,1,0,-1,-1,1,1,1,1];
SNRdb=20;


%Variables para sincronia con diferencias
Threeshold_umbral = 0.10;

%% SIMULACIONES
% SNRdb=1;
% out = sim('OFDMmodel_synch',1);
% Data_1dB = out.TxData(:,1,1:10);
% pause(2);
% 
% SNRdb=5;
% out = sim('OFDMmodel_synch',1);
% Data_5dB = out.TxData(:,1,1:10);
% pause(2);
% 
% SNRdb=10;
% out = sim('OFDMmodel_synch',1);
% Data_10dB = out.TxData(:,1,1:10);
% pause(2);
% 
% SNRdb=15;
% out = sim('OFDMmodel_synch',1);
% Data_15dB = out.TxData(:,1,1:10);
% pause(2);
% 
% SNRdb=20;
% out = sim('OFDMmodel_synch',1);
% Data_20dB = out.TxData(:,1,1:10);
% pause(2);
% 
% save('Data_Rx.mat','Data_1dB','Data_5dB','Data_10dB','Data_15dB','Data_20dB');




% SNR=0:5:20;
% BER=zeros(numel(SNR));
% K=5;
% for SNRdb = SNR
%     tempBER=0;
%     NoiserSigma2= 10^(-SNRdb/10);
%     for i=1:K
%      simput=sim('OFDMmodel');
%      tempBER=simput.yout{1}.Values.Data+tempBER;
%     end
%     BER(SNRdb == SNR)=tempBER/K;
%     SNRdb
% end
% figure
% semilogy(SNR,BER);
% 
% 
% 

%% Post morten sinch
% L =NFFT/2;
% P0 = 0;
% for d=1:L
%      P0 = P0 + OFDM_frame1(d)' * OFDM_frame1(d+L);
% end
% P = zeros(1,3360-2*L);
% 

% P(1) = P0;
% for d=1:3360-2*L
%      P(d+1) = P(d) + OFDM_frame1(d+L)' * OFDM_frame1(d + 2*L) - OFDM_frame1(d)' * OFDM_frame1(d+L);
% end
% 
% R0 = 0;
% for d=1:L
%      R0 = R0 + abs(OFDM_frame1(d+L))^2;
% end
% 
% R = zeros(1,3360-2*L);
% 
% R(1) = R0;
% for d=1:3360-2*L
%      R(d+1) = R(d) + abs(OFDM_frame1(d+2*L))^2 - abs(OFDM_frame1(d+L))^2;
% end
% 
% m = (abs(P).^2)./(R).^2;



