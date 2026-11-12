clear
close all
Fs=10e6; %Tasa de muestreo
NFFT=256; %Tamaño de FFT
CP=32; %Tamaño de prefijo cíclico
M=16; %Indice de modulación
Ng=50; %Subportadoras de Guarda
FrameL=32;
NullSimb=2;
Fdmax=50;
Lmax=CP;%dispersión maxima en taps, soportada por el estimador

SNRdB=30;
%PDP=transpose(exp(-0.6*(0:12))); %Perfil de Pot de Ret
PilotSpacing=6;
%dW=2*pi*dF/Fs;
MaxDelay=NFFT*10;

chan = comm.RayleighChannel(...
"SampleRate",          Fs,...
"PathDelays",          [0 100e-9 230e-9 490e-9],...
"AveragePathGains",    [0 -2 -8 -12],...
"NormalizePathGains",  true,...c
"MaximumDopplerShift", 50,...
"PathGainsOutputPort", false);

f1=figure(1);
pos=f1.Position;
f1.Position=[0,50 pos(3:4)];
%sf1=subplot(1,2,1);
hline = scatter(nan,nan,'.c'); hold on
hline2 = scatter(nan,nan,'.b');
xlim([-1.5 1.5]);ylim([-1.5 1.5]);grid on ;


f2=figure(2);
pos=f1.Position;
f2.Position=[0,pos(4)+150, pos(3:4)];
f2p1=plot(nan);
hold on
f2p2=plot(nan);
f2p3=plot(nan);
f2s1=scatter(nan,nan,80,'o','LineWidth',2);
f2s2=scatter(nan,nan,80,'o','LineWidth',2);
f2s3=scatter(nan,nan,80,'o','LineWidth',2);
ylim([0 1]);grid on ;

%% Calcúlo de parametros
Ts=1/Fs; %Periodo de muestreo
Nmuestras=NFFT+CP; %Tamaño total del símbolo OFDM
Tsimbolo=Nmuestras*Ts; %Periodo de simbolo OFDM
Indices=1:NFFT;

%Indices de portadora nulos y no Nulos
guardIndices=[1:Ng/2,NFFT/2+1, NFFT-Ng/2+2:NFFT]; % Indices de guardas
nonGuardIndices=Indices; 
nonGuardIndices(guardIndices)=[];%Indices que no son guardas

%Indice de portadoras con pilotos
pilotsIndices=nonGuardIndices(PilotSpacing/2:PilotSpacing:end); %Indices de pilotos
nPilotos=numel(pilotsIndices);

%Calculo de indices de portadora de datos
dataIndices=Indices; 
dataIndices([guardIndices pilotsIndices])=[];%Indices de datos
NDataInd=numel(dataIndices); 

%Valores en las portadoras piloto
%valoresPiloto=sqrt(5)*(2*(randn(nPilotos,1)>0)-1);

pnseq = comm.PNSequence('Polynomial', [7 6 0], ...
'SamplesPerFrame', nPilotos, 'InitialConditions', [0 0 0 0 0 0 1]);
valoresPiloto = sqrt(2)*(2*pnseq()-1);

Ps=(NDataInd+2*nPilotos)/NFFT;
Pn=Ps*10^(-SNRdB/10);
%PDP=PDP/(sum(PDP));

% Preambulo
%preamFrec=2*(randn(NFFT/2,1)>0)-1;

pnseq = comm.PNSequence('Polynomial', [10 7 0], ...
'SamplesPerFrame', NFFT/2, 'InitialConditions', [1 0 0 1 0 1 0 0 0 1]);
preamFrec = 1/sqrt(NFFT/2)*fft(2*pnseq()-1);

preamFrec=upsample(preamFrec,2);
preambTime=3*sqrt(NFFT)*ifft(fftshift(preamFrec));

%BEM
bemtimeN=2*ceil(Tsimbolo*FrameL*Fdmax)+3;
bemDelyN=min(2*ceil(Lmax*(NFFT-Ng+1)/(2*NFFT))+3,Lmax);

[bemDelay,bemweight]=dpss(Lmax,Lmax*(NFFT-Ng+1)/(2*NFFT),bemDelyN);
[bemTime,~]=dpss(FrameL,Tsimbolo*FrameL*Fdmax,bemtimeN);
bemFrec=1/sqrt(NFFT)*fftshift(fft([bemDelay;zeros(NFFT-size(bemDelay,1),bemDelyN)],[],1),1);
% bemFrecPilots=bemFrec(pilotsIndices,:);
% bemFrecData=bemFrec(dataIndices,:);
% bemIntMtx=bemFrecData*((bemFrecPilots'*bemFrecPilots)\bemFrecPilots');

bem2D=zeros(NFFT,FrameL,bemDelyN*bemtimeN);
bem2DdelayTime=zeros(Lmax,FrameL,bemDelyN*bemtimeN);
for ii=1:bemtimeN
     for jj=1:bemDelyN
         bem2D(:,:,jj+(ii-1)*bemDelyN)=bemFrec(:,jj)*transpose(bemTime(:,ii));
         bem2DdelayTime(:,:,jj+(ii-1)*bemDelyN)=bemDelay(:,jj)*transpose(bemTime(:,ii));
     end
end
bem2Dpilots=bem2D(pilotsIndices,:,:);
bem2DData=bem2D(dataIndices,:,:);
bem2Dpilots=reshape(bem2Dpilots,nPilotos*FrameL,bemDelyN*bemtimeN);
bem2DData=reshape(bem2DData,NDataInd*FrameL,bemDelyN*bemtimeN);
bem2DdelayTime=reshape(bem2DdelayTime,Lmax*FrameL,bemDelyN*bemtimeN);

bem2DFull=reshape(bem2D,NFFT*FrameL,bemDelyN*bemtimeN);
            
bem2DIntMtxData=bem2DData*((bem2Dpilots'*bem2Dpilots)\bem2Dpilots');
bem2DIntMtxFull=bem2DFull*((bem2Dpilots'*bem2Dpilots)\bem2Dpilots');
%bem2DIntMtxData=bem2DIntMtxData.*((20*log10(abs(bem2DIntMtxData)))>-50);

bem2DIntMtxDelayTime=bem2DdelayTime*((bem2Dpilots'*bem2Dpilots)\bem2Dpilots');
%agc1 = comm.AGC("AveragingLength",NFFT*5,"AdaptationStepSize",0.01,"MaxPowerGain",1);

f3=figure(3);
pos=f2.Position;
f3.Position=[pos(3)+10,pos(4)+150, pos(3:4)];
sf1=surf(zeros(NFFT,FrameL)); shading interp
delayrnd=randi(MaxDelay+1)-1;
idealSincDelay=delayrnd+chan.info.ChannelFilterDelay;

f4=figure(4);
pos=f3.Position;
f4.Position=[pos(1),50, pos(3:4)];
sf2=imagesc(zeros(NFFT,FrameL)); shading interp

f5=figure(5);
f5p1=semilogy(nan);
f5.CurrentAxes.YLim=[1e-7,1e-1];
ber=zeros(2000,1);
iters=1:2000;
terros=0;
k=0;
peakPw=0;
%% transmisor
for i=1:2000

     if ~mod(i,10)
         chan.reset;
     end


txBits =randi(2,NDataInd*log2(M),FrameL)-1;
txSimbols=qammod(txBits,M,'InputType','bit','UnitAveragePower',true);

%txSimbols=1/sqrt(NDataInd)*1/sqrt(FrameL)*fft2(txSimbols);
txSimbols=1/sqrt(NDataInd)*fft(txSimbols,[],1);
txSimbols=1/sqrt(FrameL)*fft(txSimbols,[],2);
%txSimbols=reshape(txSimbols,NDataInd,FrameL);

frecSimb=zeros(NFFT,FrameL);
frecSimb(dataIndices,:)=txSimbols;
frecSimb(pilotsIndices,:)=repmat(valoresPiloto,1,FrameL);

tsimbol=sqrt(NFFT)*ifft(fftshift(frecSimb,1),[],1); %IFFT normalizada
tsimbol=[zeros(NFFT,NullSimb) preambTime tsimbol];
txSignal=tsimbol([end-CP+1:end 1:end],:);  %Prefijo cíclico
txSignal=txSignal(:);
lMaxpw=(max(abs(txSignal))).^2;
if(lMaxpw>peakPw)
    peakPw=lMaxpw;
end

%% canal
txSignal=[txSignal;zeros(MaxDelay,1)]; %retardo aleatorio para prueba de sincronía

txSignalDelay=delayseq(txSignal,delayrnd);
%t=(0:numel(txSignal)-1);
%rxSignal=filter(hn,1,txSignalDelay);
rxSignal=chan(txSignalDelay);
rxSignal=rxSignal+sqrt(Pn/2)*(randn(Nmuestras*(FrameL+...
    NullSimb+1)+MaxDelay,1)+1j*randn(Nmuestras*(FrameL+NullSimb+1)+MaxDelay,1));
dW=2*pi*(2*rand-1)/(1.5*NFFT);
rxSignal=rxSignal.*transpose(exp(1j*dW*(0:numel(rxSignal)-1)));
%% Receptor

% Sincronia
%rxSignal=agc1(rxSignal);
sRxSig=delayseq(rxSignal,NFFT/2);
acorrM=conj(rxSignal).*sRxSig;
pw=conj(rxSignal).*rxSignal;
movavgWindow = dsp.MovingAverage(NFFT/2);
pD=movavgWindow(acorrM);
rD=movavgWindow(pw);
mD=abs(pD).^2./(rD.^2+.03);

detThr=mD>0.6;
if(sum(detThr))
    syncThr=abs(mD-max(mD))<(0.04*max(mD));
    syncWLen=sum(syncThr);
    meanInd=floor(sum(syncThr.*(1:numel(syncThr))')/syncWLen);
    %meanInd=min(max(floor(sum(syncThr.*(1:numel(syncThr))')/syncWLen),3*Nmuestras),numel(rxSignal)-Nmuestras*(FrameL+NullSimb+1)+3*Nmuestras-8);
    estSyncInd=meanInd-3*Nmuestras+CP/4-6;
end


cfOest=2*angle(pD(meanInd))/NFFT;
%sincInd=idealSincDelay; %usando valor real de retardo de Rx
sincInd=estSyncInd; %Usar indice estimado

rxSignal=rxSignal.*transpose(exp(1j*cfOest*(0:numel(rxSignal)-1))); %corrección de CFO

%Calculo de metrica de sincronía
rxFrame=reshape(rxSignal((1:Nmuestras*(FrameL+NullSimb+1))+sincInd),Nmuestras,FrameL+NullSimb+1); %%Reorganización de trama
rxFrameData=rxFrame(CP+1:end,2+NullSimb:end); %quitar CP
Rxsimbol=1/sqrt(NFFT)*fftshift(fft(rxFrameData,[],1),1); %%FFT normalizada
RxPilotsimbol=Rxsimbol(pilotsIndices,:); 
RxDatasimbol=Rxsimbol(dataIndices,:);

%Estimación del canal 
canalEnPiloto = RxPilotsimbol./repmat(valoresPiloto,1,FrameL);

%Interpolación BEM
%CanalenDatos=bemIntMtx*canalEnPiloto;
%Canal=bemFrec*((bemFrecPilots'*bemFrecPilots)\bemFrecPilots')*canalEnPiloto;
CanalenDatos=bem2DIntMtxData*canalEnPiloto(:);
CanalenDatos=reshape(CanalenDatos,NDataInd,FrameL);
tvCIR=bem2DIntMtxDelayTime*canalEnPiloto(:);
tvCIR=reshape(tvCIR,Lmax,FrameL);
pdpEst=mean(conj(tvCIR).*tvCIR,2);
deltaFcorr=1/sqrt(NFFT)*fftshift(fft(pdpEst,NFFT,1));

CanalFull=bem2DIntMtxFull*canalEnPiloto(:);
CanalFull=reshape(CanalFull,NFFT,FrameL);

%Ecualización
Rxeqsimbol=RxDatasimbol.*conj(CanalenDatos)./(abs(CanalenDatos).^2+Pn); %Ecualizador de un tap
%Rxeqsimbol=sqrt(NDataInd)*sqrt(FrameL)*ifft2(Rxeqsimbol);
Rxeqsimbol=sqrt(FrameL)*ifft(Rxeqsimbol,[],2);
Rxeqsimbol=sqrt(NDataInd)*ifft(Rxeqsimbol,[],1);
%Rxeqsimbol=reshape(Rxeqsimbol,NDataInd,FrameL);
rxBits=qamdemod(Rxeqsimbol,M,'OutputType','bit','UnitAveragePower',true);

%% gráficas


 set(hline, 'YData', imag(RxDatasimbol(:)));  % Update the y data of the line
 set(hline, 'XData', real(RxDatasimbol(:)));  % Update the y data of the line
 set(hline2, 'YData', imag(Rxeqsimbol(:)));  % Update the y data of the line
 set(hline2, 'XData', real(Rxeqsimbol(:)));  % Update the y data of the line
sf1.ZData=1/sqrt(NFFT)*(abs(fftshift((fft(tvCIR,NFFT,2)),2)));
sf2.CData=(abs(CanalFull));

 % figure(1)
% scatter(real(RxDatasimbol(:)),imag(RxDatasimbol(:)),'.c')
% xlim([-2.5 2.5]);ylim([-2.5 2.5])
% hold on
% scatter(real(Rxeqsimbol(:)),imag(Rxeqsimbol(:)),'.b')
% xlim([-2.5 2.5]);ylim([-2.5 2.5])
if(sincInd<=idealSincDelay)
    k=k+1;
    terros=terros+sum(sum(abs(rxBits-txBits)));
    ber(k)=terros/(k*numel(rxBits));
end
f2p1.XData=(1:numel(rxSignal));
f2p1.YData=(0.2*abs(rxSignal));
f2p2.XData=(1:numel(syncThr));
f2p2.YData=(syncThr);
f2p3.XData=(1:numel(mD));
f2p3.YData=(mD);
f2s1.XData=estSyncInd+2*Nmuestras;
f2s1.YData=0.2*abs(rxSignal(estSyncInd+2*Nmuestras));
f2s2.XData=idealSincDelay+2*Nmuestras;
f2s2.YData=0.2*abs(rxSignal(idealSincDelay+2*Nmuestras));
f2s3.XData=meanInd;
f2s3.YData=1;
f2.CurrentAxes.XLim=[meanInd-Nmuestras, meanInd+NFFT];

f5p1.XData=1:i;
f5p1.YData=ber(1:i);

drawnow
Tsimbolo*(FrameL+3)*(i-1)

% if(idealSincDelay<estSyncInd)
% input('press key to continue','s')
% end
end

% figure(3)
% surf((fftshift(abs(fft(tvCIR,[],2)),2))); shading interp

% 
% figure(4)
% subplot(2,1,1)
% plot((0:CP-1)*Ts,pdpEst(1:CP),'.-')
% subplot(2,1,2)
% plot(-Fs/2:Fs/NFFT:Fs/2-Fs/NFFT,10*log10(abs(deltaFcorr)))