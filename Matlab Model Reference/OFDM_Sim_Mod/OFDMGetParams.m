function [OFDMParams] = OFDMGetParams(OFDMConfig)
OFDMParams=OFDMConfig;

if(OFDMConfig.VTmode==1)
    OFDMParams.RUStrmap=reshape(upsample(OFDMParams.RUStrmap(:),2),2*size(OFDMParams.RUStrmap,1),size(OFDMParams.RUStrmap,2));
end


OFDMParams.RUStrsize=size(OFDMParams.RUStrmap);
OFDMParams.NRU=2*floor((OFDMParams.N-1)*OFDMParams.EffBW/(2*OFDMParams.RUStrsize(1)));
OFDMParams.M=OFDMParams.MRU*OFDMParams.RUStrsize(2); %Total number of OFDM Symbols
OFDMParams.Nu=OFDMParams.RUStrsize(1)*OFDMParams.NRU; %Total number of subcarriers in the RU
if(OFDMParams.Nu>OFDMParams.N-2) 
       error("Frame dimensions inconsistent n-subcarriers greater than NFFT")
end
OFDMParams.Ng=OFDMParams.N-OFDMParams.Nu; %Guard subcarriers

%Generate frame with tiled RUs
FrameMap=[ ...
    zeros(OFDMParams.Ng/2,OFDMParams.RUStrsize(2)); ...
    repmat(OFDMParams.RUStrmap,OFDMParams.NRU/2,1); ...
    zeros(1,OFDMParams.RUStrsize(2));...
    repmat(OFDMParams.RUStrmap,OFDMParams.NRU/2,1); ...
    zeros(OFDMParams.Ng/2-1,OFDMParams.RUStrsize(2))];

%%Per symbol Indexes%%%%%%%%%%%%%
% DataIndexesSymbol=cell(RUStr.size(2),1);
% for i=1:RUStr.size(2)
% DataIndexesSymbol{i}=
% 
% end
%%%%%%%%%%%%%%%%%%%%%%


if(OFDMParams.VTmode==1)
    FrameMap(2:OFDMParams.N/2+1,:)=FrameMap(1:OFDMParams.N/2,:);
end

%show Frame mapping Graph
%figure(1); imagesc(FrameMap)

%% RU-FrameBlock(7) mapping 
OFDMParams.Frameblock = numel(FrameMap(1,:));
OFDMParams.RUDataFBIndVec=uint32(find(FrameMap(:)==1));
OFDMParams.RUSizeDataFB=numel(OFDMParams.RUDataFBIndVec); %Total number of Data subcarriers in FrameBlock
OFDMParams.RUSizeDataBitsFB=OFDMParams.RUSizeDataFB*log2(OFDMParams.Mary); %
%%Pilots
OFDMParams.RUPilotFBIndexesVec=uint32(find(FrameMap(:)==2));
OFDMParams.RUSizePilotFB=numel(OFDMParams.RUPilotFBIndexesVec); %Total number of pilots in RU FrameBlock



%% Complete Frame
FrameMap=repmat(FrameMap,1,OFDMParams.MRU);
OFDMParams.FrameMapVec=FrameMap(:); %Linear Frame Map

%% RU-Wide carrier mapping
%%Data
OFDMParams.RUDataIndexesVec=uint32(find(OFDMParams.RUStrmap(:)==1));
OFDMParams.RUSizeData=numel(OFDMParams.RUDataIndexesVec); %Total number of Data subcarriers
OFDMParams.RUSizeDataBits=OFDMParams.RUSizeData*log2(OFDMParams.Mary);
%%Pilots
OFDMParams.RUPilotIndexesVec=uint32(find(OFDMParams.RUStrmap(:)==2));
OFDMParams.RUSizePilot=numel(OFDMParams.RUPilotIndexesVec); %Total number of pilots
%% Frame-Wide carrier mapping
%%Data
OFDMParams.DataIndexesVec=uint32(find(OFDMParams.FrameMapVec==1));
OFDMParams.SizeData=numel(OFDMParams.DataIndexesVec); %Total number of Data subcarriers
OFDMParams.SizeDataBits=OFDMParams.SizeData*log2(OFDMParams.Mary);
%%Pilots
OFDMParams.PilotIndexesVec=uint32(find(OFDMParams.FrameMapVec==2));
OFDMParams.SizePilot=numel(OFDMParams.PilotIndexesVec); %Total number of pilots
%Random Pilots 
OFDMParams.TxPilots=sqrt(OFDMParams.PilotPw)*(2*randi(2,OFDMParams.SizePilot,1)-3);
%%Hadamrd Precoding
%HdLengths=transpose(2.^(1:14))*([12 20]);
%HdLengths=sort(HdLengths(:));
%HdLengths=HdLengths(HdLengths<OFDMParams.SizeData);
%OFDMParams.HmL=HdLengths(end);
%OFDMParams.Hm=hadamard(HdLengths(end));
OFDMParams.HmL=2^(floor(log2(OFDMParams.SizeData)));
temp=OFDMParams.TxPilots;
%% Symbol-Wide mapping
PerSymbol=cell(OFDMParams.RUStrsize(2),1);
PerRU=cell(OFDMParams.RUStrsize(2),1);
for i=1:OFDMParams.RUStrsize(2)
    PerSymbol{i}=struct;
    PerSymbol{i}.DataInd=uint32(find(FrameMap(:,i)==1));
    PerSymbol{i}.nDataSym=numel(PerSymbol{i}.DataInd);
    PerSymbol{i}.nDataBits= log2(OFDMParams.Mary)*PerSymbol{i}.nDataSym;
    PerSymbol{i}.PilotInd=uint32(find(FrameMap(:,i)==2));
    PerSymbol{i}.EmptyInd=uint32(find(FrameMap(:,i)==0));
    PerSymbol{i}.Pilot=temp(1:numel(PerSymbol{i}.PilotInd));
    PerSymbol{i}.hasPilots=numel(PerSymbol{i}.Pilot)>0;
    temp(1:numel(PerSymbol{i}.PilotInd))=[];

    PerRU{i}=struct;

    PerRU{i}.hasPilots=numel(PerSymbol{i}.Pilot)>0;
end
OFDMParams.PerSymbol=PerSymbol;
%%Non-used
OFDMParams.EmptIndexesVec=uint32(find(OFDMParams.FrameMapVec==0));
OFDMParams.SizeEmpt=numel(OFDMParams.EmptIndexesVec); %Total number of unused subcarriers

%% Preamble
if OFDMParams.WIFI_PREAM
    short = generatePreambleSTF();
    long = generatePreambleLTF();
    full = [short long];
    OFDMParams.preambTime = sqrt(OFDMParams.PrePW) * full/rms(full);
    OFDMParams.preambTime = transpose(OFDMParams.preambTime);
else
    pnSeq = comm.PNSequence('Polynomial', [10 7 0], ...
    'SamplesPerFrame', OFDMParams.N/2, 'InitialConditions', [1 0 0 1 0 1 0 0 0 1]);
    %preamFrec = 1/sqrt(NFFT/2)*fft(2*pnseq()-1);
    preamFrecLT=2*pnSeq()-1;
    preamFrecLT=upsample(preamFrecLT,2);
    preamFrecLT=flipud(preamFrecLT);
    OFDMParams.preambTime=sqrt(OFDMParams.PrePW*2)*sqrt(OFDMParams.N)...
        *ifft(fftshift(preamFrecLT));
end

OFDMParams.TxPW=(OFDMParams.SizeData+OFDMParams.PilotPw...
    *OFDMParams.SizePilot)/(OFDMParams.N*OFDMParams.M);
OFDMParams.FrameL=(OFDMParams.N+OFDMParams.CP)*(OFDMParams.M+1);

end

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
    %s64   = s64 / rms(s64);                % potencia unitaria (opcional)
    preamShort = repmat(s64(1:16), 1, 10);      % 160 muestras = 8 µs
end