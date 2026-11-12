function [TxSignal] = OFDMGenTxFrame(OFDMParams,TxDataSymbols)

if OFDMParams.DFTsFrame==1
     TxDataSymbols=1/sqrt(numel(TxDataSymbols))*fft(TxDataSymbols);
   %  TxDataSymbols = matintrlv(TxDataSymbols,OFDMParams.RUSizeData,...
    %     OFDMParams.MRU*OFDMParams.NRU);
     
end

TxFrameFrec=zeros(OFDMParams.N,OFDMParams.M);
TxFrameFrec(OFDMParams.DataIndexesVec)=TxDataSymbols;

%Pïlots
TxFrameFrec(OFDMParams.PilotIndexesVec)=OFDMParams.TxPilots;

%%DFTs
if OFDMParams.DFTs==1
    for i=1:OFDMParams.M
        RUi=mod(i-1,OFDMParams.RUStrsize(2))+1;
        index=OFDMParams.PerSymbol{RUi}.DataInd;
        TxFrameFrec(index,i)=1/sqrt(numel(index))*fft(TxFrameFrec(index,i));
    end
end

%Transform to Time domain
TxFrameTime=sqrt(OFDMParams.N)*ifft(fftshift(TxFrameFrec,1),[],1);

% Preambulo
nSTFSamples = 16*10;
nLTFSamples = 64*2;

if OFDMParams.WIFI_PREAM
    TxFrameTime = [...
        OFDMParams.preambTime(1:nSTFSamples, :)...
        OFDMParams.preambTime(nSTFSamples+1:nSTFSamples+1+nLTFSamples, :)...
        TxFrameTime ];
else
    TxFrameTime = [OFDMParams.preambTime TxFrameTime];
end

%CP
GI2 = 32; % STF guard
GI  = 16; % OFDM Symbol guard

payloadStartIdx = nSTFSamples+nLTFSamples+1;

if OFDMParams.WIFI_PREAM
    TxFrameTime = [...
        TxFrameTime(1:payloadStartIdx-1, :)...
        TxFrameTime(end-GI+1:end, :)...
        TxFrameTime(payloadStartIdx:end, :) ];

else
    TxFrameTime=TxFrameTime([end-OFDMParams.CP+1:end 1:end],:);
end
TxSignal=TxFrameTime(:);
