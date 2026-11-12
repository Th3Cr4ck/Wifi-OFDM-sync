function [RxBits,RxSym] = OfdmLmmseEq (OFDMParams,RxFrecFrame,CIR,NoisePW,G2)

RxSymbols=zeros(OFDMParams.SizeData,1);
%RxBits=zeros(OFDMParams.SizeDataBits,1);
sIndex=1;

subRu=OFDMParams.NRU/2;
for i=1:OFDMParams.M
    G = G2{i,1};
    RUnumb=mod(i-1,OFDMParams.RUStrsize(2))+1;
    y=RxFrecFrame(:,i);
    %%%%LMMSE
    if OFDMParams.PerSymbol{RUnumb}.hasPilots
        tempPilotVec=zeros(OFDMParams.N,1);
        tempPilotVec(OFDMParams.PerSymbol{RUnumb}.PilotInd)=...
                                OFDMParams.PerSymbol{RUnumb}.Pilot;
        y=y-G*tempPilotVec;
    end
%     for k=1:subRu
%         GTemp=G((k-1)*OFDMParams.RUStrsize(1)+1:k*OFDMParams.RUStrsize(1),...
%             (k-1)*OFDMParams.RUStrsize(1)+1:k*OFDMParams.RUStrsize(1));
% 
%         GTemp=GTemp(:,OFDMParams.PerSymbol{RUnumb}.DataInd);
%         x=(G'*G+NoisePW.*eye(size(G,2)))\G'*y; 
%     end
    
     G=G(:,OFDMParams.PerSymbol{RUnumb}.DataInd);
     x=(G'*G+NoisePW.*eye(size(G,2)))\G'*y; 
    
    if OFDMParams.DFTs==1
        x=sqrt(numel(x))*ifft(x);%DFTs
    end
     RxSymbols(sIndex:sIndex+OFDMParams.PerSymbol{RUnumb}.nDataSym-1)=x;
     sIndex=sIndex+OFDMParams.PerSymbol{RUnumb}.nDataSym;

end
% RxSymbolsArr=transpose(reshape(RxSymbols,OFDMParams.RUSizeData,...
%     OFDMParams.NRU*OFDMParams.MRU));
% RxSymbolsArr=sqrt(OFDMParams.NRU*OFDMParams.MRU)*ifft(RxSymbolsArr,[],1);
% RxSymbols=RxSymbolsArr(:);
if OFDMParams.DFTsFrame==1
    RxSymbols=sqrt(numel(RxSymbols))*ifft(RxSymbols);
     
    %RxSymbols(1:OFDMParams.HmL)=1/sqrt(OFDMParams.HmL)*...
   % ifwht(RxSymbols(1:OFDMParams.HmL),OFDMParams.HmL);
end
    RxBits=qamdemod(RxSymbols,OFDMParams.Mary,'OutputType','bit','UnitAveragePower',true);
if nargout>1
    RxSym=RxSymbols;
end
