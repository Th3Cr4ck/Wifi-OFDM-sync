function [RxBits,RxSym] = OfdmZfoneTap (OFDMParams,RxFrecFrame,CIR,NoisePW,G)
RxSymbols=zeros(OFDMParams.SizeData,1);
%RxBits=zeros(OFDMParams.SizeDataBits,1);
sIndex=1;
for i=1:OFDMParams.M
    g = G{i,1};
    RUnumb=mod(i-1,OFDMParams.RUStrsize(2))+1;
    y=RxFrecFrame(:,i);
    %%ZF
    g = diag(g);
    g=g(OFDMParams.PerSymbol{RUnumb}.DataInd);
    y=y(OFDMParams.PerSymbol{RUnumb}.DataInd);
    x=conj(g).*y./((abs(g)).^2+NoisePW); %ZF
    if OFDMParams.DFTs==1
        x=sqrt(numel(x))*ifft(x); %DFTs
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
   % RxSymbols(1:OFDMParams.HmL)=1/sqrt(OFDMParams.HmL)*...
   % ifwht(RxSymbols(1:OFDMParams.HmL),OFDMParams.HmL);
end
RxBits=uint8(qamdemod(RxSymbols,OFDMParams.Mary,'OutputType','bit','UnitAveragePower',true));
if nargout>1
    RxSym=RxSymbols;
end
