function [TxBits,TxDataSymbols] = OFDMGenPayLoad(OFDMParams,varargin)

if nargin>1
    TxBits=varargin{1};
else
   TxBits=uint8(randi(2,OFDMParams.SizeData*log2(OFDMParams.Mary),1)-1);
end

TxDataSymbols=qammod(TxBits,OFDMParams.Mary,'InputType','bit','UnitAveragePower',true);

%TxDataSymbols(:)=1;